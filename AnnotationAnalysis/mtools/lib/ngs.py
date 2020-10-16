# ------------------------------------------------------------
#
# Matthieu Defrance ULB 2016
# lib/ngs
#
# ------------------------------------------------------------

import os
import sys
import random
import time
from base import *

# ------------------------------------------------------------
# STAR
# ------------------------------------------------------------
# 1. index
# copy mm9.fa to the directory; cd to starIndex directory (e.g mm9/star) 
# STAR  --runMode genomeGenerate --runThreadN 12 --genomeDir ./ --genomeFastaFiles mm9.fa
# 2. mapping
# STAR --genomeDir mm9-starIndex/ --runThreadN 12 --readFilesIn read1.fastq read2.fastq 
# --outFileNamePrefix Experiment1Star
# input reads in fastq.gz
# output sorted bam
def star(reads, bam, genomeDir, threads = 8):
    tmp = tmp_dir()
    run('STAR',
    '--genomeDir %s'         % genomeDir,
    '--runThreadN %d'        % threads,
    '--readFilesIn %s'       % reads,
    '--outFileNamePrefix %s' % tmp,
    '--readFilesCommand zcat',
    '--outSAMtype BAM SortedByCoordinate')
    run('cp %s/Aligned.sortedByCoord.out.bam %s' % (tmp, bam))
    clean_tmp(tmp)

# ------------------------------------------------------------
# BWA (DNA)
# ------------------------------------------------------------
# 1. index
# copy mm9.fa to the directory; cd to bwa_index directory
# 2. mapping
# "bwa index -a bwtsw mm9.fa

# def bwa(reads, sam, genome, threads = 8):
#     fa = sprintf('/illumina/runs/Data/genomes/%s/%s.fa', genome, genome)
#     tmp = tmp_file()
#     if pairedEnd:
#         action = 'sampe'
#     else:
#         action = 'samse'
#     run('bwa aln -t %d %s %s > %s' % (threads, genome, reads, tmp))
#     run('bwa %s %s %s %s > %s' % (action, fa, tmp, reads, sam))
#     clean_tmp(tmp)

# ------------------------------------------------------------
# macs2
# ------------------------------------------------------------
def macs2(treatment, control, genome, xls, bdg, narrowPeak, summits):
    genome_size = {
        'hg19'    : 'hs',
        'mm9'     : 'mm',
        'mm9-tr'  : '61400000',
        'hg19-tr' : '70100000'
    }

    tmp   = tmp_dir()
    run('macs2 callpeak',
        '-t %s'       % treatment,
        '%s'          % select(control == 'NA', '', '-c %s' % control),
        '-g %s'       % genome_size[genome],
        '--outdir=%s' % tmp,
        '-f AUTO',
        '-n MACS -B -q 0.05',
        '--call-summits')
    run('cp %s/MACS_peaks.xls %s'        % (tmp, xls))
    run('cp %s/MACS_treat_pileup.bdg %s' % (tmp, bdg))
    run('cp %s/MACS_peaks.narrowPeak %s' % (tmp, narrowPeak))
    run('cp %s/MACS_summits.bed      %s' % (tmp, summits))

    clean_tmp(tmp)

def  bdg_to_wig(bdg, wig, step = 40):
    title = wig.split('/')[-1].split('.')[0]
    run('bedgraph_to_wig.pl --bedgraph %s --wig %s --step %d --name \"%s\"' % (bdg, wig, step, title))
    run('gzip -f %s' % wig)

# ------------------------------------------------------------
# read / write peaks
# ------------------------------------------------------------
def read_peaks_xls(filename):
    f = open(filename)
    data = []
    for line in f:
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        e = line.split('\t')
        if e[0] == 'chr':
            continue
        e[1] = int(e[1])
        e[2] = int(e[2])
        e[3] = int(e[3])
        e[4] = int(e[4])
        e[5] = float(e[5])
        e[6] = float(e[6])
        e[7] = float(e[7])

        if len(e) == 8:
            e = e + [100]
        else:
            e[8] = float(e[8])
        
        data += [e]
    return data

def read_peaks_bed(filename):
    f = open(filename)       
    data = []                
    for line in f:           
        line = line.strip() # order: chr, start, stop, name, score, strand, ...
        if line.startswith('#') or line == '' or line.startswith('track name'):
            continue         
        e = line.split('\t')
        new_e = [] # order: chr, start, stop, length peaks, summits, and rest except name is in index 9
        new_e.append( str(e[0]) ) # 0: chr 
        new_e.append( int(e[1]) ) # 1: start    
        new_e.append( int(e[2]) ) # 2: stop    
        new_e.append( new_e[2] - new_e[1] ) # 3: length
        new_e.append( new_e[1] + int( float(new_e[-1]) / 2 ) ) # 4: summit
        new_e = new_e + e[3:] # 5: name, 6: score, 7: strand, ...
        
        data += [new_e]          
    return data

def write_peaks_xls(peaks, filename):
    f = open(filename, 'w')
    header = ['chr', 'start', 'end', 'length', 'abs_summit', 
            'pileup', '-log10(pvalue)', 'fold_enrichment', '-log10(qvalue)', 'name', 'annotation']
    f.write('\t'.join(header))
    f.write('\n')
    for p in peaks:
        f.write('\t'.join([str(i) for i in p]))
        f.write('\n')
    f.close()

def write_peaks_bed( peaks, filename, header_sup ):
    f = open(filename, 'w')
    header = ['chr', 'start', 'end', 'length', 'abs_summit'] + header_sup + ['annotation']
    f.write('\t'.join(header))
    f.write('\n')
    for p in peaks:
        f.write('\t'.join([str(i) for i in p]))
        f.write('\n')
    f.close()

# ------------------------------------------------------------
# RIP peaks annotation
# ------------------------------------------------------------
# read refFlat UCSC annotation
def read_genomic_annotation(filename, detectFormat=True):

    g = {}
    h = {}
    f = open(filename)

    name2_Pos = 0
    exon_starts_Pos = -2
    exon_ends_Pos = -1
    if detectFormat:
        print('Start RefFlat format detection:...')
        e = next(f)
        l = e.strip()
        if l.startswith('#'):
            l = l.split()
            print('Header:')
            print(l)
            if l[0] == '#bin' and 'name2' in l:
                while l[name2_Pos] != 'name2':
                    name2_Pos += 1
                print('Gene Name at column ' + str(name2_Pos))
            else:
				print('Gene Name not found, using legacy position: ' + str(name2_Pos))
            i = 0; foundSt = False; foundEnd = False;
            while (not foundSt) or (not foundEnd):
                if l[i] == 'exonStarts':
                    foundSt = True; exon_starts_Pos = i;
                    print('exonStarts at column ' + str(exon_starts_Pos))
                elif l[i] == 'exonEnds':
                    foundEnd = True; exon_ends_Pos = i;
                    print('exonEnds at column ' + str(exon_ends_Pos))
                i += 1
            if not foundSt:
				print('exonStarts not found, using legacy position: ' + str(exon_starts_Pos))
            if not foundEnd:
				print('exonEnds not found, using legacy position: ' + str(exon_ends_Pos))
        else:
            print('Header not found, using legacy format: geneName at first column and exons positions at the 2 lasts columns')
    else:
        print('legacy format: geneName at first column and exons positions at the 2 lasts columns')

    for e in f:
        l = e.strip()
        if l.startswith('#') or l == '':
            continue
        l                       = l.split()
        geneId, id, chr, strand = l[name2_Pos], l[1], l[2], l[3] #legacy: l[0], l[1], l[2], l[3]
        tr_a, tr_b              = int(l[4]), int(l[5])
        cds_a, cds_b            = int(l[6]), int(l[7])
        exon_starts             = [int(p) for p in l[exon_starts_Pos].split(',')[:-1]] #legacy: [int(p) for p in l[-2].split(',')[:-1]]
        exon_ends               = [int(p) for p in l[exon_ends_Pos].split(',')[:-1]] #legacy: l[-1].split(',')[:-1]]

        entries = h.get(chr, [])
        entries += [[chr, tr_a, tr_b, id, geneId, strand, cds_a, cds_b, exon_starts, exon_ends] ]
        h[chr] = entries

    f.close()

    return h

def annotate_peaks_RNA(data, annotation, useStr=False, promo='None', nearTTS='None'):
    if promo != 'None': #give promoter borders in bp comma separated (ex: 2000,1000 for -2kb to +1kb around TSS)
		promo = promo.split(',')
		promo[0] = -1*int(promo[0])
		promo[1] = int(promo[1])
    if nearTTS != 'None': #give nearTTS region borders in bp comma separated (ex: 1000,1000 for -1kb to +1kb around TTS)
		nearTTS = nearTTS.split(',')
		nearTTS[0] = -1*int(nearTTS[0])
		nearTTS[1] = int(nearTTS[1])
    for e in data:
        chr, i, j, k = e[0], e[1], e[2], e[4]
        if useStr:
			if e[7] == "+" or e[7] == "-":
				e_str = e[7] #assume it is from read_peaks_bed (read_peaks_xls do not have strand information)
			else:
				raise ValueError('useStr was True but strand value was'+e[7]+', while expecting + or - .')
			
        e_anno = []
        for t in annotation.get(chr, []):
            chr, tr_a, tr_b, id, geneId, strand = t[0], t[1], t[2], t[3], t[4], t[5]
            cds_a, cds_b = t[6], t[7]
            if promo != 'None':
				if strand == '+':
					p_a = tr_a + promo[0]
					p_b = tr_a + promo[1]
				else:
					p_a = tr_b - promo[1]
					p_b = tr_b - promo[0]
            else:
				p_a = float('Inf')
				p_b = float('-Inf')				
            if nearTTS != 'None':
				if strand == '+':
					nt_a = tr_b + nearTTS[0]
					nt_b = tr_b + nearTTS[1]
				else:
					nt_a = tr_a - nearTTS[1]
					nt_b = tr_a - nearTTS[0]
            else:
				nt_a = float('Inf')
				nt_b = float('-Inf')				
            if useStr and not (k >= min(p_a, tr_a, nt_a) and k <= max(p_b, tr_b, nt_b) and e_str == strand):
                continue
            if not (k >= min(p_a, tr_a, nt_a) and k <= max(p_b, tr_b, nt_b) ):
                continue
            # initialise exon_anno with the correct noncoding/coding status
            if cds_a == cds_b:
                exon_anno = 'ncIntron'
            else:
                exon_anno = 'intron'

            # position in transcript
            if strand == '+':
                relative_position     = (k - tr_a) / float(tr_b - tr_a + 1)
                distance_to_TSS       = k - tr_a
                distance_to_TES       = k - tr_b
                if exon_anno == 'ncIntron': #noncoding
                    distance_to_CDS_start = 0
                    distance_to_CDS_end   = 0
                else:
                    distance_to_CDS_start = k - cds_a
                    distance_to_CDS_end   = k - cds_b

            else:
                relative_position     = (tr_b - k) / float(tr_b - tr_a + 1)
                distance_to_TSS       = tr_b - k
                distance_to_TES       = tr_a - k
                if exon_anno == 'ncIntron': #noncoding
                    distance_to_CDS_start = 0
                    distance_to_CDS_end   = 0
                else:
                    distance_to_CDS_start = cds_b - k
                    distance_to_CDS_end   = cds_a - k
            # PROMOTER / near TTS
            if k >= p_a and k <= p_b: #never True when promo=None' has p_a is infinite (and p_b -infinite)
				promo_anno = 'Promoter'
            elif k >= nt_a and k <= nt_b: #never True when nearTTS=None' has nt_a is infinite (and nt_b -infinite)
				promo_anno = 'near_TTS'
            else:
				promo_anno = 'NonPromoter'
            # EXON
            exon_starts = t[-2]
            exon_ends   = t[-1]
            dex_tr = 0; pos_cds_a = 0; pos_cds_b = 0;
            dex_TSS = 0; dex_TTS = 0; dex_start = 0; dex_stop = 0;
            for i in range(len(exon_starts)):
                exon_size = exon_ends[i] - exon_starts[i] + 1
                
                if cds_a >= exon_starts[i] and cds_a <= exon_ends[i]:
                    pos_cds_a = dex_tr + cds_a - exon_starts[i] + 1
                
                if cds_b >= exon_starts[i] and cds_b <= exon_ends[i]:
                    pos_cds_b = dex_tr + cds_b - exon_starts[i] + 1
                
                if k >= exon_starts[i] and k <= exon_ends[i]:
                    if exon_anno == 'ncIntron': #noncoding
                        exon_anno = 'ncExon'
                    else:
                        exon_anno = 'exon'
                    pos_k = dex_tr + k - exon_starts[i] + 1
                
                dex_tr += exon_size

            # CDS
            if exon_anno == 'intron' or exon_anno == 'ncIntron':
                cds_anno = exon_anno
            else:
                if strand == '+':
                    if exon_anno == 'ncExon':
                        cds_anno = 'ncExon'
                    elif k < cds_a:
                        cds_anno = '5UTR'
                    elif k > cds_b:
                        cds_anno = '3UTR'
                    else:
                        cds_anno = 'CDS'
                    dex_TSS = pos_k - 1
                    dex_TTS = pos_k - dex_tr
                    if exon_anno != 'ncExon':
                        dex_start = pos_k - pos_cds_a
                        dex_stop  = pos_k - pos_cds_b
                else:
                    if exon_anno == 'ncExon':
                        cds_anno = 'ncExon'
                    elif k < cds_a:
                        cds_anno = '3UTR'
                    elif k > cds_b:
                        cds_anno = '5UTR'
                    else:
                        cds_anno = 'CDS'
                    dex_TSS = dex_tr - pos_k
                    dex_TTS = 1 - pos_k
                    if exon_anno != 'ncExon':
                        dex_start = pos_cds_b - pos_k
                        dex_stop  = pos_cds_a - pos_k
            if promo_anno == 'Promoter' and (cds_anno == 'ncExon' or cds_anno == 'ncIntron'):
				promo_anno = 'ncPromoter' # this info is needed for ncAware.get.feat in Annot_modifier.R
            elif promo_anno == 'near_TTS' and (cds_anno == 'ncExon' or cds_anno == 'ncIntron'):
				promo_anno = 'ncTTS' # this info is needed for ncAware.get.feat in Annot_modifier.R
            if promo_anno != 'NonPromoter': #(nc)Promoter and near_TTS have priority over other categories but I still to get dex_start and dex_stop info so I need to overwrite cds_anno at this step
                cds_anno = promo_anno

            e_anno += ['%s,%s,%s,%f,%d,%d,%d,%d,%d,%d,%d,%d' % 
                    (id, geneId, cds_anno, relative_position, 
                        distance_to_TSS, distance_to_TES, distance_to_CDS_start, distance_to_CDS_end,
                        dex_TSS, dex_TTS, dex_start, dex_stop)]

        if e_anno == []:
            e += ['NA']
        else:
            e += [';'.join(e_anno)]

    return data

def RNA_peaks_annotate(peaksFile, outputFile, annotationFile, promo='None', nearTTS='None'):
    useStr     = False #not available in the xls
    annotation = read_genomic_annotation(annotationFile)
    peaks      = read_peaks_xls(peaksFile)
    peaks      = annotate_peaks_RNA(peaks, annotation, useStr, promo, nearTTS)
    write_peaks_xls(peaks, outputFile)

def RNA_peaks_annotate_frombed(peaksFile, outputFile, annotationFile, header_sup, useStr=False, promo='None', nearTTS='None'):
    header_sup = header_sup.split(",")
    annotation = read_genomic_annotation(annotationFile)
    peaks      = read_peaks_bed(peaksFile)
    peaks      = annotate_peaks_RNA(peaks, annotation, useStr, promo, nearTTS)
    write_peaks_bed(peaks, outputFile, header_sup)


def extract_genes(data):
    genes = {}
    for e in data:
        if e[-1] == 'NA':
            continue
        ids = [i[1] for i in e[-1].split(';')]
        for id in ids:
            genes[id] = 1
    return genes.keys()

def extract_category(data):
    l = []
    for e in data:
        if e[-1] == 'NA':
            category = 'NA'
        else:
            features = [i.split(',')[2] for i in e[-1].split(';')]
            category = features[0]
            for f in features:
                if f != category:
                    category = 'multiple'
                    break
        l += [category]
    return l

def peaks_index(data):
    h = {}
    for e in data:
        chr = e[0]
        entries = h.get(chr, [])
        entries += [e]
        h[chr] = entries
    return h

def peaks_intersect(a, b, w = 1000, computeSummit = False):
    r = []
    h = peaks_index(b)
    for e in a:
        l = h.get(e[0], None)
        if l == None:
            continue
        for f in l:
            if computeSummit:
                se = (int(e[1]) + int(e[2])) / 2
                sf = (int(f[1]) + int(f[2])) / 2
                if abs(se - sf) < w:
                    r += [e]
                    break
            else:
                if abs(int(e[4]) - int(f[4])) < w:
                    r += [e]
                    break
    return r
