# ------------------------------------------------------------
#
# Bioinformatics Toolkit
# Matthieu Defrance ULB 
#
# ------------------------------------------------------------

# nohup python scripts/current/run.py &> log.2016.08.11.txt &

#
#SAMPLE_TABLE   = 'SEQ160817.csv'



import os
import sys
import random
import time

# ------------------------------------------------------------
# utils
# ------------------------------------------------------------
def select(cond, ifTrue, ifFalse):
    if cond:
        return ifTrue
    else:
        return ifFalse

# ------------------------------------------------------------
# log / run
# ------------------------------------------------------------
def log(msg):
    sys.stderr.write('[%s] %s\n' % (time.strftime('%Y/%m/%d %H:%M'), msg))
    sys.stderr.flush()

def run(*args):
    cmd = ' '.join(['%s' % str(v) for v in args])
    log('%s' % cmd)
    os.system(cmd)

# ------------------------------------------------------------
# csv read / write
# ------------------------------------------------------------
def read_csv(filename, sep = ',', header = True):
    f = open(filename)
    lines = f.readlines()
    lines = lines[1:]
    return [line.strip().split(sep) for line in lines]

# ------------------------------------------------------------
# tmp
# ------------------------------------------------------------
def tmp_dir(default_dir = '/illumina/runs/tmp'):
    tmp = '%s/%s%d/' % (default_dir, time.strftime('%y%m%d%H%M'), random.randint(1, 1000))
    run('mkdir -p %s' % tmp)
    return tmp

def tmp_file(default_dir = '/illumina/runs/tmp'):
    tmp = '%s/%s%d.tmp' % (default_dir, time.strftime('%y%m%d%H%M'), random.randint(1, 1000))
    return tmp

def clean_tmp(tmp):
    run('rm -rf %s' % tmp)

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
def star(reads, bam, genome, threads = 8):
    tmp = tmp_dir()
    run('STAR',
    '--genomeDir %s'         % '/illumina/runs/Data/genomes/%s/star' % genome,
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
# def macs2(treatment, control, genome, xls, bdg):
#     genome_size = {
#         'hg19'    : 'hs',
#         'mm9'     : 'mm',
#         'mm9-tr'  : '61400000',
#         'hg19-tr' : '70100000'
#     }

#     tmp   = tmp_dir()
#     run('macs2 callpeak',
#         '-t %s'       % treatment,
#         '-c %s'       % control,
#         '-g %s'       % genome_size[genome],
#         '--outdir=%s' % tmp,
#         '-f AUTO',
#         '-n MACS -B -q 0.05',
#         '--call-summits')
#     run('macs2 callpeak -t %s -c %s -f AUTO -g %s -n MACS -B -q 0.05 --outdir=%s --call-summits' %
#         (treatment, control, genome_size[genome], tmp, extra))
#     run('cp %s/MACS_peaks.xls %s'        % (tmp, xls))
#     run('cp %s/MACS_treat_pileup.bdg %s' % (tmp, bdg))
#     clean_tmp(tmp)
# run('cp %s/MACS_peaks.broadPeak %s' % (tmp, broadPeak))
# run('mkdir -p %s' % outdir)
# run('cp %s/MACS_peaks.xls %s/%s_peaks.xls'        % (tmp, outdir, name))
# run('cp %s/MACS_treat_pileup.bdg %s/%s.bdg' % (tmp, outdir, name))
# run('cp %s/MACS_peaks.narrowPeak %s/%s_peaks.narrowPeak' % (tmp, outdir, name))
# run('cp %s/MACS_summits.bed      %s/%s_summits.bed' % (tmp, outdir, name))
# run('cp %s/MACS_peaks.broadPeak %s/' % (tmp, outdir))

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
# peaks
# ------------------------------------------------------------
# def peaks_index(data):
#     h = {}
#     for e in data:
#         chr = e[0]
#         entries = h.get(chr, [])
#         entries += [e]
#         h[chr] = entries
#     return h

# def peaks_intersect(a, b, w = 1000, computeSummit = False):
#     r = []
#     h = bed_index(b)
#     for e in a:
#         l = h.get(e[0], None)
#         if l == None:
#             continue
#         for f in l:
#             if computeSummit:
#                 se = (int(e[1]) + int(e[2])) / 2
#                 sf = (int(f[1]) + int(f[2])) / 2
#                 if abs(se - sf) < w:
#                     r += [e]
#                     break
#             else:
#                 if abs(int(e[4]) - int(f[4])) < w:
#                     r += [e]
#                     break
#     return r

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

# ------------------------------------------------------------
# RIP peaks annotation
# ------------------------------------------------------------
# read refFlat UCSC annotation
def read_genomic_annotation(filename):
    g = {}
    h = {}
    f = open(filename)
    for e in f:
        l = e.strip()
        if l.startswith('#') or l == '':
            continue
        l                       = l.split()
        geneId, id, chr, strand = l[0], l[1], l[2], l[3]
        tr_a, tr_b              = int(l[4]), int(l[5])
        cds_a, cds_b            = int(l[6]), int(l[7])
        exon_starts             = [int(p) for p in l[-2].split(',')[:-1]]
        exon_ends               = [int(p) for p in l[-1].split(',')[:-1]]

        entries = h.get(chr, [])
        entries += [[chr, tr_a, tr_b, id, geneId, strand, cds_a, cds_b, exon_starts, exon_ends] ]
        h[chr] = entries

    return h

def annotate_peaks_RNA(data, annotation):
    for e in data:
        chr, i, j, k = e[0], e[1], e[2], e[4]
        e_anno = []
        for t in annotation.get(chr, []):
            chr, tr_a, tr_b, id, geneId, strand = t[0], t[1], t[2], t[3], t[4], t[5]
            cds_a, cds_b = t[6], t[7]

            if not (k >= tr_a and k <= tr_b):
                continue

            # position in transcript
            if strand == '+':
                relative_position = (k - tr_a) / float(tr_b - tr_a + 1)
                distance_to_TSS = k - tr_a
                distance_to_TES = k - tr_b
            else:
                relative_position = (tr_b - k) / float(tr_b - tr_a + 1)
                distance_to_TSS = tr_b - k
                distance_to_TES = tr_a - k

            # EXON
            exon_anno = 'intron'
            exon_starts, exon_ends = t[-2], t[-1]
            for i in range(len(exon_starts)):
                if k >= exon_starts[i] and k <= exon_ends[i]:
                    exon_anno = 'exon'

            # CDS
            if exon_anno == 'intron':
                cds_anno = 'intron'
            else:
                if strand == '+':
                    if k < cds_a:
                        cds_anno = '5UTR'
                    elif k > cds_b:
                        cds_anno = '3UTR'
                    else:
                        cds_anno = 'CDS'
                else:
                    if k < cds_a:
                        cds_anno = '3UTR'
                    elif k > cds_b:
                        cds_anno = '5UTR'
                    else:
                        cds_anno = 'CDS'

            e_anno += ['%s,%s,%s,%f,%d,%d' % 
                    (id, geneId, cds_anno, relative_position, distance_to_TSS, distance_to_TES)]

        if e_anno == []:
            e += ['NA']
        else:
            e += [';'.join(e_anno)]

    return data

def RNA_peaks_annotate(peaksFile, outputFile, annotationFile):
    annotation = read_genomic_annotation(annotationFile)
    peaks      = read_peaks_xls(peaksFile)
    peaks      = annotate_peaks_RNA(peaks, annotation)
    write_peaks_xls(peaks, outputFile)

def extract_genes(data):
    genes = {}
    for e in data:
        if e[-1] == 'NA':
            continue
        ids = [i[1] for i in e[-1].split(';')]
        for id in ids:
            genes[id] = 1
    return genes.keys()

# def extract_features(data, idx = 3):
#     genes = {}
#     for e in data:
#         features = [i[idx] for i in e[-1].split(';')]
#         for id in ids:
#             genes[id] = 1
#     return genes.keys()

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

