#! /usr/bin/env python
#-*- coding: utf-8 -*-
 
# =================================================================
#   -------------------------------------------------------------
#   Context : Extract transcript exonic sizes from a gff or a gtf
#             By : BIZET Martin
#                 Year : 2017
#   -------------------------------------------------------------
# =================================================================

##### ========================
####             
###    MAIN 
##               
# ============================

import sys

ofile = ''
ifile = ''
fID = 'transcript_id'

#Parse arguments
args = sys.argv
for i in range(1, len(args)):
	if (args[i] == '-h') or (args[i] == '-help'):
		print 'Extract transcript size from a gff or a gtf file\n\t-i or -input: gff/gtf file to read\n\t-o or -out: output file\n\t-f or -featureID: feature to extract the size from [default="transcript_id"]\n\t-h or -help: this help' 
		argTyp = ''
		argVal = ''
	elif '=' in args[i]:
		el = args[i].split('=')
		argTyp = el[0]
		argVal = el[1]
	elif i < (len(args)-1):
		argTyp = args[i]
		argVal = args[i+1]
	if (argTyp == '-i') or (argTyp == '--input'):
		ifile = argVal
	elif (argTyp == '-o') or (argTyp == '--out'):
		ofile = argVal
	elif (argTyp == '-f') or (argTyp == '--featureID'):
		fID = argVal
	
print 'input: ' + ifile
print 'output: ' + ofile
print 'featureID: ' + fID

if (ifile != '') and (ofile != ''):
	#Read input
	trExons = {}
	ifile = open(ifile)
	for line in ifile:
		if line[0] != '#':
			line = line.split('\t')
			if line[2] == 'exon':
				sta = int(line[3])
				en = int(line[4])
				tr = line[8].split(fID)[1].split(';')[0]
				tr = tr.replace('=', '')
				tr = tr.replace(' ', '')
				tr = tr.replace('\"', '')
				if tr in trExons:
					exID = len(trExons[tr])/2
					trExons[tr] += [(sta, exID, 0, 'Start'), (en, exID, 1, 'End')]
				else:
					trExons[tr] = [(sta, 0, 0, 'Start'), (en , 0, 1, 'End')]
	ifile.close()
	#compute size
	toWrite = '\t'.join([fID, 'size']) + '\n'
	for key in sorted(trExons.keys()):
		tmp = trExons[key]
		tmp.sort()
		size = 0
		opened = 0
		currentSta = 0
		for ex in tmp:
			if ex[3] == 'Start':
				if opened == 0:
					currentSta = ex[0]
				opened += 1
			else:#if ex[3] == 'End'
				opened -= 1
				if opened == 0:
					size += ex[0] - currentSta + 1		
		toWrite += key + '\t' + str(size) + '\n'
	#Write out
	ofile = open(ofile, 'w')
	ofile.write(toWrite)
	ofile.close()
