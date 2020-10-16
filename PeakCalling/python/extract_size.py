#! /usr/bin/env python
#-*- coding: utf-8 -*-
 
# =================================================================
#             ----------------------------------------------------
#            Context : Extract transcript sizes from a gff or a gtf
#                               By : BIZET Martin
#                                   Year : 2016
#            ----------------------------------------------------
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
	

if (ifile != '') and (ofile != ''):
	#Read input
	trSize = {}
	ifile = open(ifile)
	for line in ifile:
		line = line.split('\t')
		if line[2] == 'exon':
			sta = int(line[3])
			en = int(line[4])
			tr = line[8].split(fID)[1].split(';')[0]
			tr = tr.replace('=', '')
			tr = tr.replace(' ', '')
			tr = tr.replace('\"', '')
			if tr in trSize:
				trSize[tr][0] = min(trSize[tr][0], sta)
				trSize[tr][1] = max(trSize[tr][1], en)
			else:
				trSize[tr] = [sta, en]
	ifile.close()
	#compute size
	toWrite = '\t'.join([fID, 'size']) + '\n'
	for key in sorted(trSize.keys()):
		tmp = trSize[key]
		toWrite += key + '\t' + str(tmp[1]-tmp[0]+1) + '\n'
	#Write out
	ofile = open(ofile, 'w')
	ofile.write(toWrite)
	ofile.close()
