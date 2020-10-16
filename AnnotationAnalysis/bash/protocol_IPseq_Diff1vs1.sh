#!/usr/bin/env bash

###################################################################
# Description: Pipeline for differential RIP 1 sample vs 1 sample #
# Date: 2017/01                                                   #
# Author: Martin BIZET                                            #
###################################################################

### ======================
##    Default parameters
# ========================

conf=${COCKTAIL}/AnnotationAnalysis/Diff1vs1_config.sh
pathTsv=""
pathBam=""
d="Both"

### ================
##    ARGUMENTS
# ==================

HELP="
\nDESCRIPTION:
\n
\n\tTwo-conditions differential analysis with one sample by condition
\n
\n-c=|--config= config file. [default= Diff1vs1_config.sh]
\n-pF=|--peakFolder= folder containing peaks in tsv [if not provided absolute path should provided in -p1 and -p0 arguments
\n-p1=|--peaks1= name of peak 1(condition) in Tsv folder
\n-p0=|--peaks0= name of peak 0(control) in Tsv folder
\n-n=|--names=  comma separated names of IP1, IP0, INPUT1, INPUT0 in this order. If bL1_, bL0_, bI1_, bI0_ tags are not present they will be added (this step requires -n to not contain ° symbol) [default=-L1, -L0, -I1, -I0 filenames]
\n-o=|--out= ouput directory
\n-a=|--annotation= annotation file in refFlat format (absolute path)
\n-d=|--direction= direction of the differential analysis. It can be \'Both\', \'Hypo\' or \'Hyper\'. [default=Both]
\n-h=|--help= this help
"

# parse
for i in "$@"
do
case $i in
    -h*|--help*)
    echo -e ${HELP};
    exit 0
    ;;
    -c=*|--config=*)
    conf="${i#*=}"
    shift # past argument=value
    ;;
    -pF=*|--peakFolder=*)
    pathTsv="${i#*=}"
    shift # past argument=value
    ;;
    -p1=*|--peaks1=*)
    peaksf="${i#*=}"
    shift # past argument=value
    ;;
    -p0=*|--peaks0=*)
    peaksF="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--names=*)
    sampNames="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--out=*)
    pathOut="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annotation=*)
    refFlat="${i#*=}"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

### ===================
##    CHECK ARGUMENTS
# =====================

echo "Start checking arguments..."
check=1

##-c|--config=*)
if [ ! -f ${conf} ]
then
	echo "Error: the provided configuration file do not exist. Please verify"
	echo "You provided "${conf}
	check=0
fi

##-pF=*|--peakFolder=*)
if [[ ${pathTsv} != "" && ${pathTsv} != */ ]]
then
	#Add / to the folder if it is needed
	pathTsv=${pathTsv}/
fi
if [[ ${pathTsv} != "" && ! -d ${pathTsv} ]]
then
	echo "Error: the provided peakFolder do not exist. Please verify"
	echo "You provided "${pathTsv}
	check=0
fi

##-p1=*|--peaks1=*)
if [ -z ${peaksf} ]
then
	echo "Error: peaks in condition 1 not provided. Please fill -p1= or --peaks1 argument"
	check=0
else
	if [ ! -f ${pathTsv}${peaksf} ]
	then
		echo "Error: peaks file does not exist. Please verify -p1= or --peaks1 argument"
		echo "You provided "${pathTsv}${peaksf}
		check=0
	else
		#Check if the bed is made of four columns
		nbCol=$( head -n 1 ${pathTsv}${peaksf} | grep -o -P '\t' | wc -l )
		nbCol=$(( nbCol + 1 ))
		if (( ${nbCol} > 4 ))
		then
			cmd="cut -f1,2,3,4 "${pathTsv}${peaksf}" > "${pathTsv}${peaksf}".4col"
			echo ${cmd}
			eval ${cmd}
			peaksf=${peaksf}".4col"
		elif (( ${nbCol} == 3 ))
		then
			echo "awk \'BEGIN{OFS=FS=\"\t\"}{print $1\"\t\"$2\"\t\"$3\"\t\"$1\"-\"$2\"-\"$3}' "${pathTsv}${peaksf}" > "${pathTsv}${peaksf}".4col"
			awk 'BEGIN{OFS=FS="\t"}{print $1"\t"$2"\t"$3"\t"$1"-"$2"-"$3}' ${pathTsv}${peaksf} > ${pathTsv}${peaksf}.4col
			peaksf=${peaksf}".4col"
		elif (( ${nbCol} < 3 ))
		then
			echo "Error: peak file need at least chromosom, start and stop columns. Please verify "${pathTsv}${peaksf}" file"
			check=0
		fi
	fi
fi

##-p0=*|--peaks0=*)
if [ -z ${peaksF} ]
then
 	echo "Error: peaks in condition 0 not provided. Please fill -p0= or --peaks0 argument"
 	check=0
else
	if [ ! -f ${pathTsv}${peaksF} ]
	then
		echo "Error: peaks file does not exist. Please verify -p0= or --peaks0 argument"
		echo "You provided "${pathTsv}${peaksF}
		check=0
	else
		#Check if the bed is made of four columns
		nbCol=$( head -n 1 ${pathTsv}${peaksF} | grep -o -P '\t' | wc -l )
		nbCol=$(( nbCol + 1 ))
		if (( ${nbCol} > 4 ))
		then
			cmd="cut -f1,2,3,4 "${pathTsv}${peaksF}" > "${pathTsv}${peaksF}".4col"
			echo ${cmd}
			eval ${cmd}
			peaksF=${peaksF}".4col"
		elif (( ${nbCol} == 3 ))
		then
			echo "awk \'BEGIN{OFS=FS=\"\t\"}{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$1\"-\"\$2\"-\"\$3}\' "${pathTsv}${peaksF}" > "${pathTsv}${peaksF}".4col"
			awk 'BEGIN{OFS=FS="\t"}{print $1"\t"$2"\t"$3"\t"$1"-"$2"-"$3}' ${pathTsv}${peaksF} > ${pathTsv}${peaksF}.4col
			peaksF=${peaksF}".4col"
		elif (( ${nbCol} < 3 ))
		then
			echo "Error: peak file need at least chromosom, start and stop columns. Please verify "${pathTsv}${peaksF}" file"
			check=0
		fi
	fi
fi

##-n=*|--names=*)
if [ ! -z ${sampNames} ]
then
	if [ ! ${sampNames} == *bL1_* ]
	then
		echo "Adding tags to names"
		sampNames=$( echo "bL1_"${sampNames} | sed 's/,/°bL0_/'  | sed 's/,/°bI1_/' | sed 's/,/°bI0_/' | sed 's/°/,/g')
	fi
fi

##-o=*|--out=*)
if [ -z ${pathOut} ]
then
	echo "Error: output directory not provided. Please fill -o= or --out argument"
fi
if [[ ${pathOut} != "" && ${pathOut} != */ ]]
then
	#Add / to the folder if it is absent
	pathOut=${pathOut}/
fi
if [ ! -d ${pathOut} ]
then
	echo "Creating output directory"
	mkdir ${pathOut}
fi

##-a=*|--annotation=*)
if [ -z ${refFlat} ]
then
	echo "Error: refFlat annotation file not provided. Please fill -a= or --annotation argument"
else
	if [ ! -f ${refFlat} ]
	then
		echo "Error: refFlat annotation file does not exist. Please verify -a= or --annotation argument"
		echo "You provided "${refFlat}
		check=0
	fi
fi

##echo some
echo "IP_CONDITION1       = ${peaksf}"
echo "IP_CONDITION0         = ${peaksF}"
echo "OUT               = ${pathOut}"
echo "... Done"

### ====================
##    SOURCE CONF FILE
# ======================

echo "source ${conf}"
source ${conf}

### ========
##    MAIN
# ==========

if (( check > 0 ))
then
	echo "All arguments passed the checking."

	## Prepare
	mkdir ${pathOut}ComparPeaks/
	pathComp=${pathOut}ComparPeaks/
	mkdir ${pathOut}Diff/
	pathDiff=${pathOut}Diff/
	peaksf=${pathTsv}${peaksf}
	peaksF=${pathTsv}${peaksF}
	#Normally only the first sed is usefull the other are for compatibility
	n=$( basename ${f} | sed -re 's/\_[a-zA-Z0-9]*\-(map|less)//' | sed 's/\/accepted_hits//' | sed 's/\.bam//' | sed 's/lessDuplicates//' | sed 's/sort//' | sed 's/_//g')
	N=$( basename ${F} | sed -re 's/\_[a-zA-Z0-9]*\-(map|less)//' | sed 's/\/accepted_hits//' | sed 's/\.bam//' | sed 's/lessDuplicates//' | sed 's/sort//' | sed 's/_//g')
	outname=${n}_VS_${N}
	outComp=${pathComp}${outname}
	outDiff=${pathDiff}${outname}.csv
	tmpGb=${pathOut}${outname}_$( basename ${refFlat} ).gb

	## Run "compare_bed.py" to realise spatial comparison and peak union
	if [[ ${FLAG_COMPARE_BED} -eq 1 ]]
	then
		cmd="python ${COMPARE_BED} -f ${peaksf} -F ${peaksF} -o ${outComp} -n ${n} -N ${N}"
		echo ${cmd}
		${cmd}
	fi

	## Annotate the peak union
	if [[ ${FLAG_ANNOT_UNION} -eq 1 ]]
	then
		#"annotate-peaks"
		cmd="${ANNOTE_BED} ${outComp}_mergebed.tsv ${outComp}_annotatedpeaks.txt ${refFlat} name"
		echo ${cmd}
		${cmd}
		#"Annote_Modifier.R"
		cmd="${ANNOT_MODIF} -i=${outComp}_annotatedpeaks.txt -o=${outComp}_Reannotation.txt -leftStart=-200 -rightStart=200 -leftStop=-200 -rightStop=200 -draw=5UTR,near_Start,CDS,ncExon,intron,ncIntron,near_Stop,3UTR,Intergenic,Multiple"
		echo ${cmd}
		${cmd}
	fi
fi
