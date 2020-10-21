#!/usr/bin/env bash

###################################################################
# Description: Bash runner for Meme and Ame motif tools           #
# Date: 2019/03                                                   #
# Author: Martin BIZET                                            #
###################################################################

### ======================
##    Default parameters
# ========================

source ${THESE_TOOLS}/"dependencies.sh"

RunIUPAC=0
Tool="meme-chip"

#meme parameters
seed=555
nmot=10
minw=5
maxw=12

#ame parameters
score="avg" #(default on website) "totalhits" #(default of the tool [seems less efficient])
method="ranksum" #(default on website) "fisher" #(default of the tool [seems less efficient])
pv=0.05
alph="RNAalphabet.alph"

#Parameters optimised with Ame pval
bg=1
ntop=2500
elong="Resize"
ext=150
resize=500

### ================
##    ARGUMENTS
# ==================

HELP="
\nDESCRIPTION:
\n
\n\tRun Meme or Ame tools for peak positions assuming transcriptome experiment (e.g. meRIP)
\n
\n-pF=|--peakFolder= folder containing the peaks.
\n-i=|--input= file containing the peaks. Should be absolute path if -pF is not provided. /!\ caution: it expects a standard bed file (i.e. chromosom | start | stop | name | fold change or any score for ranking) /!\ [Mandatory]
\n-e=|--expected= file containing the expected/control peaks. Should be absolute path if -pF is not provided. [Mandatory with neg=peakNeg]
\n-g=|--genome= genome in fasta format. [Mandatory]
\n-a=|--annotation= list of transcripts with TSS and TTS in bed format (typically \'geneBody.bed\' file). [Mandatory]
\n-n=|--namePrefix= prefix of files generated in ouput directory. [default= basename of input]
\n-o=|--outFolder= ouput directory. [Mandatory]
\n-M=|--motif= string in degenerated alphabet (several motifs can be provided using coma-separation) or .meme file (absolute path; .meme extension mandatory) of the motif to evaluate using Ame. If provided Ame is run instead of meme-chip.
\n-b=|--backgroundFa= Fasta (typically obtain from bam file using \'samtools view -F 4 <bam file> | awk '{print \">\"$1\"\n\"$10}'\' commmand) to extract background Markov model. [Mandatory if bg is not 0]
\n-bg=|--backgroundOrder= order of the background Markov model. [default=1]
\n-neg=|--negativeControlMethod= method for the negative control. Please use \'peakNeg\', \'shuffle\' or \'memeDefault\'. [default=\'peakNeg\' if -e is provided, \'shuffle\' otherwise]
\n-top=|--maximalTopSequences= maximal number of top-ranked sequences to provide to Meme/Ame. [default=2500]
\n-elM=|--elongationMethod= How to elongate the sequence from the peak. Please use \'Resize\' or \'Extend\'. \'Resize\' generate peaks of size -reS (so peaks may be shorten) while \'Extend\' extend of -ext bases from both peak borders [default=\'Resize\']
\n-reS=|--peakResizing= how much should be the peak size (when -elM is set to \'Resize\'). [default=500]
\n-ext=|--peakExtension= how much the peak size should be extended on both size (when -elM is set to \'Extend\'). [default=150 (leading to 500bp sequences with standard 200bp peaks)]
\n-se=|--seed= random seed for Meme. [default=555]
\n-nMo=|--numberOfMotifs= maximal number of motifs for Meme and Dreme. [default=10]
\n-min=|--minWindowSize= minimal window size for Meme. [default=5]
\n-max=|--maxWindowSize= maximal window size for Meme. [default=12]
\n-sc=|--score= scoring method for Ame. Please use \'avg\', \'max\', \'sum\' or \'totalhits'. [default=\'avg\']
\n-met=|--method= statistical method for Ame. Please use \'fisher\', \'mhg\', \'4dmhg\', \'ranksum\', \'linreg\' or \'spearman\'. [default=\'ranksum\']
\n-pv=|--pvalue-report-threshold= pvalue threshold for reporting Ame result. [default=0.05]
\n-alp=|--alphabet= filename of the reference alphabet used with Ame. Extract from ALPHABET_FOLDER if absolute path is not provided. [default=\'RNAalphabet.alph\']
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
    -pF=*|--peakFolder==*)
    pathPeaks="${i#*=}"
    pathPeaks=${pathPeaks}/
    shift # past argument=value
    ;;
    -i=*|--input=*)
    peak="${i#*=}"
    shift # past argument=value
    ;;
    -e=*|--expected=*)
    peakNeg="${i#*=}"
    shift # past argument=value
    ;;
    -g=*|--genome=*)
    genome="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annotation=*)
    gb="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--outFolder=*)
    pathRes="${i#*=}"
    shift # past argument=value
    ;;
    -M=*|--motif==*)
    Motif="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--namePrefix=*)
    nam="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--backgroundFa=*)
    bdgSeq="${i#*=}"
    shift # past argument=value
    ;;
    -bg=*|--backgroundOrder=*)
    bg="${i#*=}"
    shift # past argument=value
    ;;
    -neg=*|--negativeControlMethod=*)
    neg="${i#*=}"
    shift # past argument=value
    ;;
    -top=*|--maximalTopSequences=*)
    ntop="${i#*=}"
    shift # past argument=value
    ;;
    -elM=*|--elongationMethod=*)
    elong="${i#*=}"
    shift # past argument=value
    ;;
    -reS=*|--peakResizing=*)
    resize="${i#*=}"
    shift # past argument=value
    ;;
    -ext=*|--peakExtension=*)
    ext="${i#*=}"
    shift # past argument=value
    ;;
    -se=*|--seed=*)
    seed="${i#*=}"
    shift # past argument=value
    ;;
    -nMo=*|--numberOfMotifs=*)
    nmot="${i#*=}"
    shift # past argument=value
    ;;
    -min=*|--minWindowSize=*)
    minw="${i#*=}"
    shift # past argument=value
    ;;
    -max=*|--maxWindowSize=*)
    maxw="${i#*=}"
    shift # past argument=value
    ;;
    -sc=*|--score=*)
    score="${i#*=}"
    shift # past argument=value
    ;;
	-met=*|--method=*)
    method="${i#*=}"
    shift # past argument=value
    ;;
	-pv=*|--pvalue-report-threshold=*)
    pv="${i#*=}"
    shift # past argument=value
    ;;
	-alp=*|--alphabet=*)
    alph="${i#*=}"
    shift # past argument=value
    ;;
    *)
          # unknown option
    ;;
esac
done

### ====================
##    CHECK ARGUMENTS
# ======================

echo "Start checking arguments..."
check=1
re='^[0-9]+$'

##-pF=|--peakFolder=
if [ ! -z ${pathPeaks} ]
then
	if [ ! -d ${pathPeaks} ]
	then
		echo "Error: the provided peak folder do not exist. Please verify"
		echo "You provided "${pathPeaks}
		check=0
	fi
fi

##-i=|--input==
if [ -z ${peak} ]
then
	echo "Error: mandatory argument -i is not provided"
	check=0
elif [ ! -f ${peak} ]
then
	if [ ! -f ${pathPeaks}${peak} ]
	then
		echo "Error: the provided peaks file do not exist. Please verify"
		if [ -z ${pathPeaks} ]
		then
			echo "You provided "${peak}
		else
			echo "You provided "${pathPeaks}${peak}
		fi
		check=0
	fi
fi

##-e=|--expected=
if [ -z ${neg} ]
then
	if [ -z ${peakNeg} ]
	then
		neg="shuffle"
	else
		neg="peakNeg"
	fi
fi
if [[ ${neg} == "peakNeg" ]]
then
	if [ -z ${peakNeg} ]
	then
		echo "Error: argument -e is not provided (mandatory with neg=peakNeg)"
		check=0
	elif [ ! -f ${peakNeg} ]
	then
		if [ ! -f ${pathPeaks}${peakNeg} ]
		then
			echo "Error: the provided negative-peaks file do not exist. Please verify"
			if [ -z ${pathPeaks} ]
			then
				echo "You provided "${peakNeg}
			else
				echo "You provided "${pathPeaks}${peakNeg}
			fi
			check=0
		else
			peakNeg=${pathPeaks}${peakNeg}
		fi
	fi
fi

##-g=|--genome=
if [ -z ${genome} ]
then
	echo "Error: mandatory argument -g is not provided"
	check=0
elif [ ! -f ${genome} ]
then
	echo "Error: the provided genome fasta file do not exist. Please verify"
	echo "You provided "${genome}
	check=0
fi

##-a=|--annotation=
if [ -z ${gb} ]
then
	echo "Error: mandatory argument -a is not provided"
	check=0
elif [ ! -f ${gb} ]
then
	echo "Error: the provided annotation file do not exist. Please verify"
	echo "You provided "${gb}
	check=0
fi

##-o=|--outFolder=
if [ -z ${pathRes} ]
then
	echo "Error: mandatory argument -o is not provided"
	check=0
elif [ ! -d ${pathRes} ]
then
	echo "Error: the provided output folder do not exist. Please verify"
	echo "You provided "${pathRes}
	check=0
fi

##-M=|--motif= string of the motif in degenerated alphabet to evaluate using Ame. If provided Ame is run instead of meme-chip.
if [ ! -z ${Motif} ]
then
	Tool="ame"
	if [[ ${Motif} != *.meme ]]
	then
		if [[ $( echo ${Motif} | sed 's/,//g' ) =~ ^[A-IK-NP-Ya-ik-np-y]+$ ]]
		then
			RunIUPAC=1
		else
			echo "Wrong motif: should only contain IUPAC letters"
                        check=0
		fi
	fi
fi

##-n=|--namePrefix=
if [ -z ${nam} ]
then
	nam=$(basename ${peak} .bed)
fi

##-b=|--backgroundFa=
if (( ${bg} > 0 ))
then
	if [ -z ${bdgSeq} ]
	then
		echo "Error: argument -b is not provided (mandatory with bg > 0)"
		check=0
	elif [ ! -f ${bdgSeq} ]
	then
		echo "Error: the provided background fasta file do not exist. Please verify"
		echo "You provided "${bdgSeq}
		check=0
	fi
fi

##-bg=|--backgroundOrder=
if [[ ! ${bg} =~ ${re} ]]
then
   echo "Error: -bg should be a number"
   check=0
fi

##-neg=|--negativeControlMethod= method for the negative control. Please use \'peakNeg\', \'shuffle\' or \'memeDefault\'. [default=\'peakNeg\' if -e is provided, \'shuffle\' otherwise]
if [[ ! ${neg} == "peakNeg" ]]
then
	if [[ ! ${neg} == "shuffle" ]]
	then
		if [[ ! ${neg} == "memeDefault" ]]
		then
			echo "Error: argument -neg should be peakNeg, shuffle or memeDefault"
			check=0
		fi
	fi
fi

##-top=|--maximalTopSequences=
if [[ ! ${ntop} =~ ${re} ]]
then
   echo "Error: -top should be a number"
   check=0
fi

##-elM=|--elongationMethod=
if [[ ! ${elong} == "Extend" ]]
then
	if [[ ! ${elong} == "Resize" ]]
	then
		echo "Error: argument -elM should be Extend or Resize"
		check=0
	fi
fi

##-ext=|--peakExtension=
if [[ ! ${ext} =~ ${re} ]]
then
   echo "Error: -ext should be a number"
   check=0
fi

##-reS=|--peakResizing=
if [[ ! ${resize} =~ ${re} ]]
then
   echo "Error: -ext should be a number"
   check=0
fi

##-se=|--seed=
if [[ ! ${seed} =~ ${re} ]]
then
   echo "Error: -se should be a number"
   check=0
fi

##-nMo=|--numberOfMotifs=
if [[ ! ${nmot} =~ ${re} ]]
then
   echo "Error: -nMo should be a number"
   check=0
fi

##-min=|--minWindowSize=
if [[ ! ${minw} =~ ${re} ]]
then
   echo "Error: -min should be a number"
   check=0
fi

##-max=|--maxWindowSize=
if [[ ! ${maxw} =~ ${re} ]]
then
   echo "Error: -max should be a number"
   check=0
fi

##-sc=|--score= scoring method for Ame. Please use \'avg\', \'max\', \'sum\' or \'totalhits'. [default=\'avg\']
if [[ ${score} != "avg" ]]
then
	if [[ ${score} != "max" ]]
	then
		if [[ ${score} != "sum" ]]
		then
			if [[ ${score} != "totalhits" ]]
			then
				echo "Error: argument -sc should be avg, max, sum or totalhits"
				check=0
			fi
		fi
	fi
fi

##-met=|--method= statistical method for Ame. Please use \'fisher\', \'mhg\', \'4dmhg\', \'ranksum\', \'linreg\' or \'spearman\'. [default=\'ranksum\']
if [[ ${method} != "fisher" ]]
then
	if [[ ${method} != "mhg" ]]
	then
		if [[ ${method} != "4dmhg" ]]
		then
			if [[ ${method} != "ranksum" ]]
			then
				if [[ ${method} != "linreg" ]]
				then
					if [[ ${method} != "spearman" ]]
					then
						echo "Error: argument -met should be avg, max, sum or totalhits"
						check=0
					fi
				fi
			fi
		fi
	fi
fi

##-pv=|--pvalue-report-threshold= pvalue threshold for reporting Ame result. [default=0.05]
if [[ $( echo "${pv} > 1" | bc) == 1 ]]
then
   echo "Error: -pv should be a number between 0 and 1"
   check=0
fi

##-alp=|--alphabet= filename of the reference alphabet used with Ame. Extract from ALPHABET_FOLDER if absolute path is not provided. [default=\'RNAalphabet.alph\']
if [ ! -f ${alph} ]
then
	if [ -f ${ALPHABET_PATH}${alph} ]
	then
		alph=${ALPHABET_PATH}${alph}
	else
		echo "Error: the provided alphabet file do not exist. Please verify"
		echo "You provided "${alph}
	fi
fi

###===========
##    MAIN
#=============

if [[ ${check} == 0 ]]
then
	exit 1
fi

echo "All arguments passed the checking."

#Generat .meme file
if [[ ${RunIUPAC} == 1 ]]
then
	cmd=${MEME_PATH}"iupac2meme -alph "${alph}" "$( echo ${Motif} | sed 's/,/ /g' )" > "${pathRes}"Motif.meme"
	echo ${cmd}
	eval ${cmd}
	Motif=${pathRes}"Motif.meme"
fi
fnam=${pathRes}${nam}
#Use bedtools intersect to associate the strand of the known transcripts to sequence

# /!\ uniq only supress duplicates lines WHEN THEY ARE ADJACENT !!! Be carefull to sort BEFORE doing uniq
cmd="bedtools intersect -wb -a "${gb}" -b "${pathPeaks}${peak}" | awk 'BEGIN{FS=OFS=\"\t\"}{ print \$7,\$8,\$9,\$10,\$11,\$6}' | sort -k5nr,5 | uniq | head -n "${ntop}" > "${fnam}"_STRANDED.bed"
echo ${cmd}
eval ${cmd}

#Elong
if [[ ${elong} == "Extend" ]]
then
	cmd="awk 'BEGIN{FS=OFS=\"\t\"}{s=\$2-"${ext}";if(s < 0) {s = 0};print \$1,s,\$3+"${ext}",\$4,\$5,\$6}' "${fnam}"_STRANDED.bed > "${fnam}"_ELONG.bed"
else
	cmd="awk 'BEGIN{FS=OFS=\"\t\"}{c=int((\$2+\$3)/2); s=c-"${resize}"/2;if(s < 0) {s = 0};print \$1,s,s+"${resize}",\$4,\$5,\$6}' "${fnam}"_STRANDED.bed > "${fnam}"_ELONG.bed"
fi
echo ${cmd}
eval ${cmd}

#Extract sequences:	bedtools getfasta -s
cmd="bedtools getfasta -s -fi "${genome}" -bed "${fnam}"_ELONG.bed -fo "${fnam}".fa"
echo ${cmd}
${cmd}

#Convert in RNA sequence (T -> U)
sequences=${fnam}".fa"
cmd="sed -i -e 'y/tT/uU/' "${sequences}
echo ${cmd}
eval ${cmd}

#Basic meme command
#MEME-ChiP algorithm: meme-chip	<input_FASTA_file> -bfile <background_file> -meme-nmotifs 10 -meme-maxsize 1000000 -meme-minw 5 -meme-maxw 12 -dreme-m 10 -norc -o <output_directory>
#-norc for stranded
if [[ ${Tool} == "meme-chip" ]]
then
	memecmd=${MEME_PATH}"meme-chip -oc "${pathRes}${nam}" -rna -index-name "${nam}".html -order "${bg}" -seed "${seed}" -meme-nmotifs "${nmot}" -meme-minw "${minw}" -meme-maxw "${maxw}" -dreme-m "${nmot}" -norc"
else
	amecmd=${MEME_PATH}"ame --oc "${pathRes}${nam}
fi

#Background parameter
if [[ ${bg} != "0" ]]
then
	modFile=${nam}"-INPUT_ord"${bg}".mod"
	if [[ ! -f ${pathRes}${modFile} ]]
	then
		#Prepare background model
		cmd=${MEME_PATH}"fasta-get-markov -rna -m "${bg}" "${bdgSeq}" "${pathRes}${modFile}
		echo ${cmd}
		${cmd}
	else
		echo "Warning: "${modFile}" was already present in "${pathRes}". Its re-generation was skipped."
	fi
	if [[ ${Tool} == "meme-chip" ]]
	then
		memecmd=${memecmd}" -bfile "${pathRes}${modFile}
	else
		amecmd=${amecmd}" --bgformat 2 --bgfile "${pathRes}${modFile}
	fi
fi

#Negative peak parameter
if [[ ${neg} == "shuffle" ]]
then
	#Prepare negative file
	esequences=${fnam}"_neg.fa"
	cmd=${MEME_PATH}"fasta-shuffle-letters "${fnam}".fa "${esequences}
	echo ${cmd}
	${cmd}
	#meme-chip
	if [[ ${Tool} == "meme-chip" ]]
	then
		memecmd=${memecmd}" -neg "${esequences}
	else
		amecmd=${amecmd}" --control "${esequences}
	fi
elif [[ ${neg} == "peakNeg" ]]
then
	#Prepare negative file
	enam=${pathRes}${nam}"_Expect"
	# /!\ uniq only supress duplicates lines WHEN THEY ARE ADJACENT !!! Be carefull to sort BEFORE doing uniq
	cmd="bedtools intersect -wb -a "${gb}" -b "${peakNeg}" | awk 'BEGIN{FS=OFS=\"\t\"}{ print \$7,\$8,\$9,\$10,\$11,\$6}' | sort -k5nr,5 | uniq | head -n "${ntop}" > "${enam}"_STRANDED.bed"
	echo ${cmd}
	eval ${cmd}
	negFile=${enam}"_ELONG"
	cmd="awk 'BEGIN{FS=OFS=\"\t\"}{s=\$2-"${ext}";if(s < 0) {s = 0};print \$1,s,\$3+"${ext}",\$4,\$5,\$6}' "${enam}"_STRANDED.bed > "${negFile}".bed"
	echo ${cmd}
	eval ${cmd}
	cmd="bedtools getfasta -s -fi "${genome}" -bed "${negFile}".bed -fo "${negFile}".fa"
	echo ${cmd}
	${cmd}
	esequences=${negFile}".fa"
	cmd="sed -i -e 'y/tT/uU/' "${esequences}
	echo ${cmd}
	eval ${cmd}
	#meme-chip
	if [[ ${Tool} == "meme-chip" ]]
	then
		memecmd=${memecmd}" -neg "${esequences}
	else
		amecmd=${amecmd}" --control "${esequences}
	fi
fi

#Run
if [[ ${Tool} == "meme-chip" ]]
then
	cmd=${memecmd}" "${sequences}
else
	cmd=${amecmd}" --scoring "${score}" --method "${method}" --pvalue-report-threshold "${pv}" "${sequences}" "${Motif}
fi
echo ${cmd}
eval ${cmd}

echo "+-+-+-+-+-"
