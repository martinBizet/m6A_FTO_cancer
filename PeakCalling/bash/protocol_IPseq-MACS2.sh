#!/usr/bin/env bash

# dependencies
echo "source "${THESE_TOOLS}"/dependencies.sh"
source ${THESE_TOOLS}/dependencies.sh

# ==================
# ARGUMENTS
# ==================

HELP="
\nDESCRIPTION:
\n
\n\tIP-seq analysis with MACS2
\n
\nUSAGE:
\n
\n\t${0} -I=<bam of IP> -C=<bam of INPUT> -o=<path to directory for output> -p=<prefexie of output [default=IP-seq]> -q=<MACS2 q-value [default=0.05]> -a=<ref.flat annotation> -ts=<transcriptome size [default=extract from -a keeping introns]> -s=<half size of summits-sentered peaks [default=100]> -kd[if provided MACS2 keep the duplicates] -es=<MACS2 extsize [if provided MACS2 --nomodel option is used]>
\n
\nTEST:
\n\t[romy@Pipeline ~]$ bash PeakCalling/bash/protocol_IPseq-MACS2.sh -I=<IP sorted bam file> -C=<INPUT sorted bam file> -o=<output directory> -a=<annotation refFlat format>
\n
\nreal    1m1.664s
\nuser    0m57.694s
\nsys 0m1.704s

\n
"
DATE=$( date +%Y%m%d-%H%M%S-%s%N )

#Default
prefix="IP-seq"
MACS2q=0.05
MACS2kd=""
MACS2sh=""
MACS2TSVhalfsize=100

# parse
for i in "$@"
do
case $i in
    -h*|--help*)
    echo -e ${HELP};
    exit 0
    ;;
    -I=*|--IP=*)
    f_condition="${i#*=}"
    shift # past argument=value
    ;;
    -C=*|--CONTROL=*)
    f_control="${i#*=}"
    shift # past argument=value
    ;;
    -o=*|--out=*)
    out_dir="${i#*=}"
    shift # past argument=value
    ;;
    -q=*|--MACS2qvalue=*)
    MACS2q="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--prefix=*)
    prefix="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--annotation=*)
    fannot="${i#*=}"
    shift # past argument=value
    ;;
    -ts=*|--TranscriptomeSize=*)
    size_transcriptom="${i#*=}"
    shift # past argument=value
    ;;
    -s=*|--halfSize=*)
    MACS2TSVhalfsize="${i#*=}"
    shift # past argument=value
    ;;
    -kd*|--MACS2keepDuplicates*)
    MACS2kd=" --keep-dup all"
    shift # past argument=value
    ;;
    -es*|--extsize=*)
    MACS2sh="${i#*=}"
    if [[ ${MACS2sh} =~ ^[0-9]+$ ]]
    then
        MACS2sh=" --nomodel --extsize "${MACS2sh}
    else
        echo "[ERROR] MACS2 extsize should be numerical"
        exit 1
    fi
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

echo "F_CONDITION       = ${f_condition}"
echo "F_CONTROL         = ${f_control}"
echo "PREFIX            = ${prefix}"
echo "OUT               = ${out_dir}"

dir_ori=$(pwd)
cd $out_dir

## Define name ##
echo "[LOG] Define name of output"
#if accepted_hits file use folder name
#otherwise use filename (for compatibility but normally not used anymore)
basename1=$( basename ${f_condition} )
if [[ ${basename1} == accepted_hits* ]]
then
	basename1=$( basename $(dirname ${f_condition} )).bam
fi
basename2=$( basename ${f_control} )
if [[ ${basename2} == accepted_hits* ]]
then
	basename2=$( basename $(dirname ${f_control} )).bam
fi
#commun_prefixe=$( printf "%s\n%s\n" "$basename1" "$basename2" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' )
commun_prefixe=${prefix}
name_out=${commun_prefixe}-t$( echo $basename1 | sed 's/'${commun_prefixe}'//' )-c$( echo $basename2 | sed 's/'${commun_prefixe}'//' )
name_expect=${commun_prefixe}-e$( echo $basename2 | sed 's/'${commun_prefixe}'//' )

## Compute transcriptom size ##
if [[ -z ${size_transcriptom} ]]
then
	echo "[LOG] Compute size of transcriptom WITH INTRON !!!!"
	# format annotation to bed
	commande="awk 'BEGIN{FS=OFS=\"\t\"}{ if ( \$0 !~ /^#/ ){ print \$3,\$5,\$6}}' ${fannot} | bedtools sort -i > tmp.${DATE}.bed ; bedtools merge -i tmp.${DATE}.bed | awk 'BEGIN{total=0}{total=total+\$3-\$2}END{print total}'; rm tmp.${DATE}.bed "
	echo "[LOG] ... command: ${commande}"
	size_transcriptom=$( awk 'BEGIN{FS=OFS="\t"}{ if ( $0 !~ /^#/ ){ print $3,$5,$6}}' ${fannot} | bedtools sort -i > tmp.${DATE}.bed ; bedtools merge -i tmp.${DATE}.bed | awk 'BEGIN{total=0}{total=total+$3-$2}END{print total}'; rm tmp.${DATE}.bed )
fi
echo "[LOG] ... size_transcriptom = ${size_transcriptom}"

## Peak calling ##
echo "[LOG] Detect enriched regions"
mkdir MACS2
cd MACS2
commande="${MACS2} callpeak -t ${f_condition} -c ${f_control} -g ${size_transcriptom} --outdir="$(pwd)"/ -f AUTO -n ${name_out} -B -q ${MACS2q}${MACS2sh}${MACS2kd} --call-summits"
echo "[LOG] ... command: ${commande}"
$commande
out_macs=$( ls ${out_dir}/MACS2/${name_out}*.narrowPeak )
if [[ -z $out_macs ]]
then
    echo "[ERROR] file don't exist: ${out_dir}/MACS2/${name_out}*.narrowPeak"
    cd $dir_ori
    exit 1
fi
# get stats
nb_peaks=$( awk '{print $1":"$2"-"$3}' $out_macs | sort | uniq |  wc -l ); 
nb_summits=$( wc -l $out_macs | awk '{print $1"\t"$2}' ); echo -e "${d1}\t${d2}\t${nb_peaks}\t${nb_summits}"
echo "[LOG] ... nb_peaks: ${nb_peaks}"
echo "[LOG] ... nb_summits: ${nb_summits}"

## Expected ##
echo "[LOG] Detect expected enriched regions"
out_macs=$( ls ${out_dir}/MACS2/${name_expect}*.narrowPeak )
if [[ -z $out_macs ]]
then
	commande="${MACS2} callpeak -t ${f_control} -g ${size_transcriptom} --outdir="$(pwd)"/ -f AUTO -n ${name_expect} -B -q ${MACS2q}${MACS2sh}${MACS2kd} --call-summits"
	echo "[LOG] ... command: ${commande}"
	$commande
	out_macs=$( ls ${out_dir}/MACS2/${name_expect}*.narrowPeak )
	if [[ -z $out_macs ]]
	then
		echo "[ERROR] file don't exist: ${out_dir}/MACS2/${name_expect}*.narrowPeak"
		cd $dir_ori
		exit 1
	fi
else
	echo "[INFO] ${out_macs} was already present and has not been recomputed"
fi
# get stats
nb_peaks=$( awk '{print $1":"$2"-"$3}' $out_macs | sort | uniq |  wc -l );
nb_summits=$( wc -l $out_macs | awk '{print $1"\t"$2}' ); echo -e "${d1}\t${d2}\t${nb_peaks}\t${nb_summits}"
echo "[LOG] ... nb_peaks: ${nb_peaks}"
echo "[LOG] ... nb_summits: ${nb_summits}"
cd $out_dir

## MACS2TSV ##
echo "[LOG] Define regions around summits"
mkdir MACS2TSV
cd MACS2TSV
#peaks
commande="${MACS2TSV} -i=${out_dir}/MACS2/${name_out} -o="$(pwd)"/ -s="${MACS2TSVhalfsize}
echo "[LOG] ... command: ${commande}"
${commande}
#expected
commande="${MACS2TSV} -i=${out_dir}/MACS2/${name_expect} -o="$(pwd)"/ -s="${MACS2TSVhalfsize}
echo "[LOG] ... command: ${commande}"
${commande}
cd $out_dir

## Annotate ##
#peaks
echo "[LOG] Annote regions"
out=$(pwd)/ANNOTATION/${name_out}_annoted.xls
mkdir -p $( dirname $out )
cd $( dirname $out )
commande="time ${ANNOTE} ${out_dir}/MACS2/${name_out}_peaks.xls ${out} ${fannot}"
echo "[LOG] ... command: ${commande}"
$commande
cd $out_dir
#expected
echo "[LOG] Annote expected regions"
outEx=$(pwd)/ANNOTATION/${name_expect}_annoted.xls
mkdir -p $( dirname $outEx )
cd $( dirname $outEx )
commande="time ${ANNOTE} ${out_dir}/MACS2/${name_expect}_peaks.xls ${outEx} ${fannot}"
echo "[LOG] ... command: ${commande}"
$commande
cd $out_dir

## Make annotation report ##
echo "[LOG] Re-annote regions and generate graph"
f=$out
e=$outEx
out=$(pwd)/RE-ANNOTATE/$( basename $f | sed 's/_annoted.xls/_re-annoted/g' )
mkdir -p $( dirname $out )
cd $( dirname $out )
commande="${ANNOT_MODIF} -i=${f} -e=${e} -o=${out} -leftStart=-200 -rightStart=200 -leftStop=-200 -rightStop=200 -draw=5UTR,near_Start,CDS,ncExon,intron,ncIntron,near_Stop,3UTR,Intergenic,Multiple"
echo "[LOG] ... command: ${commande}"
$commande
cd $out_dir

cd $dir_ori


