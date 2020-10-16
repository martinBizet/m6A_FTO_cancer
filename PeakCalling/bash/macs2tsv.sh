#!/usr/bin/env bash                                                                                                                                                          
 
############################################################
# Descrition: convert .xls from MACS2 to tsv
# Date: 2016/06
# Author: Romy CHEN-MIN-TAO
############################################################
 
# ==================
# VARIABLES
# ==================
 
DATE=$( date +%Y%m%d-%H%M%S-%s%N )
DIR_SCRIPT=$( dirname $0 )
DIR_ORI=$(pwd)
SIZ=100
 
# ==================
# DEPENDENCIES    
# ==================
 
source ${DIR_SCRIPT}/multithread_management.sh
source ${DIR_SCRIPT}/check.sh
 
# ==================
# FUNCTIONS
# ==================

function macs2tsv(){
    text_help="
    \n# DESCRIPTION:
    \n\tconvert macs2 .xls file to tsv 
    \n# ARGUMENTS:
    \n\t-h or --help     : get help
    \n\t-i or --input    : path to directory containing .xls from MACS2
    \n\t-o or --out      : path to directory where you want storage results
    \n\t-s or --halfSize : half-size of sumits-centered peaks in bp [default=100]
    \n"
    # Parse aguments
    first_command="$@"
    for i in "$@"
    do
    case $i in
        -h*|--help*)                 
        echo -e ${text_help};                                           
        exit 0
        ;;
        -i=*|--input=*)
        DIR_XLS="${i#*=}"
        shift # past argument=value
        ;;
        -o=*|--out=*)  
        OUT_DIR="${i#*=}"                                          
        shift # past argument=value                                
        ;;
        -s=*|--halfSize=*)  
        SIZ="${i#*=}"                                          
        shift # past argument=value                                
        ;;
        *)
        echo "[ERROR] option unknow: '${i}'"
        echo -e ${text_help}; # I can't print it -_-
        exit 0
        # unknown option
        ;;
    esac
    done
    # run
    # - check if OUT_DIR exist
    if [[ ! -d $OUT_DIR ]]
    then
        mkdir -p $OUT_DIR
    fi
    cd $OUT_DIR
    # scan xls files
    for f in $( ls ${DIR_XLS}*.xls ) 
    do           
        echo "[LOG] process ${f}" >&2 
        basename_f=$( basename $f .xls )
        # remove comment lines and header 
        commande="sed '/^[ \t]*$/d;s/#.*//g' $f | strings | sed '1d' > ${basename_f}_nocomments.tsv"
        echo "[LOG] ... command: ${commande}"
        sed '/^[ \t]*$/d;s/#.*//g' $f | strings | sed '1d' > ${basename_f}_nocomments.tsv
        awk 'BEGIN{FS=OFS="\t"; antecedent=""}{
                if( antecedent == "" || antecedent!=$1":"$2"-"$3 ){ 
                    print $1,$2,$3,$NF,$8; 
                    antecedent=$1":"$2"-"$3
                }
            }' ${basename_f}_nocomments.tsv  > ${basename_f}_foldchange_nocomments.bed
        # if exist get file with only summits
        flag_summit=$( grep abs_summit ${f} )
        if [[ ! -z $flag_summit ]]
        then     
            awk 'BEGIN{FS=OFS="\t"}{start=$2; $2=$5; $5=start; $3=$2+1; print $0}' ${basename_f}_nocomments.tsv > ${basename_f}_nocomments_summits.tsv;
            awk 'BEGIN{FS=OFS="\t"; antecedent=""}{
                if( antecedent == "" || antecedent!=$1":"$2"-"$3 ){ 
                    print $1,$2,$3,$NF,$8; 
                    antecedent=$1":"$2"-"$3
                }                  
            }' ${basename_f}_nocomments_summits.tsv | uniq > ${basename_f}_foldchange_nocomments_summits.bed
	else
	    awk 'BEGIN{FS=OFS="\t"}{printf "%s\t%.0f\t%.0f\n", $1,($2+$3)/2,(($2+$3)/2)+1}' ${basename_f}_nocomments.tsv | uniq > ${basename_f}_foldchange_nocomments_center.bed
        fi      
        # if exist get file with summits +/- 100bp
        flag_summit=$( grep abs_summit ${f} )
        if [[ ! -z $flag_summit ]]
        then     
            awk 'BEGIN{FS=OFS="\t"}{start=$2; $2=$5-'${SIZ}'; $3=$5+'${SIZ}'; $5=start; if( $2 < 0 ){ $2=0 }; print $0}' ${basename_f}_nocomments.tsv > ${basename_f}_nocomments_summits${SIZ}.tsv;
            awk 'BEGIN{FS=OFS="\t"; antecedent=""}{
                if( antecedent == "" || antecedent!=$1":"$2"-"$3 ){ 
                    print $1,$2,$3,$NF,$8; 
                    antecedent=$1":"$2"-"$3
                }                  
            }' ${basename_f}_nocomments_summits${SIZ}.tsv | uniq > ${basename_f}_foldchange_nocomments_summits${SIZ}.bed
        fi
    done         
    cd ${DIR_ORI}
}
    
# ==================
# MAIN
# ==================

argv=${@}
echo "[LOG] Your command line: ${argv}" >&2
flag=$( echo ${0} | grep macs2tsv.sh )
if [[ ! -z $flag ]]
then
    macs2tsv ${argv}
fi

