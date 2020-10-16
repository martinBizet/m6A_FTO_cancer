############################################################
# Description: List of tools/scripts/executable
# Date: 2016/01
############################################################

export COCKTAIL="<path to here>"

# ======================
# Annotation
# ======================

export MACS2="<path to macs2 executable>"
export MACS2TSV="PeakCalling/bash/macs2tsv.sh"
export ANNOTE="<path to python2.7 executable> AnnotationAnalysis/mtools/annotate-peaks"
export ANNOTE_BED="<path to python2.7 executable> AnnotationAnalysis/mtools/annotate-peaks-bed"
export ANNOT_MODIF="PeakCalling/R/Annot_Modifier.R"
export IPMACS2="PeakCalling/bash/protocol_IPseq-MACS2.sh"

# ======================
# Differential
# ======================

export COMPARE_BED="PeakCalling/python/compare_bed.py"
export DUPLIC="Utils/R/duplic_lines.R"
export DIFF_RIPORT_1_1="PeakCalling/R/Diff_1vs1_RIPort.R"
export BEDTOOLS="<path to bedtools executable>"
export FEATURE_COUNTS="/illumina/runs/Tools/bin/featureCounts"

# =====================
# Other Analyses
# =====================

export MEME_PATH="<path to meme-suite executables>"
export ALPHABET_PATH="Data/Alphabets/"