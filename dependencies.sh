############################################################
# Description: List of tools/scripts/executable
# Date: 2016/01
############################################################

# ======================
# Annotation
# ======================

export MACS2="<path to macs2 executable>"
export MACS2TSV=${THESE_TOOLS}/"PeakCalling/bash/macs2tsv.sh"
export ANNOTE="<path to python2.7 executable> "${THESE_TOOLS}/"AnnotationAnalysis/mtools/annotate-peaks"
export ANNOTE_BED="<path to python2.7 executable> "${THESE_TOOLS}/"AnnotationAnalysis/mtools/annotate-peaks-bed"
export ANNOT_MODIF=${THESE_TOOLS}/"PeakCalling/R/Annot_Modifier.R"
export IPMACS2=${THESE_TOOLS}/"PeakCalling/bash/protocol_IPseq-MACS2.sh"

# ======================
# Differential
# ======================

export COMPARE_BED=${THESE_TOOLS}/"PeakCalling/python/compare_bed.py"
export DUPLIC=${THESE_TOOLS}/"Utils/R/duplic_lines.R"
export DIFF_RIPORT_1_1=${THESE_TOOLS}/"PeakCalling/R/Diff_1vs1_RIPort.R"
export BEDTOOLS="<path to bedtools executable>"

# =====================
# Other Analyses
# =====================

export MEME_PATH="<path to meme-suite executables>"
export ALPHABET_PATH=${THESE_TOOLS}/"Data/Alphabets/"

#source dependenciesMartin.sh
