# m6A_FTO_cancer

The scripts run without any installation as long as the required dependencies and languages are properly installed.
(Tested on Ubuntu 16.04 linux system)

1° m6A-seq analysis
-------------------
**0) prepare input files and set-up the tool**

0.1) Get bam files and bedgraphs for IP and INPUT using classifical preprocessing tools (see Jeschke at al, 2020). Examples are in "Data/RNAi/".

0.2) Add the path to your local copy of this git folder in the bash variable "THESE_TOOLS" (Add the following line in your.bashrc for a permanent usage or run on your terminal before any usage of the m6A-seq analyses tools) 
```
export THESE_TOOLS=<path to your local git folder>
```
0.3) Set up "dependencies.sh": Update with the appropriate paths

**1) generate peaks and expected-peaks**

1.1) Run "protocol_IPseq-MACS2.sh"

Description:
bash tool to generate and annotate peak lists and expected-peaks lists using MACS2.

Expected output:
ANNOTATION,  MACS2,  MACS2TSV and  RE-ANNOTATE folders containing respectively peak annotations, MACS2 results, re-sized peaks (tsv and bed formats) and improved annotation with illustrations (pdf format). The peaks lists are called "IPseq-t\<IP bam filename\>-c\<INPUT bam filename\>_peaks_foldchange_nocomments_summits100.bed" (while expected peaks lists are called "IPseq-e\<INPUT bam filename\>_peaks_foldchange_nocomments_summits100.bed" [see Material and Methods of Jeschke et al, 2020 for more details])

Typical usage (for a description of all the parameters please use -h) (Expected runtime: less than 1 day):
```
./PeakCalling/bash/protocol_IPseq-MACS2.sh -I=<IP sorted bam> -C=<INPUT sorted bam> -o=<Output directory> -a=<annotation RefFlat file>
```
Example (note the -es parameter is used here to set-up the extsize MACS2 parameter because of the small size to the example, with normal size bam, extsize is determined automatically by MACS2 [see MACS2 documentation for more information]) (Expected runtime: few minutes):

```
./PeakCalling/bash/protocol_IPseq-MACS2.sh -I=${THESE_TOOLS}/Data/RNAi/RNAiFTO_m6A_chr2_50M_sort.bam -C=${THESE_TOOLS}/Data/RNAi/RNAiFTO_INP_chr2_50M_sort.bam -o=${THESE_TOOLS}/Data/RNAi -a=${THESE_TOOLS}/Data/Annotations/hg19_RefSeq_chr2.refFlat -es=230
./PeakCalling/bash/protocol_IPseq-MACS2.sh -I=${THESE_TOOLS}/Data/RNAi/RNAiCtrl_m6A_chr2_50M_sort.bam -C=${THESE_TOOLS}/Data/RNAi/RNAiCtrl_INP_chr2_50M_sort.bam -o=${THESE_TOOLS}/Data/RNAi -a=${THESE_TOOLS}/Data/Annotations/hg19_RefSeq_chr2.refFlat -es=230

```
Dependencies:
- awk (version 4.1.3)
- MACS2 (version v2.1.0.20150731)
- the provided python code "annotate-peaks.py" and its dependencies:
  - python 2 (version 2.7.12)
  - python modules: argparse (version 1.2.1)
- the provided R code "Annot_Modifier.R" and its dependencies:
  - R (version 3.5.1)

**2) Differential peaks analysis**

2.1) Run "protocol_IPseq-Diff1-vs-1.sh"

Description:
Compare Case and Control peak lists to prepare for differential. "ComparPeaks" cotaining merge peak list from Case and Control on which Differential analysis will be done. Depending on the FLAG_ANNOT_UNION parameter in the configuration file annotation of this file merged list will be done (=1) or not (=0).

Expected output:
"merged.tsv" file, a bed file merging case and control peak list (and eventually "Reannotation.txt" file if FLAG_ANNOT_UNION=1)

Typical usage (for a description of all the parameters please use -h) (Expected runtime: few minutes):
```
./AnnotationAnalysis/bash/protocol_IPseq_Diff1vs1.sh -c=<Configuration file> -pF=<peaks folder> -p1=<peak in Case condition> -p0=<peak in Control condition> -o=<Output directory> -a=<annotation RefFlat file>
```
Example (Note -pF -p1 -p0 are generated at step 1) (Expected runtime: few minutes):
```
./AnnotationAnalysis/bash/protocol_IPseq_Diff1vs1.sh -c=${THESE_TOOLS}/AnnotationAnalysis/Diff1vs1_config.sh -pF=${THESE_TOOLS}/Data/RNAi/MACS2TSV/ -p1=IP-seq-tRNAiFTO_m6A_chr2_50M_sort.bam-cRNAiFTO_INP_chr2_50M_sort.bam_peaks_foldchange_nocomments_summits100.bed -p0=IP-seq-tRNAiCtrl_m6A_chr2_50M_sort.bam-cRNAiCtrl_INP_chr2_50M_sort.bam_peaks_foldchange_nocomments_summits100.bed -o=Data/RNAi -a=Data/Annotations/hg19_RefSeq_chr2.refFlat
```

Dependencies:
- awk (version 4.1.3)
- sed (version 4.2.2)
- the provided python code "compare_bed.py" and its dependencies:
  - python 2 (version 2.7.12)
  - python modules: pybedtools (version 0.7.10), numpy (version 1.11.0), pandas (version 0.17.1)
  - bedtools (version 2.25.0)
- the provided python code "annotate-peaks.py" and its dependencies:
  - python modules: argparse (version 1.2.1)
- the provided R code "Annot_Modifier.R" and its dependencies:
  - R (version 3.5.1)

2.2) Run "Diff_1vs1_RIPfromBdg.R"

Description:
Compute the actual differential using overlap peaks list and bedgraphs. Note that three options are possible, we strongly recommend to use option 2 (used in Jeschke et al, 2020).
Expected output:
"AlternOption2.txt" file, a text file computing Log-odd-Ratio (LOR) differential value between case and control for each peaks and "annot.txt" file showing annotation of the differenial peaks.


Typical usage (Expected runtime: few hours):
```
./AnnotationAnalysis/bash/protocol_IPseq_Diff1vs1.sh -c=${THESE_TOOLS}/AnnotationAnalysis/Diff1vs1_config.sh -pF=${THESE_TOOLS}/Data/RNAi/MACS2TSV/ -p1=IP-seq-tRNAiFTO_m6A_chr2_50M_sort.bam-cRNAiFTO_INP_chr2_50M_sort.bam_peaks_foldchange_nocomments_summits100.bed -p0=IP-seq-tRNAiCtrl_m6A_chr2_50M_sort.bam-cRNAiCtrl_INP_chr2_50M_sort.bam_peaks_foldchange_nocomments_summits100.bed -o=Data/RNAi -a=Data/Annotations/hg19_RefSeq_chr2.refFlat

```
Example (Expected runtime: few minutes):
```
./AnnotationAnalysis/R/Diff_1vs1_RIPfromBdg.R -i=${THESE_TOOLS}/Data/RNAi/ComparPeaks/RNAiFTOm6Achr250M_VS_RNAiCtrlm6Achr250M_mergebed.tsv -o=${THESE_TOOLS}/Data/RNAi/RNAi -a=${THESE_TOOLS}/Data/RNAi/ComparPeaks/RNAiFTOm6Achr250M_VS_RNAiCtrlm6Achr250M_Reannotation.txt -bF=Data/RNAi/ -bL1=RNAiFTO_m6A_chr2_50M.bedgraph -bL0=RNAiCtrl_m6A_chr2_50M.bedgraph -bI1=RNAiFTO_INP_chr2_50M.bedgraph -bI0=RNAiCtrl_INP_chr2_50M.bedgraph -opt=2 -gS=2897310462 -gtf=Data/Annotations/hg19_RefSeq_chr2.gtf
```
Dependencies:
- R (version 3.5.1)
- bedtools (version 2.25.0)
- awk (version 4.1.3)

2° Public Cohorts Analysis
-----------------------------------------------------------------------
**2.1) Boxplot_FTO_GDC.R**

Descripiton:
Ready-to-use code to reproduce Jeschke et al, 2020 results on GDC

Expected output:
an illutration file (pdf format) and p-values (tsv format) for FTO differential analysis Cancer vs Normal and EMT high versus EMT low.
Note that the code internally compute EMT metascore using ROKAVEC et al genes.

Usage (Expected runtime: few minutes):
```
./Rscript Boxplot_FTO_GDC.R
```

Dependencies:
- R (version 3.5.1)


**2.2) Boxplot_FTO_prostateGDS1439.R**

Descripiton:
Ready-to-use code to reproduce Jeschke et al results on GDS1439

Expected output:
an illutration file (pdf format) for FTO lvels in 'Benign', 'Localized' and 'Metastatic' tumours.

Usage (Expected runtime: few minutes):
```
./Rscript Boxplot_FTO_prostateGDS1439.R
```

Dependencies:
- R (version 3.5.1)
