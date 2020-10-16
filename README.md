# m6A_FTO_cancer

The scripts run without any installation as long as the required dependencies and languages are properly installed.
(Tested on Ubuntu 16.04 linux system)

1° m6A-seq analysis
-------------------
**0) prepare input files and set-up the tool**

- get bam files and bedgraph for IP and INPUT using classifical preprocessing tools (see Jeschke at al, 2020)
- set up "dependencies.sh": Update with the appropriate paths and source it on a bash terminal

**1) generate peaks and expected-peaks**

- run "protocol_IPseq-MACS2.sh"

bash tool to generate and annotate peak lists and expected-peaks list.

expected output:
ANNOTATION,  MACS2,  MACS2TSV and  RE-ANNOTATE folders containing respectively peak annotations, MACS2 results, re-sized peaks (tsv and bed formats) and improved annotation with illustrations (pdf format).

usage:
```
bash PeakCalling/bash/protocol_IPseq-MACS2.sh -I=<IP sorted bam file> -C=<INPUT sorted bam file> -o=<output directory> -a=<annotation refFlat format>
```
**2) Differential peaks analysis**

- run "protocol_IPseq-Diff1-vs-1.sh"
Compare Case and Control peak list to prepare for differential. See help
usage:
```
bash AnnotationAnalysis/bash/protocol_IPseq-Diff1-vs-1.sh -h
```

- run "Diff_1vs1_RIPfromBdg.R"
Compute differential using overlap peak list and bedgraphs. See help
usage:
```
bash AnnotationAnalysis/R/Diff_1vs1_RIPfromBdg.R -h
```
- run "Diff_1vs1_RIPfromBdg.R"
Visualise results. See help
usage:
```
bash AnnotationAnalysis/R/Diff_1vs1_RIPort.R -h
```

**3) Motif analysis**
- run "Motif_AMeme.sh"
Identify motif in peak list or evaluate the significance of a known motif (using Meme-Chip or AME respectively). See help
usage:
```
bash AnnotationAnalysis/bash/Motif_AMeme.sh -h
```

2° Public Cohorts Analysis
-----------------------------------------------------------------------
**1) Boxplot_FTO_GDC.R**
Ready-to-use code to reproduce Jeschke et al results on GDC

expected output:
an illutration file (pdf format) and p-values (tsv format) for FTO differential analysis Cancer vs Normal and EMT high versus EMT low.
Note that the code internally compute EMT metascore using ROKAVEC et al genes.

usage:
```
./Rscript Boxplot_FTO_GDC.R
```
**2) Boxplot_FTO_prostateGDS1439.R**
Ready-to-use code to reproduce Jeschke et al results on GDS1439

expected output:
an illutration file (pdf format) for FTO lvels in 'Benign', 'Localized' and 'Metastatic' tumours.

usage:
```
./Rscript Boxplot_FTO_prostateGDS1439.R
