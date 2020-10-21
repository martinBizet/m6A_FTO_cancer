#!/usr/bin/env Rscript

argu <- commandArgs(TRUE)

# ================================================
# Title: Differential of RIP using bedgraph
# Author: Martin BIZET
# Date: 2017/10
# ================================================


###-------------------
##	ARGUMENTS & HELP
#---------------------

##Defaults
opt0 <- FALSE; opt1 <- TRUE; opt2 <- TRUE;
ggtf <- NULL; gannot <- NULL;
geno5Size <- NA
peaks <- NULL; output <- NULL; pannot <- NULL;
bF <- ''
bL1 <- bL0 <- bI1 <- bI0 <- NULL;
form <- "auto";
nam <- NULL;
cid <- "name";
nthreads <- 4;
print.help <- FALSE

##Treat argu
argu <- strsplit(argu, split='=')
argu.1 <- argu[which(sapply(argu, length) == 1)]
if(length(argu.1) > 0) {
	for(el in argu.1) {
		if(el[1] == '-help') { print.help <- TRUE }
	}
}
if(print.help) {
	cat('Diff_1vs1_frombedgraph computes the differential between peaks from 2 conditions using IPs and INPUT .bedgraph files and peaks list.\n
		\t-i=|-inputFile=  name of the peak .bed file containing the union of peaks region from the two conditions (output from peak_union.py tool). [mandatory]\n
		\t-a=|-annot=  name of the reannotated peak file (output from Annot_Modifier) [mandatory (except if opt is 0)]\n
		\t-o=|-outputFile=  name of the differential table file (output from this tool). [mandatory]\n
		\t-bF=|-bedgraphFolder=  folder containing the .bedgraph files (2 IP and 2 INPUT) [default is the root]\n
		\t-bL1=|-bdgIP1=  name of the .bedgraph file from IP condition test [mandatory]\n
		\t-bL0=|-bdgIP0=  name of the .bedgraph file from IP condition control [mandatory]\n
		\t-bI1=|-bdgINPUT1=  name of the .bedgraph file from INPUT condition test [mandatory]\n
		\t-bI0=|-bdgINPUT0=  name of the .bedgraph file from INPUT condition control [mandatory]\n
		\t-opt=|-option= how to account for whole gene expression (several options can be coma-separated). 0=only at peak location; 1=ratio peak/whole gene at IP and Input; 2=max(peak, whole gene) at Input [default=1,2]\n
		\t-gS=|-genomeSize=  genome size value. [mandatory] [use 2897310462 for hg19 and 2620345972 for mm9]\n		
		\t-gtf=|-geneAnnotation=  name of the gene annotation file (in .gtf format) [mandatory (except if opt is 0)]\n
		\t-n=|-names=  comma separated names of IP1, IP0, INPUT1, INPUT0 in this order. [default=bL1, bL0, bI1, bI0 names]\n
		\t-help  print this help\n
		Example\n
		-------\n
		Rscript Diff_1vs1_RIPfrombdg.R -i=</path/to/union/peaks/xls/>union_peaks.txt -a=</path/to/union/peaks/xls/>union_peaks.txt-o=</path/to/output/differential/table/>diff_table.csv -bF=</path/to/bedgraph> -bL1=<bL1.bedgraph> -bL0=<bL0.bedgraph> -bI1=<bI1.bedgraph> -bI0=<bI0.bedgraph> -opt=1,2 -gS=2897310462 -gtf=</path/to/annotation>geneAnnotation.gtf-n=TKO,WT,ITKO,IWT\n')
} else {
	argu <- argu[which(sapply(argu, length) > 1)]
	argu.keys <- sapply(argu, '[[', 1)
	argu.vals <- sapply(argu, '[[', 2)
	if(length(intersect(c('-i', '-inputFile'), argu.keys)) > 0) {
		peakfile <- argu.vals[which(is.element(argu.keys, c('-i', '-inputFile')))]
	} else {
		stop('Please provide an input file (-i or -inputFile argument)')
	}
	if(length(intersect(c('-a', '-annot'), argu.keys)) > 0) {
		pannot <- argu.vals[which(is.element(argu.keys, c('-a', '-annot')))]
	}
	if(length(intersect(c('-o', '-outputFile'), argu.keys)) > 0) {
		output <- argu.vals[which(is.element(argu.keys, c('-o', '-outputFile')))]
	} else {
		stop('Please provide an output file (-o or -outputFile argument)')
	}
	if(length(intersect(c('-bF', '-bedgraphFolder'), argu.keys)) > 0) {
		bF <- argu.vals[which(is.element(argu.keys, c('-bF', '-bedgraphFolder')))]
	}	
	if(length(intersect(c('-bL1', '-bdgIP1'), argu.keys)) > 0) {
		bL1 <- argu.vals[which(is.element(argu.keys, c('-bL1', '-bdgIP1')))]
	} else {
		stop('Please provide an test IP bedgraph (-bL1 or -bdgIP1 argument)')
	}
	if(length(intersect(c('-bL0', '-bdgIP0'), argu.keys)) > 0) {
		bL0 <- argu.vals[which(is.element(argu.keys, c('-bL0', '-bdgIP0')))]
	} else {
		stop('Please provide an control IP bedgraph (-bL0 or -bdgIP0 argument)')
	}
	if(length(intersect(c('-bI1', '-bdgINPUT1'), argu.keys)) > 0) {
		bI1 <- argu.vals[which(is.element(argu.keys, c('-bI1', '--bdgINPUT1')))]
	} else {
		stop('Please provide an test INPUT bedgraph (-bI1 or -bdgINPUT1 argument)')
	}
	if(length(intersect(c('-bI0', '-bdgINPUT0'), argu.keys)) > 0) {
		bI0 <- argu.vals[which(is.element(argu.keys, c('-bI0', '-bdgINPUT0')))]
	} else {
		stop('Please provide an control INPUT bam (-bI0 or -bdgINPUT0 argument)')
	}
	if(length(intersect(c('-opt', '-option'), argu.keys)) > 0) {
		opt <- argu.vals[which(is.element(argu.keys, c('-opt', '-option')))]
		opt <- strsplit(opt, split=',')[[1]]
		#OPION0: just TPM at peak region on IP and INPUT: do NOT take whole-gene into account
		opt0 <- is.element('0', opt)
		#OPTION1: calculer l'enrichment ratio (peak regions / whole-gene)
		opt1 <- is.element('1', opt)
		#OPTION2: only modify the counting of INPUT: minimum between level at the windows and in the whole gene
		opt2 <- is.element('2', opt)
	}
	if(length(intersect(c('-gS', '-genomeSize'), argu.keys)) > 0) {
		genoSize <- as.numeric(argu.vals[which(is.element(argu.keys, c('-gS', '-genomeSize')))])
	} else {
		stop('Please provide an genome size value (-gS or -genomeSize argument)')
	}
	if(length(intersect(c('-gtf', '-geneAnnotation'), argu.keys)) > 0) {
		ggtf <- argu.vals[which(is.element(argu.keys, c('-gtf', '-geneAnnotation')))]
	} else {
		stop('Please provide an genome size value (-gS or -genomeSize argument)')
	}
	if(length(intersect(c('-n', '-names'), argu.keys)) > 0) {
		nam <- strsplit(argu.vals[which(is.element(argu.keys, c('-n', '-names')))], split=',')[[1]]
		nam <- nam[c(2,1,4,3)] #Matthieu starts by the control so 0 then 1, while I ask to provide 1 then 0...
	} else {
		f <- function(el) {
			#suppress path
			el <- strsplit(el, split='/')[[1]]
			if(length(grep('accepted_hits', el[length(el)])) > 0) {
				el <- el[length(el)-1]
			} else {
				el <- el[length(el)]
			}
			#suppress .bedgraph extension if any
			if(substr(el, nchar(el)-nchar('.bedgraph')+1, nchar(el)) == '.bedgraph') { return(substr(el, 1, nchar(el)-nchar('.bedgraph'))) }
			if(substr(el, nchar(el)-nchar('.bdg')+1, nchar(el)) == '.bdg') { return(substr(el, 1, nchar(el)-nchar('.bdg'))) }
			return(el)
		}
		nam <- paste(c('bL0', 'bL1', 'bI0', 'bI1'), sapply(c(bL0, bL1, bI0, bI1), f), sep='_')
	}

	
	print(paste('Bedgraph file:', file.path(bF, c(bL0, bL1, bI0, bI1))))
	reads <- file.path(bF, c(bL0, bL1, bI0, bI1)) #Matthieu starts by the control so 0 then 1

	###--------------
	##	FUNCTIONS  --
	#----------------

	source(file.path(Sys.getenv("THESE_TOOLS"),'AnnotationAnalysis/mtools/lib.R'))

	count <- function(peaks, bdg){
		tmp1 <- pipe(paste("bedtools intersect -wb -a", bdg, "-b", peaks))
		peaks <- as.data.frame( read.table(tmp1, sep='\t') )
		peaks$width <- peaks$V3 - peaks$V2
		peaks$count <- peaks$width*peaks$V4
		f <- function(v) { sum(as.numeric(v)) }
		counts <- tapply(peaks$count, peaks$V8, f) #sum)
		names(counts) <- levels(peaks$V8)
		# USELES IF NOT INTERESTED BY INTERGENIC then return(counts)
                tmp2 <- pipe(paste("awk '{s+=($3-$2)*$4}END{print s}'", bdg))
                rest <- as.numeric(readLines(tmp2)) - sum(as.numeric(counts))
                res <- c(counts, rest)
                names(res)[length(res)] <- 'rest'
		return(res)
	}
	
	combine.peak.transcript <- function(exon.count, ge.peak, peak.count, conv, tSize, longest) {
		total.read <- sum(peak.count) #with rest it is the real library size
		peak.count <- peak.count[which(names(peak.count) != 'rest')]
		
		#Gene counts
		nam <- sapply(strsplit(names(exon.count), split='|', fixed=TRUE), '[[', 1)
		ge.count <- tapply(exon.count, nam, sum)
		ge.peak <- ge.count[ge.peak]
		ge.peak <- pmax(ge.peak - peak.count, 0) #if the peak is intronic it can be higher than the full exonic transcript...
		ge.peak <- ge.peak[which(!is.na(ge.peak))]
		ge.peak <- ge.peak[order(ge.peak, decreasing=TRUE)]
		ge.peak <- ge.peak[which(!duplicated(names(ge.peak)))]
		ge.count <- ge.count[which(!is.element(names(ge.count), c(names(ge.peak), 'rest')))]
		res <- c(peak.count, ge.peak, ge.count)
		res <- c(res, max(total.read - sum(res), 0)); names(res)[length(res)] <- 'Intergenic';
		return(res)

	}
		
	get.tpm <- function(counts, pSize) {
		############################################################
		#Source : http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
		############################################################
		#read per kilobase
		rpk <- apply(counts, 2, function(x) (x/pSize)/1000) #pSize should be a vector
		#per million scaling factor
		pmsf <- colSums(rpk)/1e6
		#Calculating tpm by dividing rpk by respective pmsf 
		tpm <- t(apply(rpk,1,function(x) x/pmsf))
		colnames(tpm) <- paste(colnames(counts), 'TPM', sep='_')
		return(tpm)
	}

	get.feat <- function(v, idx=1, as.num=FALSE) {
		if( (length(v) == 1) && (is.na(v)) ) { return(NA) }
		categ <- sapply(strsplit(v, split=','), '[[', idx)
		if(as.num) { categ <- as.numeric(categ) }
		categ <- unique(categ)
		return(categ)
	}
	
	###---------
	##  MAIN  --
	#-----------

	outputC=paste(output, "RawCountsAndSizes.txt", sep='-')
	outputA=paste(output, "Peak2Gene.txt", sep='-')
	output0=paste(output, "AlternOption0.txt", sep='-')
	output1=paste(output, "AlternOption1.txt", sep='-')
	output2=paste(output, "AlternOption2.txt", sep='-')

	print('Loading')
	print(peakfile)
	peaks <- read.table(peakfile, header=FALSE, row.names=NULL, sep='\t', as.is=TRUE)
	colnames(peaks) <- c('#Chr', 'Start', 'End', 'name')
	rownames(peaks) <- peaks$name
	peaks$Length <- peaks$End - peaks$Start
	peaks$Strand <- rep('+', nrow(peaks)) #Strand insensitive analysis
	if(!is.null(pannot)) {
		print(pannot)
		pannot <- read.table(pannot, sep='\t', header=TRUE, row.names=NULL, check.names=FALSE, as.is=TRUE, quote='', comment.char='')
		rownames(pannot) <- pannot$name
	}
	
	if(!is.null(ggtf)) {
		print(ggtf)
		gannot <- read.table(ggtf, header=FALSE, row.names=NULL, sep='\t', as.is=TRUE, quote="\'")
		gannot <- gannot[which(gannot[,3] == 'exon'),]
		tmp.tr <- sapply(strsplit(gannot[,9], split=';'), '[[', 2)
		tmp.tr <- sapply(strsplit(tmp.tr, split='\"'), '[[', 2)
		tmp.ge <- sapply(strsplit(gannot[,9], split=';'), '[[', 1)
		tmp.ge <- sapply(strsplit(tmp.ge, split='\"'), '[[', 2)
		gannot <- cbind(gannot[, c(1, 4, 5)], paste(tmp.ge, rownames(gannot), sep='|'), tmp.tr)
		conv <- tmp.ge; names(conv) <- tmp.tr;
		conv <- conv[which(!duplicated(names(conv)))]
		colnames(gannot) <- c('#Chr', 'Start', 'End', 'name', 'tr')
		rownames(gannot) <- gannot$name
		gannot$Length <- gannot$End - gannot$Start
		gannot$Strand <- rep('+', nrow(gannot))
		tSize <- tapply(gannot$Length, tmp.tr, sum)
		conv <- conv[names(tSize)]
		f <- function(v) { return(names(v)[which.max(v)]) }
		longest <- tapply(tSize, conv, f)
		gannot <- gannot[which(is.element(gannot$tr, longest)),]
		gannot <- gannot[-which(colnames(gannot) == 'tr')]
		genefile <- paste(ggtf, "exon.bed", sep='_')
		print(paste('writing:', genefile))
		write.table(gannot, genefile, col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
	}
	
	print('peaks counts')
	counts.ac <- mat.or.vec(nr=nrow(peaks)+1, nc=4); dimnames(counts.ac) <- list(c(as.character(peaks$name), 'rest'), nam);
	bdg.l <- c(bL0, bL1, bI0, bI1)
	#count in peaks
	for(i in 1:ncol(counts.ac)) {
		tmp <- file.path(bF, bdg.l[i]); print(tmp);
		tmp <- count(peakfile, tmp); counts.ac[names(tmp), i] <- tmp;
	}	
	
	if(opt1 || opt2) {
		if(is.null(gannot)) { stop('-gtf must be provided if for option 1 and 2') }
		if(is.null(pannot)) { stop('-a must be provided if for option 1 and 2') }
		print('gene counts')
		counts <- mat.or.vec(nr=nrow(peaks)+length(longest)+1, nc=4); dimnames(counts) <- list(c(as.character(peaks$name), names(longest), 'Intergenic'), nam);
		corresp <- list()
		#get transcript from pannot and extract gene count
		f <- function(el) { return(sapply(strsplit(el, split=','), '[[', 1)) }
		tr <- sapply(strsplit(pannot$annotation, split=';'), f)
		f <- function(v, conv, tSize, longest) { 
			if(length(v) == 1 && is.na(v)) { return(NA) } #intergenic
			ge <- conv[v]
			ge <- ge[which.max(tSize[longest[ge]])] #if several genes the longest is selected
			if(length(ge) == 0) { return(NA) } #problem with version of RefSeq; it is not referred in gtf => consider it as intergenic
			return(ge)
		}
		ge.peak <- sapply(tr, f, conv=conv, tSize=tSize, longest=longest)
		ge.peak[is.na(ge.peak)] <- 'Intergenic'
	        write.table(cbind(peak=rownames(pannot), Gene=ge.peak), outputA, row.names=FALSE, quote=FALSE, sep = "\t")
		#count in genes
		for(i in 1:ncol(counts)) {
			tmp <- file.path(bF, bdg.l[i]); print(tmp);
			tmp <- count(genefile, tmp)
			tmp <- combine.peak.transcript(tmp, ge.peak, counts.ac[,i], conv, tSize, longest)
			counts[names(tmp), i] <- tmp
		}
	}

	pSize <- peaks$End - peaks$Start
	if(opt1 || opt2) {
		size <- c(pSize, tSize[longest], genoSize-sum(pSize))
	} else {
		size <- c(pSize, genoSize-sum(pSize))
		counts <- counts.ac
	}
	write.table(cbind(counts, size=size), outputC, quote=FALSE, sep = "\t")

        print('TPM normalisation')
	tpm <- get.tpm(counts, size)

	#Pseudocount: 1 TPM
	tpm <- tpm+1
	tpm.ac <- tpm[as.character(peaks$name),]
	if(opt1 || opt2) {
		tpm.bd <- tpm[ge.peak,] 
		rownames(tpm.bd) <- rownames(tpm.ac)
	}

	print('Compute differential')
	# compute the differential analysis and write it
	# LOR: log odd ratio (IP / input) / (IP / input))
	# LR: log ratio IP / IP
	# IR: logratio Input / Input
	nam <- paste(nam, 'TPM', sep='_')	
	if(opt0) {
		enrich <- log2(tpm.ac)
		dif.opt0 <- differential.table(enrich, nam[2], nam[1], nam[4], nam[3])
		dif.opt0 <- cbind(peaks, dif.opt0[rownames(peaks),])
		write.table(dif.opt0, output0, quote = F, row.names = F, sep = "\t")	
	}
	
	if(opt1) {
		#OPTION1: calculer l'enrichment score et l'enrichment ratio sur les regions definies par MACS2...
		#compute log2 enrich
		enrich <- log2(tpm.ac/tpm.bd)
		dif.opt1 <- differential.table(enrich, nam[2], nam[1], nam[4], nam[3])
		dif.opt1 <- cbind(peaks, dif.opt1[rownames(peaks),])
		write.table(dif.opt1, output1, quote = F, row.names = F, sep = "\t")
	}
	
	if(opt2) {
		#OPTION2: only modify the counting of INPUT: minimum between level at the windows and in the whole gene
		enrich <- log2(tpm.ac)
		tpm.bd <- log2(tpm.bd)
		enrich[,3] <- pmax(enrich[,3], tpm.bd[,3])
		enrich [,4] <- pmax(enrich[,4], tpm.bd[,4])
		dif.opt2 <- differential.table(enrich, nam[2], nam[1], nam[4], nam[3])
		dif.opt2 <- cbind(peaks, dif.opt2[rownames(peaks),])
		write.table(dif.opt2, output2, quote = F, row.names = F, sep = "\t")	
	}

	###-------------------
	##	CBIND ANNOTATION
	#---------------------

	if(!is.null(pannot)) {
		print('Add annotation')
		pannot <- pannot[peaks$name,]
		colnames(pannot)[ncol(pannot)] <- 'category'
		pannot$gene <- sapply(sapply(strsplit(pannot$annotation, split=';'), get.feat, idx=2), paste, collapse=';')
		if(opt0) {
			comb <- cbind(dif.opt0, abs_summit=pannot$abs_summit, annotation=pannot$annotation, gene=pannot$gene,category=pannot$category)
			write.table(comb, paste(output0, 'annot.txt', sep='_'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
		}
		if(opt1) {
			comb <- cbind(dif.opt1, abs_summit=pannot$abs_summit, annotation=pannot$annotation, gene=pannot$gene,category=pannot$category)
			write.table(comb, paste(output1, 'annot.txt', sep='_'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
		}
		if(opt2) {
			comb <- cbind(dif.opt2, abs_summit=pannot$abs_summit, annotation=pannot$annotation, gene=pannot$gene,category=pannot$category)
			write.table(comb, paste(output2, 'annot.txt', sep='_'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
		}
	}
}
