#!/usr/bin/env Rscript
argu <- commandArgs(TRUE)

###-------------
##	ARGUMENTS & HELP
#---------------

##Defaults
infile <- outfile <- expect <- NULL
leftSta <- rightSta <- leftStop <- rightStop <- leftTTS <- rightTSS <- 0
print.help <- FALSE
draw.annot <- FALSE
draw.promo <- FALSE
draw.TTSout <- FALSE
conflict.mode <- 'Multiple'
def.prior <- c('near_Start', 'near_Stop', 'near_TSS', 'near_TTS', '5UTR', '3UTR', 'CDS', 'ncExon', 'intron', 'ncIntron', 'Intergenic')
extra.annot <- NULL

##Treat argu
argu <- strsplit(argu, split='=')
argu.1 <- argu[which(sapply(argu, length) == 1)]
if(length(argu.1) > 0) {
	for(el in argu.1) {
		if(el[1] == '-help') {
			print.help <- TRUE
		} else if (el[1] == '-draw') {
			draw.annot <- TRUE
			sel.categ <- NULL
		} else if (el[1] == '-includePromo') {
			draw.promo <- TRUE
			def.prior <- c('Promoter', def.prior)
		} else if (el[1] == '-includeTTSout') {
			draw.TTSout <- TRUE
			i <- which(def.prior == 'near_TTS')
			if(is.element('Promoter', def.prior)) {
				tmp <- c('Promoter', 'near_TTS', def.prior[2:(i-1)])
			} else {
				tmp <- c('near_TTS', def.prior[1:(i-1)])
			}
			def.prior <- c(tmp, def.prior[(i+1):length(def.prior)])
		}
	}
}
if(print.help) {
	cat('Annot_Modifier resumes the annotation to one category by peak and generate a one-line-by-peak-annotation version of the annotation (.duplic file). Depending of the parameters it will also refine the annotation and/or draw the repartition.\n
		\t-i|-inputFile=  name of the Annotated peak file to use as input (output from annotate-peaks mtool) [mandatory]\n
		\t-e|-expectFile=  name of the Annotated peak file to use as expected values (output from annotate-peaks mtool on expected peaks) [ignored if not provided]\n
		\t-o|-outputFile=  name of the Re-Annotated peak file (output from this tool) [mandatory]\n
		\t-m|-solvingConflictMode= mode to solve annotation conflicts. If -exF is provided, this mode is used to solve the remaining conflicts. Could be \'Multiple\' or \'Prioritise\'. [default=\'Multiple\']\n
		\t-dp|-defaultPriority= category priority when \'Prioritise\' is used for -m. [default=\'near_Start+near_Stop+near_TSS+near_TTS+5UTR+3UTR+CDS+ncExon+intron+ncIntron+Intergenic\']\n
		\t-exF|externalAnnotationFiles= path to external annotation file used to solve annotation conflicts. The first column (rownames) should be transcript ID (same identifier used in -i file). If omit, no extra file is used\n
		\t-exC|externalAnnotationColumns= comma-separated columnnames to use for solving conflicts. the order of the columnd determine their importance in solving the confilcts. Unused if -exF is omit and mandatory is -exF is provided\n
		\t-exP|externalAnnotationPriority= Priority attributed to each category of a -exC column. Categories are separated by \'+\', priority of each -exC columns are coma-separated (same order then -exC).\n
		\t-leftStart=  left border for the near_Start category which surround the start codon [default=0]\n
		\t-rightStart=  right border for the near_Start category which surround the start codon [default=0]\n
		\t-leftStop=  left border for the near_Stop category which surround the stop codon [default=0]\n
		\t-rightStop=  right border for the near_Stop category which surround the stop codon [default=0]\n
		\t-leftTTS=  left border for the near_TTS category which surround the TTS [default=0]\n
		\t-rightTSS=  right border for the near_TSS category which surround the TSS [default=0]\n
		\t-draw|-draw=  draw a pdf with category repartitions. When no argument or TRUE is provided all categories are drawn. Alternatively a comma-separated list of categories can be provided to restrict the representation. [default=FALSE]\n
		\t-includePromo  Include promoter category in drawing and reannoting (for ChIP analysis only). [default=FALSE]\n
		\t-includeTTSout  Include nearTTS category in drawing and reannoting (for ChIP analysis only). This category overwrite near_TTS inclusively intra-transcript and -leftTTS parameter is ignored [default=FALSE]\n
		\t-help  print this help\n
		Example\n
		-------\n
		Rscript Annot_Modifier.R -i=</path/to/annotated/peaks/>1608-ES-WT-hmC-D_annotatedpeaks.txt -o=</output/dir/>1608-ES-WT-hmC-D_Reannotation.txt -leftStop=-500 -rightStop=500 -draw=5UTR,CDS,ncExon,intron,ncIntron,near_Stop,3UTR')
} else {
	argu <- argu[which(sapply(argu, length) > 1)]
	argu.keys <- sapply(argu, '[[', 1)
	argu.vals <- sapply(argu, '[[', 2)
	if(length(intersect(c('-i', '-inputFile'), argu.keys)) > 0) {
		infile <- argu.vals[which(is.element(argu.keys, c('-i', '-inputFile')))]
	} else {
		stop('Please provide an input file (-i or -inputFile argument)')
	}
	if(length(intersect(c('-e', '-expectFile'), argu.keys)) > 0) {
		expect <- argu.vals[which(is.element(argu.keys, c('-e', '-expectFile')))]
	}
	if(length(intersect(c('-o', '-outputFile'), argu.keys)) > 0) {
		outfile <- argu.vals[which(is.element(argu.keys, c('-o', '-outputFile')))]
	} else {
		stop('Please provide an output file (-o or -outputFile argument)')
	}
	if(length(intersect(c('-m', '-solvingConflictMode'), argu.keys)) > 0) {
		conflict.mode <- argu.vals[which(is.element(argu.keys, c('-m', '-solvingConflictMode')))]
		if(!is.element(conflict.mode, c('Multiple', 'Prioritise'))) { stop('wrong -m argument. Please use \'Multiple\'or \'Prioritise\'') }
	}
	if(length(intersect(c('-dp', '-defaultPriority'), argu.keys)) > 0) {
		def.prior <- argu.vals[which(is.element(argu.keys, c('-dp', '-defaultPriority')))]
		def.prior <- strsplit(def.prior, split='\\+')[[1]]
	}
	if(length(intersect(c('-exF', '-externalAnnotationFiles'), argu.keys)) > 0) { extra.annot <- argu.vals[which(is.element(argu.keys, c('-exF', '-externalAnnotationFiles')))] }
	if(length(intersect(c('-exC', 'externalAnnotationColumns'), argu.keys)) > 0) {
		extra.cols <- argu.vals[which(is.element(argu.keys, c('-exC', 'externalAnnotationColumns')))]
		extra.cols <- strsplit(extra.cols, split=',')[[1]]
	}
	if(length(intersect(c('-exP', 'externalAnnotationPriority'), argu.keys)) > 0) {
		extra.prior <- argu.vals[which(is.element(argu.keys, c('-exP', '-externalAnnotationPriority')))]
		extra.prior <- strsplit(strsplit(extra.prior, split=',')[[1]], split='\\+')
	}
	if(is.element('-leftStart', argu.keys)) { leftSta <- as.numeric(argu.vals[which(argu.keys == '-leftStart')]) }
	if(is.element('-rightStart', argu.keys)) { rightSta <- as.numeric(argu.vals[which(argu.keys == '-rightStart')]) }
	if(is.element('-leftStop', argu.keys)) { leftStop <- as.numeric(argu.vals[which(argu.keys == '-leftStop')]) }
	if(is.element('-rightStop', argu.keys)) { rightStop <- as.numeric(argu.vals[which(argu.keys == '-rightStop')]) }
	if(is.element('-leftTTS', argu.keys)) { leftTTS <- as.numeric(argu.vals[which(argu.keys == '-leftTTS')]) }
	if(is.element('-rightTSS', argu.keys)) { rightTSS <- as.numeric(argu.vals[which(argu.keys == '-rightTSS')]) }
	if(is.element('-draw', argu.keys)) {
		draw.annot <- argu.vals[which(argu.keys == '-draw')]
		if(draw.annot != 'FALSE'){
			if(draw.annot == 'TRUE') {
				sel.categ <- NULL
			} else {
				sel.categ <- strsplit(draw.annot, split=',')[[1]]
			}
			draw.annot <- TRUE			
		}
	}
	
	###-------------
	##	FUNCTIONS
	#---------------

	get.feat <- function(v, idx=1, as.num=FALSE, uniq=TRUE) {
		if( (length(v) == 1) && (is.na(v)) ) { return(NA) }
		categ <- sapply(strsplit(v, split=','), '[[', idx)
		if(as.num) { categ <- as.numeric(categ) }
		if(uniq) { categ <- unique(categ) }
		return(categ)
	}
        ncAware.get.feat <- function(v, idx=1, as.num=FALSE, uniq=TRUE) {
		ncidx <- 3
                if( (length(v) == 1) && (is.na(v)) ) { return(NA) }
                categ <- sapply(strsplit(v, split=','), '[[', idx)
		nc <- is.element( sapply(strsplit(v, split=','), '[[', ncidx), c('ncPromoter', 'ncTTS', 'ncIntron', 'ncExon') )
		categ <- categ[which(!nc)]
                if(as.num) { categ <- as.numeric(categ) }
                if(uniq) { categ <- unique(categ) }
                return(categ)
        }
	modif.annot <- function(v, def.prior, extra.prior, leftSta, rightSta, leftStop, rightStop, rightTSS, leftTTS, TTSout) {
		if( (length(v) == 1) && (is.na(v)) ) { return('Intergenic') }
		modif.one <- function(v, prior, leftSta, rightSta, leftStop, rightStop) {
			categ <- v[3]
			if(categ == 'ncPromoter') { categ='Promoter' } #regroup promoters of coding and non-coding
			if(categ == 'ncTTS') { categ='near_TTS' } #regroup promoters of coding and non-coding
			dists <- as.numeric(v[5:8])
			tmp <- c( (dists[1] >= 0 & dists[1] < rightTSS), (dists[2] > leftTTS), (dists[3] > leftSta & dists[3] < rightSta), (dists[4] > leftStop & dists[4] < rightStop) ) #can be negative in case of Promoters usage
			if( !is.na(tmp[1]) && tmp[1] == TRUE ) {
				categ <- 'near_TSS'
			} else if( ! TTSout && !is.na(tmp[2]) && tmp[2] == TRUE ) { #in case of TTSout near_TTS category already exist and should not be recomputed using leftTTS
				categ <- 'near_TTS'
			} else if( is.element(categ, c('ncIntron', 'ncExon')) ) {
			    if(sum(tmp[1:2], na.rm=TRUE) > 1) {
					if(is.null(prior)) {
						categ <- 'Multiple'
					} else {
						categ <- categ[which.min(factor(categ, levels=prior))]
					}
				}
			} else {
			    if(sum(tmp, na.rm=TRUE) > 1) {
					if(is.null(prior)) {
						categ <- 'Multiple'
					} else {
						categ <- categ[which.min(factor(categ, levels=prior))]
					}
				} else if( !is.na(tmp[3]) && tmp[3] == TRUE ) {
					categ <- 'near_Start'
				} else if( !is.na(tmp[4]) && tmp[4] == TRUE ) {
					categ <- 'near_Stop'
				}
			}
			return(c(categ, v[1]))
		}
		categ.tab <- t(sapply(strsplit(v, split=','), modif.one, prior=def.prior, leftSta=leftSta, rightSta=rightSta, leftStop=leftStop, rightStop=rightStop))
		categ <- unique(categ.tab[,1])
		if(length(categ) > 1) {
			if(is.null(def.prior)) { return('Multiple') } #solve conflict using Multiple
			if(is.null(extra.prior)) { return(categ[which.min(factor(categ, levels=def.prior))]) } #solve conflict using only def.prior
			for(i in 1:length(extra.prior)) { #solve conflict using the extra.prior table...
				if(length(categ) > 1) {
					tmp <- as.numeric(extra.prior[categ.tab[,2], i])
					if(sum(is.na(tmp)) < length(tmp)) { #if a column is full of NA everything is deleted...
						categ.tab <- categ.tab[which(tmp == min(tmp, na.rm=TRUE)),,drop=FALSE]
						categ <- unique(categ.tab[,1])
					}
				}
			}
			if(length(categ) > 1) { return(categ[which.min(factor(categ, levels=def.prior))]) } #... and solve the remaining ones (if any) using def.prior
		}
		return(categ)
	}
	reNorm <- function(v) {
		dTSS <- v[1]
		dTTS <- v[2]
		dSta <- v[3]
		dSto <- v[4]
		if(sum(is.na(v) > 0)) { return(NA) }
		if(sum(abs(v[3:4]) == 0)) { return(NA) }
		if(dSto > 0) { 
			#3UTR
			dReNorm <- 2+(dSto/(dSto-dTTS+1))
		} else if(dSta > 0) { 
			#CDS
			dReNorm <- 1+(dSta/(dSta-dSto+1))
		} else { 
			#5UTR
			dReNorm <- 0+(dTSS/(dTSS-dSta+1))
		}
		return(dReNorm/3)
	}
	Norm <- function(v) {
		dTSS <- v[1]
		dTTS <- v[2]
		if(sum(is.na(v) > 0)) { return(NA) }
		if(sum(abs(v) == 0)) { return(NA) }
		return(dTSS/(dTSS-dTTS+1))
	}
	transform.annot <- function(annot.unl) {
		names(annot.unl) <- paste(names(annot.unl), '.nGum5eVr8o.', sep='')
		annot.unl <- unlist(annot.unl)
		annot.unl <- sapply(annot.unl, strsplit, split=',')
		nb.feat <- sapply(annot.unl, length)
		for(i in which(nb.feat == 1)) { annot.unl[[i]] <- rep(NA, max(nb.feat)) }
		annot.unl <- data.frame(t(as.matrix(data.frame(annot.unl, stringsAsFactors=FALSE, check.names=FALSE))), stringsAsFactors=FALSE, check.names=FALSE)
		annot.unl$PeakID <- sapply(strsplit(rownames(annot.unl), split='.nGum5eVr8o.'), '[[', 1)
                return(annot.unl)
	}
	prepare <- function(annot, annot.splitted, annot.unl, borderTSS, borderTTS, borderCDS, bins, binsNorm, stato, dex, colo, draw.promo=FALSE, draw.TTSout=FALSE) {
		dNorm <- unlist(sapply(annot.splitted, get.feat, idx=4, as.num=TRUE))
		dTSS <- unlist(sapply(annot.splitted, get.feat, idx=5, as.num=TRUE)); dTTS <- unlist(sapply(annot.splitted, get.feat, idx=6, as.num=TRUE));
		if(draw.promo) {
			dTSS.cut <- dTSS[which(abs(dTSS) < borderTSS)]
			if(draw.TTSout) {
				dTTS.cut <- dTTS[which(abs(dTTS) < abs(borderTTS))]
			} else { #There are case when the gene is so small that promoter go outside TTS, I want to ignore those case if TSSout is not TRUE
				dTTS.cut <- dTTS[which(dTTS > borderTTS & dTTS <= 0)]
			}
		} else if(draw.TTSout) { #if user want TTSout but not promoter... I think it will be rare, it is just in case
			dTSS.cut <- dTSS[which(dTSS < borderTSS)]
			dTTS.cut <- dTTS[which(abs(dTTS) < abs(borderTTS))]
		} else {
			dTSS.cut <- dTSS[which(dTSS < borderTSS)]
			dTTS.cut <- dTTS[which(dTTS > borderTTS)]
		}
		dReNorm <- NULL; dsta.cut=NULL; dstop.cut=NULL;
		if(stasto) {
			dsta <- unlist(sapply(annot.splitted, ncAware.get.feat, idx=7, as.num=TRUE)); dstop <- unlist(sapply(annot.splitted, ncAware.get.feat, idx=8, as.num=TRUE));
			dsta.cut <- dsta[which(abs(dsta) < borderCDS)]; dstop.cut <- dstop[which(abs(dstop) < borderCDS)];
			dReNorm <- apply(apply(annot.unl[,c('dist_TSS', 'dist_TTS', 'dist_start', 'dist_stop')], 2, as.numeric), 1, reNorm)
		}
		dexNorm <- NULL; dexReNorm <- NULL;
		if(dex) {
			tmp <- apply(annot.unl[,c('exonic_dist_TSS', 'exonic_dist_TTS', 'exonic_dist_start', 'exonic_dist_stop')], 2, as.numeric)
			dexNorm <- apply(tmp[,1:2], 1, Norm)
			dexReNorm <- apply(tmp, 1, reNorm)
		}
		count.categ <- rep(0, length(colo)) #to deal with empty categories
		names(count.categ) <- names(colo)
		tmp <- table(annot$Reannot)
		count.categ[names(tmp)] <- tmp	
# 		if(is.null(sel.categ)) { #sel.categ is never NULL anymore at this step (define at the beginning of the "DRAW" section)
# 			sel.categ <- names(colo) #use names(colo) order by default; empty categ not drawn by default;[which(is.element(names(colo), names(tmp)))] 
# 		} else if(sum(is.element(sel.categ, names(count.categ))) < length(sel.categ)) {
		if(sum(is.element(sel.categ, names(count.categ))) < length(sel.categ)) {
			stop('one of the category is not available. Please verify -draw argument')
		}
		count.categ <- count.categ[sel.categ]
		colo <- colo[sel.categ]
		return(list(count.categ=count.categ, colo=colo, dNorm=dNorm, dReNorm=dReNorm, dexNorm=dexNorm, dexReNorm=dexReNorm, dTSS.cut=dTSS.cut, dTTS.cut=dTTS.cut, dsta.cut=dsta.cut, dstop.cut=dstop.cut))
	}
	plot.dens <- function(d, d.e=NULL, stasto, main) {
		if(sum(is.na(d)) == length(d)) { 
			dens <- list(x=c(0,100), y=c(0,0))
		} else {
			dens <- density(d*100, bw=1, na.rm=TRUE)
		}
		if(!is.null(d.e)) {
			if(sum(is.na(d.e)) == length(d.e)) {
				dens.e <- list(x=c(0,100), y=c(0,0))        
			} else {
				dens.e <- density(d.e*100, bw=1, na.rm=TRUE)
			}
			ymax <- max(c(0.0001, dens$y, dens.e$y))
			plot(dens.e, xlim=c(0, 100), ylim=c(0, ymax*1.05), main=main, type='l', lwd=2, col='grey', xlab='Relative position in transcript', ylab='Density')
			lines(dens, lwd=2, col='blue')
		} else {
			ymax <- max(c(0.0001, dens$y))
			plot(dens, xlim=c(0, 100), ylim=c(0, ymax*1.05), main=main, type='l', lwd=2, col='blue', xlab='Relative position in transcript', ylab='Density')				
		}
		if(stasto) { 
			borders <- 100*(1:2)/3
			abline(v=borders, lty=2)
			text(x=c(borders, 100) - borders[1]/2, y=ymax, labels=c('5UTR', 'CDS', '3UTR'))
		}
	}
	###-------------
	##	REANNOTATE
	#---------------
	
	annot <- read.table(infile, sep='\t', header=TRUE, row.names=NULL, check.names=FALSE, as.is=TRUE, quote='', comment.char='')
	if(conflict.mode == 'Multiple') {
		def.prior <- NULL
	} else if(!is.null(extra.annot)) {
		extra.annot <- read.table(extra.annot, sep='\t', header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE, quote='', comment.char='')
		extra.annot <- extra.annot[,extra.cols]
		for(i in 1:length(extra.prior)) { extra.annot[,i] <- factor(extra.annot[,i], levels=extra.prior[[i]]) }
	}
	rownames(annot) <- annot$name
	annot.splitted <- strsplit(annot$annotation, split=';'); names(annot.splitted) <- rownames(annot);
	annot$Reannot <- sapply(annot.splitted, modif.annot, def.prior=def.prior, extra.prior=extra.annot, leftSta=leftSta, rightSta=rightSta, leftStop=leftStop, rightStop=rightStop, leftTTS=leftTTS, rightTSS=rightTSS, TTSout=draw.TTSout)
	annot.df <- data.frame(chrom=annot$chr, chromStart=annot$start, chromEnd=annot$end, name=rownames(annot), annot[, 4:ncol(annot)], stringsAsFactors=FALSE, check.names=FALSE)
	colnames(annot.df)[1] <- paste('#', colnames(annot.df)[1], sep='')
	write.table(annot.df, outfile, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
	with.expect <- !is.null(expect)
	if(with.expect) {
		expect <- read.table(expect, sep='\t', header=TRUE, row.names=NULL, check.names=FALSE, as.is=TRUE, quote='', comment.char='')
		rownames(expect) <- expect$name
		expect.splitted <- strsplit(expect$annotation, split=';'); names(expect.splitted) <- rownames(expect);
		expect$Reannot <- sapply(expect.splitted, modif.annot, def.prior=def.prior, extra.prior=extra.annot, leftSta=leftSta, rightSta=rightSta, leftStop=leftStop, rightStop=rightStop, leftTTS=leftTTS, rightTSS=rightTSS)
	}
	
	###-------------
	##  REORGANISE
	#---------------
	
	annot.unl <- transform.annot(annot.splitted)
	tmp <- c('transcript_ID', 'HugoName', 'Region', 'Normalised_Position', 'dist_TSS', 'dist_TTS', 'dist_start', 'dist_stop')
	stasto <- TRUE
	dex <- FALSE
	if(length(tmp)+1 > ncol(annot.unl)) {
		print('[WARNING] distance to start/stop codon not provided !!')
		tmp <- tmp[1:(ncol(annot.unl)-1)]
		stasto <- FALSE
	} else if (length(tmp)+1 < ncol(annot.unl)) {
		tmp <- c('transcript_ID', 'HugoName', 'Region', 'Normalised_Position', 'dist_TSS', 'dist_TTS', 'dist_start', 'dist_stop', 'exonic_dist_TSS', 'exonic_dist_TTS', 'exonic_dist_start', 'exonic_dist_stop')
		dex <- TRUE
	}
	colnames(annot.unl) <- c(tmp, 'PeakID')
	if(with.expect) {
		expect.unl <- transform.annot(expect.splitted)
		if(ncol(annot.unl) != ncol(expect.unl)) { stop('-i and -e should have the same structure') }
		colnames(expect.unl) <- colnames(annot.unl)
	}	
	annot.reord <- annot.df[annot.unl$PeakID,]
	annot.reord <- cbind(annot.reord[, 1:(ncol(annot.reord)-2)], annot.unl[, 1:(ncol(annot.unl)-1)], Reannot=annot.reord[, ncol(annot.reord)])
	write.table(annot.reord, paste(outfile, 'duplic', sep='.'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
	if(draw.annot) {
		
		###-------------
		##	DRAW
		#---------------
		
		##Prepare
		borderCDS <- 3000; bins <- 25; binsNorm <- 0.005;
		if(draw.promo) {
			borderTSS <- 3000
		} else {
			borderTSS <- 5000
		}
		if(draw.TTSout) {
			borderTTS <- 3000
		} else {
			borderTTS <- -5000
		}
		colo <- c('royalblue4', 'darkorange', 'springgreen4', 'cornflowerblue', 'cyan', 'darkmagenta', 'mediumorchid1', 'red3', 'yellow', 'pink', 'gray', 'grey60')
		names(colo) <- c('near_TSS', '5UTR', 'near_Start', 'CDS', 'ncExon', 'intron', 'ncIntron', 'near_Stop', '3UTR', 'near_TTS', 'Multiple', 'Intergenic')
		if(draw.promo) { colo <- c('purple', colo); names(colo)[1] <- 'Promoter'; }
		if(is.null(sel.categ)) { #take all categ except the ones that are impossible to achieve because their sizes == 0
			if(rightTSS == 0) { colo <- colo[-1*which(names(colo) == 'near_TSS')] }		
			if(rightSta-leftSta == 0) { colo <- colo[-1*which(names(colo) == 'near_Start')] }
			if(rightStop-leftStop == 0) { colo <- colo[-1*which(names(colo) == 'near_Stop')] }
			if(leftTTS == 0 && ! draw.TTSout) { colo <- colo[-1*which(names(colo) == 'near_TTS')] }
			sel.categ <- names(colo)
		}
		#annot
		prep <- prepare(annot=annot, annot.splitted=annot.splitted, annot.unl=annot.unl, borderTSS=borderTSS, borderTTS=borderTTS, borderCDS=borderCDS, bins=bins, binsNorm=binsNorm, stato=stato, dex=dex, colo=colo, draw.promo=draw.promo, draw.TTSout=draw.TTSout)
		count.categ <- prep$count.categ; colo <- prep$colo; dNorm <- prep$dNorm; dReNorm <- prep$dReNorm; dexNorm <- prep$dNorm; dexReNorm <- prep$dexReNorm;
		dTSS.cut <- prep$dTSS.cut; dTTS.cut <- prep$dTTS.cut; dsta.cut <- prep$dsta.cut; dstop.cut <- prep$dstop.cut;
		dNorm.e <- dReNorm.e <- dexNorm.e <- dexReNorm.e <- NULL
		#expect
		if(with.expect) {
			colo.e <- rep('grey', 12)
			names(colo.e) <- c('near_TSS', '5UTR', 'near_Start', 'CDS', 'ncExon', 'intron', 'ncIntron', 'near_Stop', '3UTR', 'near_TTS', 'Multiple', 'Intergenic')
			if(draw.promo) { colo <- c('grey', colo); names(colo)[1] <- 'Promoter'; }
			prep <- prepare(annot=expect, annot.splitted=expect.splitted, annot.unl=expect.unl, borderTSS=borderTSS, borderTTS=borderTTS, borderCDS=borderCDS, bins=bins, binsNorm=binsNorm, stato=stato, dex=dex, colo=colo.e, draw.promo=draw.promo, draw.TTSout=draw.TTSout)
			expect.categ <- prep$count.categ; colo.e <- prep$colo; dNorm.e <- prep$dNorm; dReNorm.e <- prep$dReNorm; dexNorm.e <- prep$dexNorm; dexReNorm.e <- prep$dexReNorm;
		}
		pdf(paste(outfile, 'pdf', sep='.'))
			##Barplot catego
			opar <- par(mar=c(8.1,4.1,4.1,1.1), las=2, lwd=2)
				cnam <- names(count.categ)
				if(is.element('near_TSS', cnam)) { cnam[which(cnam == 'near_TSS')] <- 'near_TSS (cod+ncRNA)' }
				if(is.element('near_TTS', cnam)) { cnam[which(cnam == 'near_TTS')] <- 'near_TTS (cod+ncRNA)' }
				if(is.element('ncIntron', cnam)) { cnam[which(cnam == 'ncIntron')] <- 'ncRNA intron' }
				if(is.element('ncExon', cnam)) { cnam[which(cnam == 'ncExon')] <- 'ncRNA exon' }
				if(with.expect) {
					perc <- rbind(observed=count.categ*100/sum(count.categ), expected=expect.categ*100/sum(expect.categ)); colnames(perc) <- cnam;
					barplot(perc, col=rbind(colo, colo.e), beside=TRUE, main='Transcriptomic regions', ylab='percentage of peaks')
					text(x=(1:length(perc))*3-1.5, y=perc, labels=paste(round(perc, 2), '%'), pos=3, cex=0.8)
					tmp <- rbind(observed=count.categ, expected=expect.categ); colnames(tmp) <- cnam;
					barplot(tmp, col=rbind(colo, colo.e), beside=TRUE, main='Transcriptomic regions', ylab='number of peaks')
					text(x=(1:length(tmp))*3-1.5, y=tmp, labels=paste(round(perc, 2), '%'), pos=3, cex=0.8)
				} else {
					perc <- count.categ*100/sum(count.categ); names(perc) <- cnam;
					barplot(perc, col=colo, main='Transcriptomic regions', ylab='percentage of peaks', ylim=c(0, max(perc)*1.1))
					text(x=(1:length(perc))*1.2-0.5, y=perc, labels=paste(round(perc, 2), '%'), pos=3)
					tmp <- count.categ; names(tmp) <- cnam;
					barplot(tmp, col=colo, main='Transcriptomic regions', ylab='number of peaks', ylim=c(0, max(tmp)*1.1))
					text(x=(1:length(tmp))*1.2-0.5, y=tmp, labels=paste(round(perc, 2), '%'), pos=3)
				}
				print(tmp)
			par(opar)
			##profile
			#hist(dNorm, breaks=((min(dNorm, na.rm=TRUE)/binsNorm-1):(max(dNorm, na.rm=TRUE)/binsNorm+1))*binsNorm, main='Normalised distribution', xlab='Distance')
			plot.dens(d=dNorm, d.e=dNorm.e, stasto=FALSE, main='Normalised profile')
			if(stasto) { plot.dens(d=dReNorm, d.e=dReNorm.e, stasto=TRUE, main='Normalised profile') }
			if(dex) {
				plot.dens(d=dexNorm, d.e=dexNorm.e, stasto=FALSE, main='Exonic-normalised profile')
				if(stasto) { plot.dens(d=dexReNorm, d.e=dexReNorm.e, stasto=TRUE, main='Exonic-normalised profile') }
			}
			##hist
			hTSS <- hist(dTSS.cut, breaks=((ifelse(draw.promo, (-1*borderTSS/bins-1), 0)):(borderTSS/bins+1))*bins, plot=FALSE)
			hTTS <- hist(dTTS.cut, breaks=((-1*abs(borderTTS)/bins-1):(ifelse(draw.TTSout, (borderTTS/bins+1), 0)))*bins, plot=FALSE)
			if(stasto) {
				hsta <- hist(dsta.cut, breaks=((-1*borderCDS/bins-1):(borderCDS/bins+1))*bins, plot=FALSE)
				hstop <- hist(dstop.cut, breaks=((-1*borderCDS/bins-1):(borderCDS/bins+1))*bins, plot=FALSE)
				ymax <- max(c(hTSS$counts, hTTS$counts, hsta$counts, hstop$counts), na.rm=TRUE)
			} else {
				ymax <- max(c(hTSS$counts, hTTS$counts), na.rm=TRUE)
			}
			plot(hTSS, main='Distribution surrounding the TSS', xlab='Distance', ylim=c(0, ymax))
			plot(hTTS, main='Distribution surrounding the TTS', xlab='Distance', ylim=c(0, ymax))
			if(stasto) {
				plot(hsta, main='Distribution surrounding the start codon', xlab='Distance', ylim=c(0, ymax))
				plot(hstop, main='Distribution surrounding the stop codon', xlab='Distance', ylim=c(0, ymax))
			}
		dev.off()
	}
}
