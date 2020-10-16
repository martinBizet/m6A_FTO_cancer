#!/usr/bin/env Rscript
argu <- commandArgs(TRUE)

###-------------
##	ARGUMENTS & HELP
#---------------

##Defaults
infile <- NULL; ipfile <- NULL; outfile <- NULL;
filtOn <- NA; direc <- 'Both';
thres <- log2(2); thresH <- log2(4); thresS <- 0.05;
print.help <- FALSE
draw.diff <- FALSE; shorten <- TRUE;
sel.categ <- c('near_TSS', '5UTR', 'near_Start', 'CDS', 'intron', 'near_Stop', '3UTR', 'near_TTS', 'Multiple', 'Intergenic'); maskTotal <- TRUE;

##Treat argu
argu <- strsplit(argu, split='=')
argu.1 <- argu[which(sapply(argu, length) == 1)]
if(length(argu.1) > 0) {
	for(el in argu.1) {
		if(el[1] == '-help') {
			print.help <- TRUE
		} else if (el[1] == '-draw') {
			draw.diff <- TRUE
		} else if (is.element(el[1], c('-ns', '-noNameShortening'))) {
			shorten <- FALSE
		} 
	}
}
if(print.help) {
	cat('Diff_RIP.R computes the dLIR and add localisation of the peaks in the IP vs INPUT plot. It can also generate a list of peaks moving between the two conditions and an differential analysis report.\n
		\t-i=|-inputFile=  name of the differential table file to use as input (output from differential mtool). [mandatory]\n
		\t-o=|-outputFile=  name of the Re-computed differential table file (output from this tool). [mandatory]\n
		\t-f=|-filterOn=  feature to use for filtering. Possible values are dLIR, LOR, LR or Discord. If provided, a .moving file containing only the regions passing the threshold is provided in addition to the outputFile\n
		\t-d=|-direction=  direction of the analysis. Can be Both, Hypo or Hyper. [default=Both]\n
		\t-thres=|-threshold=  threshold on the filterOn value. [default=2]\n
		\t-thresH=|-thresholdHigh= High stringency threshold (for visualisation purpose only) [default=4]\n
		\t-thres2=|-threshold2=  threshold on the nonfilterOn value (LR if filterOn=dLIR or LOR; dLIR if filterOn=LR or Discord). Use on for visualisation. [default=thres]\n
		\t-thresS=|-thresholdSignificativity=  threshold on the significativity value (used only if sLR_ and sIR_ tagged column are available in -i file. [default=0.05]\n
		\t-ip=|ipFile= file with IP for all sample (for heatmap when Diff_1vs1_RIPort is called by Diff_multi_RIPort)
		\t-draw |-draw= provide a differential analysis report into a .pdf file. The annotation categories to show can be specified as comma separated names [default=all]\n
		\t-ns |-noNameShortening avoid name shortening in the plots. Typical used with when names in inputFile are shorts\n
		\t-help  print this help\n
		Example\n
		-------\n
		Rscript Diff_RIP.R -i=</path/to/differential/table/>MCF12A_vs_SKBR3_at_peaks.txt -o=</path/to/output/differential/table/>dLIR_MCF12A_vs_SKBR3_at_peaks.txt -f=dLIR -thres=2 -thres2=1.5 -draw')
} else {
	argu <- argu[which(sapply(argu, length) > 1)]
	argu.keys <- sapply(argu, '[[', 1)
	argu.vals <- sapply(argu, '[[', 2)
	if(length(intersect(c('-i', '-inputFile'), argu.keys)) > 0) {
		infile <- argu.vals[which(is.element(argu.keys, c('-i', '-inputFile')))]
	} else {
		stop('Please provide an input file (-i or -inputFile argument)')
	}
	if(length(intersect(c('-o', '-outputFile'), argu.keys)) > 0) {
		outfile <- argu.vals[which(is.element(argu.keys, c('-o', '-outputFile')))]
	} else {
		stop('Please provide an output file (-o or -outputFile argument)')
	}
	if(length(intersect(c('-f', '-filterOn'), argu.keys)) > 0) { filtOn <- argu.vals[which(is.element(argu.keys, c('-f', '-filterOn')))] }
	if(length(intersect(c('-d', '-direction'), argu.keys)) > 0) { direc <- argu.vals[which(is.element(argu.keys, c('-d', '-direction')))] }
	if(length(intersect(c('-thres', '-threshold'), argu.keys)) > 0) { thres <- log2( as.numeric(argu.vals[which(is.element(argu.keys, c('-thres', '-threshold')))]) ) }
	if(length(intersect(c('-thresS', '-thresholdSignificativity'), argu.keys)) > 0) { thresS <- as.numeric(argu.vals[which(is.element(argu.keys, c('-thresS', '-thresholdSignificativity')))]) }
	if(length(intersect(c('-thresH', '-thresholdHigh'), argu.keys)) > 0) { thresH <- log2( as.numeric(argu.vals[which(is.element(argu.keys, c('-thresH', '-thresholdHigh')))]) ) }
	if(length(intersect(c('-thres2', '-threshold2'), argu.keys)) > 0) {
		thres2 <- log2( as.numeric(argu.vals[which(is.element(argu.keys, c('-thres2', '-threshold2')))]) )
	} else {
		thres2 <- thres
	}
	if(length(intersect(c('-ip', '-ipFile'), argu.keys)) > 0) {
		ipfile <- argu.vals[which(is.element(argu.keys, c('-ip', '-ipFile')))]
	}
	if(is.element('-draw', argu.keys)) {
		draw.diff <- argu.vals[which(argu.keys == '-draw')]
		sel.categ <- strsplit(draw.diff, split=',')[[1]]
		draw.diff <- TRUE
	}
	
	###-------------
	##	FUNCTIONS
	#---------------

	if(shorten) { require('Rlibstree') }
	transparency <- function(colo, alpha=0.3) { return(rgb(t(col2rgb(colo)/255), alpha=alpha)) }
	draw.ablines <- function(thres.dLIR, thres.IR, thres.LR) {
		abline(a=0, b=1, lty=2)
		abline(a=-1*thres.dLIR, b=1, lty=1)
		abline(a=thres.dLIR, b=1, lty=1)
		abline(v=-1:1*thres.IR, lty=c(1, 2, 1))
		abline(h=-1:1*thres.LR, lty=c(1, 2, 1))
	}
	
	###-------------------------------------
	##	Compute Fenrich_1, Fenrich_2, dLIR
	#---------------------------------------

	##Loading
	dec <- c('.', ',')[length(grep(',', strsplit(readLines(con=infile), split='\t')[[2]][1], fixed=TRUE))+1]
	diff.tab <- read.table(infile, sep='\t', header=TRUE, row.names=NULL, check.names=FALSE, as.is=TRUE, quote='', dec=dec, comment.char='')
	if(!is.null(ipfile)) { ip.tab <- read.table(ipfile, sep='\t', header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE, quote='', dec=dec, comment.char='') }
	##Detect sampNames (names of the compared columns)
	tmp <- colnames(diff.tab)[which(nchar(colnames(diff.tab)) >= 4)]
	tmp.2 <- substr(tmp, 1, 4)
	tagged <- FALSE
	withSig <- FALSE
	if(sum(is.element(c('sLR_', 'sIR_'), tmp.2)) == 2) {
		sigNames <- tmp[c(which(tmp.2 == 'sLR_'), which(tmp.2 == 'sIR_'))]
		withSig <- TRUE
	}
	if(sum(is.element(c('bL0_', 'bL1_', 'bI0_', 'bI1_'), tmp.2)) == 4) {
		inpNames <- tmp[c(which(tmp.2 == 'bI1_'), which(tmp.2 == 'bI0_'))]
		sampNames <- tmp[c(which(tmp.2 == 'bL1_'), which(tmp.2 == 'bL0_'))]
		tagged <- TRUE
	} else {
		warning('colnames do not contain bL0, bl1, bI0 and bI1 tags. The analysis will continue assuming that inputs contains the word Input and that condition1 is before condition0 in IP and Input')
		columnsFromInfile <- c('#Chr', 'Start', 'End', 'name', 'Length', 'Strand', 'LOR', 'LR', 'IR', 'abs_summit', 'annotation', 'gene', 'category')
		sampNames <- colnames(diff.tab)[which(!is.element(colnames(diff.tab), columnsFromInfile))]
		i.tmp <- grep('Input', sampNames, ignore.case=TRUE)
		if(length(i.tmp) == 2) {
			inpNames <- sampNames[i.tmp]
		} else {
			stop('Can not identify input samples. Please please provided bL0, bl1, bI0 and bI1 tags or in colnames or colnames containing the word input')
		}
		sampNames <- setdiff(sampNames, inpNames)
	}

	##Computes new features
	dLIR <- diff.tab$LR - diff.tab$IR
	Fen1 <- diff.tab[, sampNames[1]] - diff.tab[, inpNames[1]]
	Fen2 <- diff.tab[, sampNames[2]] - diff.tab[, inpNames[2]]
	i.LOR <- which(names(diff.tab) == 'LOR')
	if(withSig) { sigLR <- diff.tab[, sigNames[1]]; sigIR <- diff.tab[, sigNames[2]]; }	
	diff.tab <- cbind(diff.tab[, 1:(i.LOR-1)], Fenrich_1=Fen1, Fenrich_2=Fen2, dLIR=dLIR, diff.tab[, i.LOR:ncol(diff.tab)])
	
	###-------------------------
	##	Identify moving peaks
	#---------------------------
	
	##thres
	thres.IR <- log2(1.2)
	if(is.element(filtOn, c('dLIR', 'LOR'))) {
		thres.dLIR <- thres
		thres.LR <- thres2
	} else {
		thres.dLIR <- thres2
		thres.LR <- thres
	}
	
	#Add plot-regions
	for(ths in list(list(thres.dLIR, thres.LR, thres.IR, paste('Plot_Regions_LR', round(2^thres.LR, 2), 'dLIR', round(2^thres.dLIR, 2), sep='_')), 
					list(0, 0, 0, 'Plot_Regions_noThres'))) {
		thdLIR <- ths[[1]]; thLR <- ths[[2]]; thIR <- ths[[3]]; prNam <- ths[[4]]; 
		dLIR.sig <- apply(sapply( list(diff.tab$dLIR < -1*thdLIR, diff.tab$dLIR > thdLIR), as.numeric ) , 1, paste, collapse='')
		LR.sig <- apply(sapply( list(diff.tab$LR < -1*thLR, diff.tab$LR > thLR), as.numeric ) , 1, paste, collapse='')
		IR.sig <- apply(sapply( list(diff.tab$IR < -1*thIR, diff.tab$IR > thIR), as.numeric ) , 1, paste, collapse='')
		lvl <- c(
				'10_10_10','10_10_00','10_10_01',  '10_00_10','10_00_00','10_00_01',  '10_01_10','10_01_00','10_01_01',
				'00_10_10','00_10_00','00_10_01',  '00_00_10','00_00_00','00_00_01',  '00_01_10','00_01_00','00_01_01',
				'01_10_10','01_10_00','01_10_01',  '01_00_10','01_00_00','01_00_01',  '01_01_10','01_01_00','01_01_01'
				)
		names(lvl) <-  c(
						'IP<INPUT_IP-_INPUT-','IP<INPUT_IP-_INPUT0','IP<INPUT_IP-_INPUT+',  'IP<INPUT_IP0_INPUT-','IP<INPUT_IP0_INPUT0','IP<INPUT_IP0_INPUT+',  'IP<INPUT_IP+_INPUT-','IP<INPUT_IP+_INPUT0','IP<INPUT_IP+_INPUT+',
						'IP=INPUT_IP-_INPUT-','IP=INPUT_IP-_INPUT0','IP=INPUT_IP-_INPUT+',  'IP=INPUT_IP0_INPUT-','IP=INPUT_IP0_INPUT0','IP=INPUT_IP0_INPUT+',  'IP=INPUT_IP+_INPUT-','IP=INPUT_IP+_INPUT0','IP=INPUT_IP+_INPUT+',
						'IP>INPUT_IP-_INPUT-','IP>INPUT_IP-_INPUT0','IP>INPUT_IP-_INPUT+',  'IP>INPUT_IP0_INPUT-','IP>INPUT_IP0_INPUT0','IP>INPUT_IP0_INPUT+',  'IP>INPUT_IP+_INPUT-','IP>INPUT_IP+_INPUT0','IP>INPUT_IP+_INPUT+'
						)
		reg <- apply(cbind(dLIR.sig, LR.sig, IR.sig), 1, paste, collapse='_')
		reg <- factor(reg, levels=lvl)
		levels(reg) <- names(lvl)
		diff.tab[, prNam] <- as.character(reg)
	}
	diff.tab$Double_0 <- rep('OK', nrow(diff.tab))
	diff.tab$Double_0[which( apply(diff.tab[, inpNames], 1, sum) == 0 )]  <- 'No_Input'
	write.table(diff.tab, outfile, sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
	
	#filter
	if(!is.na(filtOn)) {
		if(filtOn == 'Discord') {
			if(withSig) {
				is.IR.stab <- ((sigIR >= thresS) & (diff.tab$IR > -1*thres.IR) & (diff.tab$IR < thres.IR)) #INPUT stable				
				#neg
				is.LR.mov <- ((sigLR < thresS) & (diff.tab$LR < -1*thres.LR)) #IP down
				is.IR.inv <- (diff.tab$IR >= thres.IR) #INPUT up
				i.neg <- which( is.LR.mov & (is.IR.inv | is.IR.stab) )
				#pos
				is.LR.mov <- ((sigLR < thresS) & (diff.tab$LR > thres.LR)) #IP up
				is.IR.inv <- (diff.tab$IR <= thres.IR) #INPUT down
				i.pos <- which( is.LR.mov & (is.IR.inv | is.IR.stab) )
			} else {		
				i.neg <- which( (diff.tab$LR < -1*thres) & (diff.tab$IR > -1*thres.IR) ) #IP down and INPUT up or stable
				i.pos <- which( (diff.tab$LR > thres) & (diff.tab$IR < thres.IR) ) #IP up and INPUT down or stable
			}
			filtOn <- 'LR'
		} else {
			if(withSig) {
				if(is.element(filtOn, c('dLIR', 'LOR'))) {
					i.neg <- which(( (sigLR < thresS) | (sigIR < thresS) ) & (diff.tab[, filtOn] < -1*thres))
					i.pos <- which(( (sigLR < thresS) | (sigIR < thresS) ) & (diff.tab[, filtOn] > thres))				
				} else {
					i.neg <- which((sigLR < thresS) & (diff.tab[, filtOn] < -1*thres))
					i.pos <- which((sigLR < thresS) & (diff.tab[, filtOn] > thres))
				}
			} else {
				i.neg <- which(diff.tab[, filtOn] < -1*thres)
				i.pos <- which(diff.tab[, filtOn] > thres)			
			}
		}
		diff.tab.pos <- diff.tab[i.pos,]
		diff.tab.neg <- diff.tab[i.neg,]
		if(direc == 'Both') {
			diff.sig <- diff.tab[c(i.pos, i.neg),]
			i.nsig <- setdiff(1:nrow(diff.tab), c(i.neg, i.pos))
		} else if (direc == 'Hypo') {
			diff.sig <- diff.tab[i.neg,]
			i.nsig <- setdiff(1:nrow(diff.tab), i.neg)
		} else {# if (direc == 'Hyper') {
			diff.sig <- diff.tab[i.pos,]
			i.nsig <- setdiff(1:nrow(diff.tab), i.pos)
		}
		diff.sig <- diff.sig[order(diff.sig[, filtOn], decreasing=FALSE),]
		diff.nsig <- diff.tab[i.nsig,]
		write.table(diff.sig, paste(outfile, 'moving', sep='.'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
	}
	
	###--------
	##	DRAW
	#----------
	
	if(draw.diff) {
		sampNa <- sampNames
		if(tagged) { sampNa <- substr(sampNa, 5, nchar(sampNa)) }
		if(shorten) {
			sampNa <- gsub(getLongestCommonSubstring(sampNa), '', sampNa) #can not be enough...
			tmp <- getLongestCommonSubstring(sampNa) #... let's run it twice ...
			if(!is.null(tmp) && sum(sampNa == tmp) == 0) { sampNa <- gsub(tmp, '', sampNa) } #... and change if it is possible
		}
		pdf(paste(outfile, 'pdf', sep='.'))
			lim <- range(c(diff.tab$LR, diff.tab$IR), na.rm=TRUE)
			lvl <- union(c('5UTR', 'CDS', 'intron', '3UTR', 'multiple'), names(table(diff.sig$category))) #at least that category more if it exists
			diff.sig$category <- factor(diff.sig$category, levels=lvl)
			diff.tab$category <- factor(diff.tab$category, levels=lvl)
			if(!is.na(filtOn) && nrow(diff.sig) > 0) {
				if(direc == 'Both') { #TRUE = Hyper; FALSE = Hypo;
					count <- cbind(table(diff.sig[, filtOn] > 0), table(diff.sig[, filtOn] > 0, diff.sig$category, useNA='i'))
					if(nrow(count) == 1) { #only one case...
						if(rownames(count) == "TRUE") { count <- rbind(rep(0,nrow(count)), count) } #...Add FALSE case before
						if(rownames(count) == "FALSE") { count <- rbind(count, rep(0,nrow(count))) } #...Add TRUE case after
					}
					rownames(count) <- c('Hypo', 'Hyper')
					denom <- nrow(diff.sig)
				} else if (direc == 'Hypo') { #TRUE = Reduced; FALSE = 'Not Reduced';
					count <- cbind(table(is.element(1:nrow(diff.tab), i.neg)), table(is.element(1:nrow(diff.tab), i.neg), diff.tab$category, useNA='i'))
					if(nrow(count) == 1) { #only one case
						if(rownames(count) == "TRUE") { count <- rbind(rep(0,nrow(count)), count) }
						if(rownames(count) == "FALSE") { count <- rbind(count, rep(0,nrow(count))) }
					}
					rownames(count) <- c('Not Reduced', 'Reduced')
					denom <- nrow(diff.tab)
				} else {# if (direc == 'Hyper') { #TRUE = Enhanced; FALSE = Not Enhanced;
					count <- cbind(table(is.element(1:nrow(diff.tab), i.pos)), table(is.element(1:nrow(diff.tab), i.pos), diff.tab$category, useNA='i'))
					if(nrow(count) == 1) { #only one case
						if(rownames(count) == "TRUE") { count <- rbind(rep(0,nrow(count)), count) }
						if(rownames(count) == "FALSE") { count <- rbind(count, rep(0,nrow(count))) }
					}
					rownames(count) <- c('Not Enhanced', 'Enhanced')
					denom <- nrow(diff.tab)
				}
				count <- rbind(count, 'expected'=c(nrow(diff.tab), table(diff.tab$category, useNA='i')))
				colnames(count)[1] <- 'Total'
				if(is.na(colnames(count)[ncol(count)])) {
					colnames(count)[ncol(count)] <- 'intergenic'
				} else if (is.element('Intergenic', colnames(count))) {
					colnames(count)[which(colnames(count) == 'Intergenic')] <- 'intergenic'
				}
				ord <- c('Total', '5UTR', 'CDS', 'intron', '3UTR', 'multiple', 'intergenic')
				if(is.element('near_TSS', colnames(count))) { ord <- c(ord[1], 'near_TSS', ord[2:length(ord)]) }
				if(is.element('near_Start', colnames(count))) { ord <- c(ord[1:2], 'near_Start', ord[3:length(ord)]) }
				if(is.element('near_Stop', colnames(count))) { i.tmp <- which(ord == '3UTR'); ord <- c(ord[1:(i.tmp-1)], 'near_Stop', ord[i.tmp:length(ord)]) }
				if(is.element('near_TTS', colnames(count))) { i.tmp <- which(ord == '3UTR'); ord <- c(ord[1:i.tmp], 'near_TTS', ord[(i.tmp+1):length(ord)]) }
				count <- count[, ord]
				count.2 <- count[1:2,]
			}
			##Compare Fold enrichment
			#boxplot
			pv1 <- t.test(x=diff.tab[, sampNames[1] ], y=diff.tab[, sampNames[2] ])$p.value
			pv2 <- t.test(x=diff.tab$Fenrich_2, y=diff.tab$Fenrich_1)$p.value
			boxplot(diff.tab[, sampNames[2:1]], names=sampNa[2:1], col=c('grey20', 'springgreen4'), notch=TRUE, outline=FALSE, ylab='Normalised reads (log2 FPKM)',
					main=paste('Peak heights (from', nrow(diff.tab), 'peaks) [t-test pv=', signif(pv1, 3), ']'))
			boxplot(diff.tab$Fenrich_2, diff.tab$Fenrich_1, names=sampNa[2:1], col=c('grey20', 'springgreen4'), notch=TRUE, outline=FALSE, ylab='Normalised reads (log2 fold enrichment)',
					main=paste('Number of reads (from', nrow(diff.tab), 'peaks) [t-test pv=', signif(pv2, 3), ']'))
			
			##Differential 
			#scatterplot
			if(!is.na(filtOn)) {
				colo <- sapply(c('red3', 'springgreen4', 'grey'), transparency, alpha=0.3)
				plot(diff.nsig$IR, diff.nsig$LR, pch=18, xlim=lim, ylim=lim, col=colo[3], cex=1.3, xlab='Input (log2 fold change)', ylab='IP (log2 fold change)',
					main=paste(paste(sampNa, collapse='.vs.'), '(from', nrow(diff.tab), 'peaks)'))
				if(direc != 'Hyper') { points(diff.tab.neg$IR, diff.tab.neg$LR, pch=18, col=colo[2], cex=1.3) }
				if(direc != 'Hypo') { points(diff.tab.pos$IR, diff.tab.pos$LR, pch=18, col=colo[1], cex=1.3) }
			} else {
				colo <- transparency('red3', alpha=0.3)
				plot(diff.tab$IR, diff.tab$LR, pch=18, xlim=lim, ylim=lim, col=colo, cex=1.3, xlab='Input (log2 fold change)', ylab='IP (log2 fold change)',
						main=paste(sampNa, collapse='.vs.'))
			}
			draw.ablines(thres.dLIR, thres.IR, thres.LR)
			#pie
			if(!is.na(filtOn) && nrow(diff.sig) > 0) {
				opar <- par(mfrow=c(1,1))
				if(!is.na(thresH)) { opar <- par(mfrow=c(1,2)) }
					if(direc == 'Both') {
						colo <- c('springgreen4', 'red3'); coloH <- c('lightslateblue', 'navyblue'); ylab <- 'Differential peaks (%)';
						leg1 <- 'Hypo'; leg2 <- 'Hyper';
					} else if(direc == 'Hypo') {
						colo <- c('springgreen4', 'grey20'); coloH <- c('seagreen2', 'darkgreen'); ylab <- 'Reduced peaks (%)';
						leg1 <- 'Reduced'; leg2 <- 'Not Reduced';
					} else { #if(direc == 'Hyper') {
						colo <- c('red3', 'grey20'); coloH <- c('salmon', 'firebrick'); ylab <- 'Enhanced peaks (%)';
						leg1 <- 'Enhanced'; leg2 <- 'Not enhanced';
					}
					pie(count.2[,1], col=colo, main=colnames(count.2)[1], radius=0.9, init.angle=90)
					legend('topleft', c(paste(leg1, ': ', count.2[1,1], ' peaks (', max(0, round((count.2[1,1]/denom)*100, 1), na.rm=TRUE), '%)', sep=''), 
								paste(leg2, ': ', count.2[2,1], ' peaks (', max(0, round((count.2[2,1]/denom)*100, 1), na.rm=TRUE), '%)', sep='')), cex=0.8, fill=colo)
					if(!is.na(thresH)) { 
						countH <- 100*table(abs(diff.sig[, filtOn]) > thresH)/nrow(diff.sig)
						plot(NA, NA, xlim=c(0, 0.8), ylim=c(0, 100), ylab=ylab, xaxt='n', xlab='')
						rect(xleft=0.1, ybottom=c(0,countH[1]), xright=0.3, ytop=c(countH[1], 100), col=coloH, lwd=2)
						legend('topright', paste(c('>', '<='), 2^thresH), cex=0.8, fill=coloH[2:1])
					}
				par(opar)
				
				##Categories
				colo.cat <- c('royalblue4', 'darkorange', 'springgreen4', 'cornflowerblue', 'darkmagenta', 'red3', 'yellow', 'pink', 'gray', 'grey60')
				names(colo.cat) <- c('near_TSS', '5UTR', 'near_Start', 'CDS', 'intron', 'near_Stop', '3UTR', 'near_TTS', 'Multiple', 'Intergenic')
				#scatterplot with categorical colours
				plot(diff.tab$IR, diff.tab$LR, pch=18, xlim=lim, ylim=lim, col=sapply(colo.cat, transparency, alpha=0.3)[as.character(diff.tab$category)], cex=1.3, xlab='Input (log2 fold change)', 
						ylab='IP (log2 fold change)', main=paste(paste(sampNa, collapse='.vs.'), '(from', nrow(diff.tab), 'peaks)'))
				draw.ablines(thres.dLIR, thres.IR, thres.LR)
				colo.kept <- colo.cat[levels(as.factor(diff.tab$category))]
				legend('topleft', names(colo.kept), col=colo.kept, pch=18)
				#Barplot
				ndiff <- sum(count.2[, 1])
				sel.categ <- intersect(sel.categ, colnames(count))
				if(direc == 'Both' && !maskTotal) {
					sel.categ <- c('Total', sel.categ)
				} else { #if Hyper or Hypo, we do not need the Total column
					count <- count[, 2:ncol(count)] 
				}
				count <- count[,sel.categ]
				count.2 <- count.2[,sel.categ]
				count.3 <- rbind(count.2, Total=apply(count.2, 2, sum))
				nshow <- sum(count.3[3,])
				count.p <- 100*count.3/nshow #differential
				count.e <- 100*count[3,]/sum(count[3,]) #expected (whole peak repartition)
				opar <- par(las=3, mar=c(6.1, 4.1, 4.1, 2.1))
					plot(NA, NA, xlim=c(0, ncol(count.p)+1), ylim=c(0, 1.1*max(rbind(count.p, count.e))), xaxt='n', main=paste('Repartition of the', ndiff, 'differential peaks [', ndiff-nshow, 'not shown]'), xlab='', ylab='Percentage of peaks')
					if(direc == 'Both') {
						if(!maskTotal) { rect(xleft=par('usr')[1], ybottom=par('usr')[3], xright=1.5, ytop=par('usr')[4], col='grey95', border='transparent') }
						rect(xleft=(1:ncol(count.p))-0.35, ybottom=rep(0, ncol(count.p)), xright=(1:ncol(count.p))-0.05, ytop=count.p[1,], col='springgreen4')
						rect(xleft=(1:ncol(count.p))-0.35, ybottom=count.p[1,], xright=(1:ncol(count.p))-0.05, ytop=count.p[3,], col='red3')
					} else {
						rect(xleft=(1:ncol(count.p))-0.35, ybottom=rep(0, ncol(count.p)), xright=(1:ncol(count.p))-0.05, ytop=count.p[1,], col=colo.cat[colnames(count.p)])
					}
					rect(xleft=(1:length(count.e))+0.05, ybottom=rep(0, length(count.e)), xright=(1:length(count.e))+0.35, ytop=count.e, col='grey') #expected
					axis(2, at=(0:5)*20, labels=(0:5)*20)
					axis(1, at=1:ncol(count.p), labels=colnames(count.p))
					if(direc == 'Both') { legend('topright', c('Hyper', 'Hypo'), fill=c('red3', 'springgreen4')) }	
					oo <- options("scipen"=-0.5, "digits"=2)
						sizeText <- 5; sizeSpace <- sizeText/2;
						text(x=(1:ncol(count.p))-(0.35/2), y=count.p[3,]+sizeSpace, labels=count.3[3,], pos=3, cex=0.7)
						text(x=(1:ncol(count.p))+(0.35/2), y=count.e+sizeSpace, labels=count[3,], pos=3, cex=0.7)
					options(oo)
				par(opar)
				#Pies by category
				opar <- par(mfrow=c(2,2), mar=c(0.6, 1.1, 3.1, 1.1))
					for(nam in setdiff(colnames(count.2), c('Total', 'multiple', 'intergenic'))) {
						tmp=sum(count.2[,nam])
						if(tmp > 0) {
							pie(count.2[,nam], col=colo, radius=0.67, init.angle=90, main=nam)						
							legend('topleft', c(paste(leg1, ': ', count.2[1, nam], ' peaks (', max(0, round((count.2[1, nam]/tmp)*100, 1), na.rm=TRUE), '%)', sep=''), 
								paste(leg2, ': ', count.2[2, nam], ' peaks (', max(0, round((count.2[2, nam]/tmp)*100, 1), na.rm=TRUE), '%)', sep='')), cex=0.8, fill=colo)
						} else {
							plot(NA, NA, xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1), col='white', main=nam, bty='n')
							text(x=0.5, y=0.5, labels='No candidate available', cex=2)
						}
					}
				par(opar)
			
				##Heatmap
				require('gplots')
				tmp <- diff.sig[, sampNames]
				tmp <- as.matrix(tmp[, sampNames])
				colnames(tmp) <- sampNa
				if(nrow(tmp) > 1) {
					#color.pal <- colorRampPalette(c('white', 'rosybrown1', 'red1'), bias=0.9) #white to red palette (a lot of pink)
					color.pal <- colorRampPalette(c('darkblue', 'white', 'red1'), bias=0.9) #blue to red palette
					name <- paste(sampNa, collapse=' \nvs ')
					cexCol <- max(0.01, round(1.5/((ncol(tmp)/40)+1), 2))
					distfun <- function(x, m='euclidean') { return(dist(x, method=m)) }
					hclustfun <- function(d, m='ward') { return(hclust(d, method=m)) }
					tmp <- tmp[order(tmp[,1], decreasing=TRUE),]
					heatmap.2(tmp, distfun=distfun, hclustfun=hclustfun, cexRow=max(0.01, round(1/((nrow(tmp)/40)+1), 2)), cexCol=cexCol, scale="none", Rowv=TRUE, Colv=NA, dendrogram='row', col=color.pal,
								margins=c((1.0 + max(nchar(colnames(tmp)))/2)*cexCol, 8), labRow=diff.sig$gene, main=paste(ndiff, 'differential peaks \n(', name, ')'), key=TRUE, keysize=1.0,
								density.info='none', trace='none')
					if(!is.null(ipfile)) {
						tmp <- as.matrix(ip.tab[diff.sig$name,])
						color <- rep('grey', ncol(tmp)); color[which(substr(colnames(tmp), 1, 4) == 'IP1_')] <- 'purple';
						heatmap.2(tmp, distfun=distfun, hclustfun=hclustfun, ColSideColors=color, cexRow=max(0.01, round(1/((nrow(tmp)/40)+1), 2)), cexCol=cexCol, scale="none", Rowv=TRUE, Colv=TRUE,
								dendrogram='row', col=color.pal, margins=c((1.0 + max(nchar(colnames(tmp)))/2)*cexCol, 8), labRow=diff.sig$gene, main=paste(ndiff, 'differential peaks \n(', name, ')'),
								key=TRUE, keysize=1.0, density.info='none', trace='none')
					}
				} else {
					plot(NA, NA, xaxt='n', xlab='', yaxt='n', ylab='', xlim=c(0, 1), ylim=c(0, 1), col='white', main=nam, bty='o')
					text(x=0.5, y=0.5, labels='The heatmap needs at least 2 candidates', cex=2)
				}
			}
		dev.off()
	}
}
