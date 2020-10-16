#!/usr/bin/env Rscript
argu <- commandArgs(TRUE)

###-------------
##	ARGUMENTS & HELP
#---------------

##Defaults
print.help=FALSE
outfile=NA
col.split='annotation'
col.del=c('gene')
new.names=c('RefSeq_ID', 'GeneSymbol', 'tr_Position', 'rel_dist', 'dTSS', 'dTTS', 'dStart', 'dStop', 'dExTSS', 'dExTTS', 'dExStart', 'dExStop')
sep.new.row=';'
sep.new.col=','
sep.oth.col='\t'

##Treat argu
argu <- strsplit(argu, split='=')
argu.1 <- argu[which(sapply(argu, length) == 1)]
if(length(argu.1) > 0) {
	for(el in argu.1) {
		if(el[1] == '-help') {
			print.help <- TRUE
		} else if (el[1] == '-draw') {
			draw.diff <- TRUE
		}
	}
}
if(print.help) {
	cat('Script to duplicate lines according to one column.\n
		\t-i=|-inputFile=  name of the file to be duplicated. [mandatory]\n
		\t-o=|-outputFile=  name of the duplicated file. [default= < -inputFile >_duplic.txt]\n
		\t-col=|-columnToSplit=  name or index of the column to be splitted in multiple rows and columns. [default=annotation]\n
		\t-n=|-names=  coma separated names of the new columns [default=RefSeq_ID,GeneSymbol,tr_Position,rel_dist,dTSS,dTTS,dStart,dStop,dExTSS,dExTTS,dExStart,dExStop]\n
		\t-sRow=|-separatorNewRows=  separator to use for the generation of new rows. [default=;]\n
		\t-sCol=|-separatorNewColumns= separator to use for the generation of new columns. [default=,]\n
		\t-s=|-generalSeparator= general separator to define columns in the whole file. [default=\\t]\n
		\t-d=|-deleteColumns= coma separated names or indices (mix of names and indices are not allowed) of columns to suppress. [default=gene]\n
		\t-help  print this help\n
		Example\n
		-------\n
		Rscript duplic_lines.R -i=</path/to/differential/table/>.txt.moving')
} else {
	argu <- argu[which(sapply(argu, length) > 1)]
	argu.keys <- sapply(argu, '[[', 1)
	argu.vals <- sapply(argu, '[[', 2)
	if(length(intersect(c('-i', '-inputFile'), argu.keys)) > 0) {
		fnam <- argu.vals[which(is.element(argu.keys, c('-i', '-inputFile')))]
	} else {
		stop('Please provide an input file (-i or -inputFile argument)')
	}
	if(length(intersect(c('-o', '-outputFile'), argu.keys)) > 0) { outfile <- argu.vals[which(is.element(argu.keys, c('-o', '-outputFile')))] }
	if(length(intersect(c('-col', '-columnToSplit'), argu.keys)) > 0) {
		col.split <- argu.vals[which(is.element(argu.keys, c('-col', '-columnToSplit')))]
		if(length(grep("^[[:digit:]]*$", col.split)) == 1) { col.split <- as.numeric(col.split) }
	}
	if(length(intersect(c('-d', '-deleteColumns'), argu.keys)) > 0) {
		col.del <- strsplit(argu.vals[which(is.element(argu.keys, c('-d', '-deleteColumns')))], split=',')[[1]]
		if(length(grep("^[[:digit:]]*$", col.del)) == length(col.del)) { col.del <- as.numeric(col.del) }
	}
	if(length(intersect(c('-n', '-names'), argu.keys)) > 0) { new.names <- strsplit(argu.vals[which(is.element(argu.keys, c('-n', '-names')))], split=',')[[1]] }
	if(length(intersect(c('-sRow', '-separatorNewRows'), argu.keys)) > 0) { sep.new.row <- argu.vals[which(is.element(argu.keys, c('-sRow', '-separatorNewRows')))] }
	if(length(intersect(c('-sCol', '-separatorNewColumns'), argu.keys)) > 0) { sep.new.col <- argu.vals[which(is.element(argu.keys, c('-sCol', '-separatorNewColumns')))] }
	if(length(intersect(c('-s', '-generalSeparator'), argu.keys)) > 0) { sep.oth.col <- argu.vals[which(is.element(argu.keys, c('-ip', '-generalSeparator')))] }
	
	###-------------
	##	FUNCTIONS
	#---------------

	if(is.na(outfile)) { outfile <- paste(fnam, 'duplic.txt', sep='_') }
	a <- read.table(fnam, row.names=NULL, comment.char='', header=TRUE, as.is=TRUE, check.names=FALSE, sep=sep.oth.col)
	rownames(a) <- paste('p', 1:nrow(a), sep='')
	if(is.numeric(col.split)) { col.split <- colnames(a)[col.split] }
	if(is.numeric(col.del)) {
		a <- a[, -col.del]
	} else {
		a <- a[, which(!is.element(colnames(a), col.del))]
	}
	annot <- strsplit(a[,col.split], split=sep.new.row)
	names(annot) <- paste(rownames(a), '.nGum5eVr8o.', sep='')
	annot <- unlist(annot)
	annot <- sapply(annot, strsplit, split=sep.new.col)
	nb.col <- as.numeric(levels(as.factor(sapply(annot, length))))
	nb.col <- nb.col[order(as.numeric(nb.col), decreasing=TRUE)[1]]	
	f <- function(v, nb.col) {
		if(length(v) == 1 && is.na(v)) { return(rep(NA, nb.col)) }
		return(v)
	}
	annot <- t(as.data.frame(sapply(annot, f, nb.col=nb.col), stringsAsFactors=FALSE))
	colnames(annot) <- c('RefSeq_ID', 'GeneSymbol', 'tr_Position', 'rel_dist', 'dTSS', 'dTTS', 'dStart', 'dStop', 'dExTSS', 'dExTTS', 'dExStart', 'dExStop')
	a.dup <- a[sapply(strsplit(rownames(annot), split='.nGum5eVr8o.'), '[[', 1),]
	i <- which(colnames(a.dup) == col.split)
	a.dup <- cbind(a.dup[, 1:(i-1), drop=FALSE], annot, a.dup[,(i+1):ncol(a.dup), drop=FALSE])
	write.table(a.dup, outfile, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
}
