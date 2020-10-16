path.rsem <- '../../Data/GDC'
path.gmt <- '../../Data/EMTsignature'
path.leg <- getwd()

###Functions###
read.tabsep <- function (file, path = getwd(), sep = "\t", header = TRUE, row.names = 1, 
    check.names = FALSE, as.is = TRUE, quote = "", comment.char = "", ...) {
    fi <- file.path(path, file)
    read.table(fi, sep = sep, header = header, row.names = row.names, 
        check.names = check.names, as.is = as.is, quote = quote, 
        comment.char = comment.char, ...)
}

write.tabsep <- function (x, file, path = getwd(), first = "", sep = "\t", col.names = TRUE, 
    row.names = TRUE, quote = FALSE, ...) {
    if (row.names && col.names) {
        colnames(x)[1] <- paste(first, colnames(x)[1], sep = sep)
    }
    write.table(x, file.path(path, file), sep = sep, col.names = col.names, 
        row.names = row.names, quote = quote, ...)
}

###Normals vs Cancers###

#Load and Prepare
rsem <- read.tabsep('FTO.txt', path.rsem, header=FALSE, sep='', row.names=NULL)
rsem$V3 <- sapply(strsplit(rsem$V3, split=','), '[[', 1)
rsem$V4 <- sapply(strsplit(rsem$V4, split=','), '[[', 1)
rsem$sampType <- substr(rsem$V4, nchar(rsem$V4)-nchar('11A'), nchar(rsem$V4)-1)
rsem <- rsem[which(is.element(rsem$sampType, c('-01', '-11'))),]
tmp <- table(rsem$V3, rsem$sampType)
f <- function(v) { return(sum(v == 0) == 0) }
tmp <- rownames(tmp)[which(apply(tmp, 1, f))]
rsem <- rsem[which(is.element(rsem$V3, tmp)),]
rsem$log2RSEM <- log2(rsem$V2)

#fuse LUNG
rsem$V3 <- as.factor(rsem$V3)
levels(rsem$V3) <- gsub('LUSC', 'LUNG', gsub('LUAD', 'LUNG', levels(rsem$V3)))
rsem$V3 <- as.character(rsem$V3)
rsem$categ <- factor(paste(rsem$V3, rsem$sampType, sep=''), levels=paste.comb(levels(as.factor(rsem$V3)), c('11', '01'), sep='-'))
levels(rsem$categ) <- gsub('LUSC', 'LUNG', gsub('LUAD', 'LUNG', levels(rsem$categ)))

#Stat
lvl <- levels(as.factor(rsem$categ))
nb.cat <- length(lvl)
t.dof <- t.efSize <- t.ci <- t.stat <- pv <- nb.canc <- nb.norm <- med.canc <- med.norm <- rep(NA, nb.cat/2)


cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d
}

for(i in 1:length(pv)) {
	canc <- rsem$log2RSEM[which(rsem$categ == lvl[2*i])]
	norm <- rsem$log2RSEM[which(rsem$categ == lvl[(2*i)-1])]
	nb.canc[i] <- length(canc); nb.norm[i] <- length(norm);
	med.canc[i] <- median(canc); med.norm[i] <- median(norm);
	if(length(canc) > 1 && length(norm) > 1) {
		t.obj <- t.test(canc, norm)
		pv[i] <- t.obj$p.value
		t.stat[i] <- t.obj$statistic
		t.ci[i] <- paste(t.obj$conf.int, collapse=':')
		t.efSize[i] <- cohens_d(canc, norm)
		t.dof[i] <- t.obj$parameter
	} else {
		pv[i] <- 1
	}
}
thres <- c(0.05, 0.01, 0.001)
f <- function(el, thres) { return(paste(rep('*', sum(thres > el)), collapse='')) }
star <- sapply(pv, f, thres)
res <- cbind(nb.canc, nb.norm, med.canc, med.norm, pv, star, t.stat, t.ci, t.efSize, t.dof)
rownames(res) <- levels(as.factor(rsem$V3));
colnames(res) <- c('Number_Cancers', 'Number_Normals', 'Median_Cancers', 'Median_Normals', 't.test_pvalue', '0.05|0.01|0.001', 't.test_stat', 't.test_CI', 't.test_efSize', 't.test_df');
write.tabsep(res, 'GDC_ttest.txt', path.rsem)

palettes <- list(c('#29C710', 'darkgreen', '#EB1323', 'red4')) #try all the colours you want

#Draw
for(i in 1:length(palettes)) {
	print(i)
	colnorm <- palettes[[i]][1]; colnorm.border <- palettes[[i]][2]; colcanc <- palettes[[i]][3]; colcanc.border <- palettes[[i]][4];
	pdf(file.path(path.rsem, paste('Boxplot_GDC_ttest_palette-', i,'.pdf', sep='')), height=10, width=17)
		opar <- par(mar=c(12.1, 5.1, 4.1, 2.1), las=3)
			idx <- which( pv < thres[1] & med.norm > med.canc)
			cat.down <- rownames(res)[idx]
			rsem.down <- rsem[which(is.element(rsem$V3, cat.down)),]
			rsem.down$categ <- factor(paste(rsem.down$V3, rsem.down$sampType, sep=''), levels=paste.comb(levels(as.factor(rsem.down$V3)), c('11', '01'), sep='-'))
			nb.cat.down <- length(levels(rsem.down$categ))
			plot(NA, NA, xlim=c(0, nb.cat.down), ylim=range(rsem.down$log2RSEM), xlab='', ylab='Expression (log2 FPKM)', main='FTO expression', xaxt='n', yaxt='n', cex.lab=2, cex.main=2)
			axis(1, at=2*(1:(nb.cat.down/2))-0.5, labels=levels(as.factor(rsem.down$V3)), cex.axis=1.9)		
			rect(xleft=c(par('usr')[1], 2*(1:(nb.cat.down/2 -1))+0.5), ybottom=par('usr')[3], xright=c(2*(2:(nb.cat.down/2))+0.5, par('usr')[2]), ytop=par('usr')[4], col=c('grey80', 'white'), border=FALSE)
			b <- boxplot(rsem.down$log2RSEM ~ rsem.down$categ, add=TRUE, lwd=2, boxwex=0.6, col=c(colnorm, colcanc), staplecol=c(colnorm.border, colcanc.border), xaxt='n',
							whiskcol=c(colnorm.border, colcanc.border), medcol=c(colnorm.border, colcanc.border), boxcol=c(colnorm.border, colcanc.border), outpch=NA, cex.axis=2) #Fuks style
			tmp <- matrix(b$stat[5,], nr=2)
			text(x=2*(1:(nb.cat.down/2))-0.5, y=apply(tmp, 2, max)+0.7, label=star[idx], col='darkred', cex=2)
			legend('topright', c('Normal', 'Cancer'), fill=c(colnorm, colcanc), border=c(colnorm.border, colcanc.border, lwd=2), bg='white', cex=2.5)		
		par(opar)
	dev.off()
}

##EMT##

#Load and Prepare
rsem <- read.tabsep('ROKAVECandFTO_GDC.txt', path.rsem, header=FALSE, sep='', row.names=NULL)
rsem$V3 <- sapply(strsplit(rsem$V3, split=','), '[[', 1)
rsem$V4 <- sapply(strsplit(rsem$V4, split=','), '[[', 1)
rsem$sampType <- substr(rsem$V4, nchar(rsem$V4)-nchar('11A'), nchar(rsem$V4)-1)
rsem <- rsem[which(is.element(rsem$sampType, c('-01', '-11'))),]
tmp <- table(rsem$V3, rsem$sampType)
f <- function(v) { return(sum(v == 0) == 0) }
tmp <- rownames(tmp)[which(apply(tmp, 1, f))]
rsem <- rsem[which(is.element(rsem$V3, tmp)),]
rsem$log2RSEM <- log2(rsem$V2)

#fuse LUNG
rsem$V3 <- as.factor(rsem$V3)
levels(rsem$V3) <- gsub('LUSC', 'LUNG', gsub('LUAD', 'LUNG', levels(rsem$V3)))
rsem$V3 <- as.character(rsem$V3)
rsem$categ <- factor(paste(rsem$V3, rsem$sampType, sep=''), levels=paste.comb(levels(as.factor(rsem$V3)), c('11', '01'), sep='-'))
levels(rsem$categ) <- gsub('LUSC', 'LUNG', gsub('LUAD', 'LUNG', levels(rsem$categ)))

#Define EMTlow and EMT high
thres.pathoi <- c(0.2, 0.8)
pathw <- 'ROKAVEC_2017_SCIREP'
chip <- read.tabsep('ENSEMBL_human_gene_With-MiNuS-from-EMTsignas.chip', path.gmt)
chip <- chip[which(rownames(chip) != 'ENSG00000268439'),] #Deprecated ID
gmt <- file.path(path.gmt, 'EMTbiblio_Signed.gmt')
gmt.l <- readLines(gmt)
gmt.l <- strsplit(gmt.l, split='\t')
names(gmt.l) <- sapply(gmt.l, '[[', 1)
signa <- gmt.l[[pathw]][3:length(gmt.l[[pathw]])]
emt.guys <- gsub('MiNuS_', '', signa) #use the full signature in correlation plots
signa <- rownames(chip)[which(is.element(chip[, 'Gene Symbol'], signa))]
rsem$categ <- as.character(rsem$categ)
rsem.prim <- rsem[which(as.character(rsem$sampType) == '-01'),] #restrict to primary
correct <- rep(1, length(signa)); names(correct) <- signa;
i.tmp <- which(substr(signa, 1, nchar('MiNuS')) == 'MiNuS')
correct[i.tmp] <- -1; names(correct)[i.tmp] <- substr(names(correct), nchar('MiNuS_')+1, nchar(names(correct)))[i.tmp];
rsem.FTO <- rsem.prim[which(rsem.prim$V1 == 'ENSG00000140718.17'),]
rsem.meta <- rsem.prim[which(rsem.prim$V1 != 'ENSG00000140718.17'),]
rsem.meta$log2RSEM <- correct[sapply(strsplit(rsem.meta$V1, split='\\.'), '[[', 1)]*rsem.meta$log2RSEM
get.categ <- function(idx, mat, thres.pathoi) {
	mat <- mat[idx,]
	mat$log2RSEMscaled <- mat$log2RSEM
	i.NNA <-  which(!is.na(mat$log2RSEM))
	for(g in unique(mat$V1)) { #Bias to high expressed genes according to reviewer... So let scale the genes first
		print(g)
		i.tmp <- intersect(i.NNA, which(mat$V1 == g))
		tmp <- mat$log2RSEM[i.tmp]
		tmp[which(is.infinite(tmp) & tmp > 0)] <- max(tmp[which(!is.infinite(tmp))]) # replace Inf by maximal value
		tmp[which(is.infinite(tmp) & tmp < 0)] <- min(tmp[which(!is.infinite(tmp))]) # replace -Inf by minimal value
		mat$log2RSEMscaled[i.tmp] <- scale(tmp)
	}
	#meta <- tapply(mat$log2RSEM, mat$V4, mean, na.rm=TRUE) #not scaled
	meta <- tapply(mat$log2RSEMscaled, mat$V4, mean, na.rm=TRUE) #scaled
	meta.thres <- quantile(meta, thres.pathoi)
	categMETA <- rep(as.character(unique(mat$categ)), length(meta)); names(categMETA) <- names(meta);
	i.low <- which(meta < meta.thres[1]); i.high <- which(meta > meta.thres[2]);
	categMETA[i.low] <- paste(categMETA[i.low], 'Low', sep='-'); categMETA[i.high] <- paste(categMETA[i.high], 'High', sep='-');
	return(categMETA)
}
categMETA <- tapply(1:nrow(rsem.meta), rsem.meta$V3, get.categ, mat=rsem.meta, thres.pathoi=thres.pathoi, simplify=FALSE)
# i.test <- which(is.element(rsem.meta$V3, c('TCGA-BLCA', 'TCGA-BRCA'))) #TEST
# categMETA <- tapply(i.test, rsem.meta$V3[i.test], get.categ, mat=rsem.meta, thres.pathoi=thres.pathoi, simplify=FALSE) #TEST
rsem.FTO$categMETA <- rep(NA, nrow(rsem.FTO))
for(i in 1:length(categMETA)) {
	idx <- which(is.element(rsem.FTO$V4, names(categMETA[[i]])))
	rsem.FTO$categMETA[idx] <- categMETA[[i]][rsem.FTO[idx,]$V4]
}
rsem.FTO <- rsem.FTO[grep('Low|High', rsem.FTO$categMETA),]
tmp <- tapply(rsem.FTO$log2RSEM, rsem.FTO$V4, mean); rsem.FTO$log2RSEM <- tmp[rsem.FTO$V4];
tmp <- tapply(rsem.FTO$V2, rsem.FTO$V4, mean); rsem.FTO$V2 <- tmp[rsem.FTO$V4];
rsem.FTO <- rsem.FTO[which(!duplicated(rsem.FTO$log2RSEM)),]
rsem.FTO$categMETA <- factor(rsem.FTO$categMETA, levels=paste.comb(levels(as.factor(rsem.FTO$categ)), c('Low', 'High'), sep='-'))
	
#Stat
lvl <- levels(rsem.FTO$categMETA)
nb.cat <- length(lvl)
t.dof <- t.efSize <- t.ci <- t.stat <- pv <- nb.EMThigh <- nb.EMTlow <- med.EMThigh <- med.EMTlow <- rep(NA, nb.cat/2)
for(i in 1:length(pv)) {
	EMThigh <- rsem.FTO$log2RSEM[which(rsem.FTO$categMETA == lvl[2*i])]
	EMTlow <- rsem.FTO$log2RSEM[which(rsem.FTO$categMETA == lvl[(2*i)-1])]
	nb.EMThigh[i] <- length(EMThigh); nb.EMTlow[i] <- length(EMTlow);
	med.EMThigh[i] <- median(EMThigh); med.EMTlow[i] <- median(EMTlow);
	if(length(EMThigh) > 1 && length(EMTlow) > 1) {
		t.obj <- t.test(EMThigh, EMTlow)
		pv[i] <- t.obj$p.value
		t.stat[i] <- t.obj$statistic
		t.ci[i] <- paste(t.obj$conf.int, collapse=':')
		t.efSize[i] <- cohens_d(EMThigh, EMTlow)
		t.dof[i] <- t.obj$parameter		
	} else {
		pv[i] <- 1
	}
	names(pv) <- levels(as.factor(rsem.FTO$V3))
}
thres <- c(0.05, 0.01, 0.001)
f <- function(el, thres) { return(paste(rep('*', sum(thres > el)), collapse='')) }
star <- sapply(pv, f, thres)
res <- cbind(nb.EMThigh, nb.EMTlow, med.EMThigh, med.EMTlow, pv, star, t.stat, t.ci, t.efSize, t.dof)
rownames(res) <- levels(as.factor(rsem.meta$V3)); colnames(res) <- c('Number_EMThigh', 'Number_EMTlow', 'Median_EMThigh', 'Median_EMTlow', 't.test_pvalue', '0.05|0.01|0.001', 't.test_stat', 't.test_CI', 't.test_efSize', 't.test_df');
write.tabsep(res, 'GDC_EMTscaled_ttest.txt', path.rsem)

#Draw
palettes <- list(c('#fc9272', 'black', '#de2d26', 'black')) #try all the colours you want

source(file.path(path.leg, 'legendx.R'))
for(i in 1:length(palettes)) {
	print(i)
	colEMTlow <- palettes[[i]][1]; colEMTlow.border <- palettes[[i]][2]; colEMThigh <- palettes[[i]][3]; colEMThigh.border <- palettes[[i]][4];
	pdf(file.path(path.rsem, paste('Boxplot_GDC-EMTscaled_ttest_palette-', i,'.pdf', sep='')), height=10, width=17)
		opar <- par(mar=c(12.1, 5.1, 4.1, 2.1), las=3)
			idx <- which( pv < thres[1])
			cat.down <- rownames(res)[idx]
			rsem.down <- rsem.FTO[which(is.element(rsem.FTO$V3, cat.down)),]
			rsem.down$categ <- factor(paste(rsem.down$V3, rsem.down$sampType, sep=''), levels=paste(levels(as.factor(rsem.down$V3)), '01', sep='-'))
			rsem.down$categMETA <- factor(as.character(rsem.down$categMETA), levels=paste.comb(levels(as.factor(rsem.down$categ)), c('Low', 'High'), sep='-'))
			nb.cat.down <- length(levels(rsem.down$categMETA))
			plot(NA, NA, xlim=c(0, nb.cat.down)+0.5, ylim=range(rsem.down$log2RSEM), xlab='', ylab='Expression (log2 FPKM)', main='FTO expression', xaxt='n', yaxt='n', cex.lab=2, cex.main=2)
			axis(1, at=2*(1:(nb.cat.down/2))-0.5, labels=levels(as.factor(rsem.down$V3)), cex.axis=1.9)
			rect(xleft=c(par('usr')[1], 2*(1:(nb.cat.down/2 -1))+0.5), ybottom=par('usr')[3], xright=c(2*(2:(nb.cat.down/2))+0.5, par('usr')[2]), ytop=par('usr')[4], col=c('grey80', 'white'), border=FALSE)
			b <- boxplot(rsem.down$log2RSEM ~ rsem.down$categMETA, add=TRUE, lwd=2, boxwex=0.6, col=c(colEMTlow, colEMThigh), staplecol=c(colEMTlow.border, colEMThigh.border), xaxt='n',
							whiskcol=c(colEMTlow.border, colEMThigh.border), medcol=c(colEMTlow.border, colEMThigh.border), boxcol=c(colEMTlow.border, colEMThigh.border), outpch=NA, cex.axis=2)
			tmp <- matrix(b$stat[5,], nr=2)
			text(x=2*(1:(nb.cat.down/2))-0.5, y=apply(tmp, 2, max)+0.7, label=star[idx], col='darkred', cex=2)
			legend('topright', c(expression(paste('EMT ', bold('low'))), expression(paste('EMT ', bold('high')))), fill=c(colEMTlow, colEMThigh), box.cex=c(1.1,0.8), border=c(colEMTlow.border, colEMThigh.border, lwd=2), bg='white', cex=2.5)	
		par(opar)
	dev.off()
}

