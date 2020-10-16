path.res <- path.data <- '../../Data/GDS1439'

###Functions###
read.tabsep <- function (file, path = getwd(), sep = "\t", header = TRUE, row.names = 1, 
    check.names = FALSE, as.is = TRUE, quote = "", comment.char = "", ...) {
    fi <- file.path(path, file)
    read.table(fi, sep = sep, header = header, row.names = row.names, 
        check.names = check.names, as.is = as.is, quote = quote, 
        comment.char = comment.char, ...)
}
cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d
}
boxplot.clin <- function(v, categ, nam, lvl=NULL, col.box=NULL, col.dots=NULL, dot.cex=1.5, ylim=NULL) {
	if(is.null(lvl)) {
		categ <- as.factor(categ)
		lvl <- levels(categ)
	} else {
		categ <- factor(categ, levels=lvl)
	}
	min.shapi <- 1; thres.norm <- 0.05; thres.homo <- 0.05;
	for(z in 1:length(lvl)) { #test normality
		print(lvl[z])
		tmp <- shapiro.test(v[which(categ == lvl[z])])
		print(tmp)
		min.shapi <- min(min.shapi, tmp$p.value)
	}
	lvl <- paste(lvl, '(', table(categ), ')')
	levels(categ) <- lvl
#~ 	multi.lvl.test <- oneway.test; multi.lvl.name <- 'oneway.test';
	multi.lvl.test <- kruskal.test; multi.lvl.name <- 'kruskal.test';
	two.lvl.test <- t.test; two.lvl.name <- 't.test';
	if(min.shapi < thres.norm) { # at least one is significant (i.e. not normal)
		multi.lvl.test <- kruskal.test; multi.lvl.name <- 'kruskal.test';
		two.lvl.test <- wilcox.test; two.lvl.name <- 'wilcox.test';
	} else { #if normal check heteroscedaticity
		require('lmtest')	
		tmp <- bptest(lm(expr$logExpr ~ expr$case))
		print(tmp)
		if(tmp$p.value >= thres.homo) {
			multi.lvl.test <- function(f) { return(oneway.test(f, var.equal=TRUE)) }; multi.lvl.name <- 'oneway.test with equal variance';
			two.lvl.test <- function(f) { return(t.test(f, var.equal=TRUE)) }; two.lvl.name <- 't.test with equal variance';			
		}
	}
	
	if(length(lvl) > 2) {
		print(multi.lvl.name)
		tmp <- multi.lvl.test(v ~ categ)
		pv <- tmp$p.value
		print(tmp)
		for(z in 1:length(lvl)) { for(w in (z+1):length(lvl)) { print(lvl[c(w,z)]); print(cohens_d(v[which(categ == levels(categ)[z])], v[which(categ == levels(categ)[w])])); }}
	} else {
		print(two.lvl.name)
		tmp <- two.lvl.test(v ~ categ)
		pv <- tmp$p.value
		print(tmp)
		print(cohens_d(v[which(categ == levels(categ)[1])], v[which(categ == levels(categ)[2])]))
	}
	if(is.null(ylim)) { ylim <- range(v) }
	if(!is.null(col.box)) {
		boxplot(v ~ categ, ylim=ylim, main=paste(nam[1], '(pv=', signif(pv, 3), ')'), ylab=nam[2], col=col.box, outline=FALSE) #is.null(col.dots)
	} else {
		boxplot(v ~ categ, ylim=ylim, main=paste(nam[1], '(pv=', signif(pv, 3), ')'), ylab=nam[2], outline=FALSE) #is.null(col.dots)
	}
	if(!is.null(col.dots)) {
		stripchart(v ~ categ, method='jitter', vertical=TRUE, add=TRUE, jitter=0.2, pch=19, col=col.dots, cex=dot.cex)
	}
}

## Loading
expr <- read.tabsep('GDS1439-209702_at.txt', path.data, header=FALSE)
colnames(expr) <- c('SampDescr', 'Expr', 'Rnk')
expr$logExpr <- log2(expr$Expr)
expr$case <- as.factor(sapply(strsplit(expr$SampDescr, split=' '), '[[', 1)); levels(expr$case) <- c('Benign', 'Localized', 'Metastatic');
pdf(file.path(path.res, 'Boxplot_FTO_prostateGDS1439.pdf'))
	boxplot.clin(v=expr$logExpr, categ=expr$case, nam=c('', 'FTO expression (log2)'), col.box=c('mistyrose', 'salmon', 'red3'), col.dots=rep('black', 3), dot.cex=1, ylim=c(10.0, 11.5))
dev.off()
