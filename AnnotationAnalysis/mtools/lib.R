# ------------------------------------------------------------
#
# Matthieu Defrance - ULB 2016
# lib.R
# 
# ------------------------------------------------------------

# ------------------------------------------------------------
# external libraries
# ------------------------------------------------------------
suppressMessages(library(Rsubread))
fc <- featureCounts

# ------------------------------------------------------------
# CSV
# ------------------------------------------------------------
read.CSV <- function(filename)
{
    data <- read.csv(filename, sep = ',', header = T, row.names = 1, check.names = F)
    return(data)
}

write.CSV <- function(data, filename)
{
    write.csv(data, filename, quote = F, row.names = T)
}

# ------------------------------------------------------------
# bed
# ------------------------------------------------------------
read.bed <- function(filename)
{
    d <- read.table(filename, check.names = F, header = F, sep = "\t")
    colnames(d)[1] <- 'chr'
    colnames(d)[2] <- 'start'
    colnames(d)[3] <- 'end'

    if (ncol(d) >= 4) #bed with name (ex: merge from compare_bed)
    {
        colnames(d)[4] <- 'name'
    }
    
    if (ncol(d) >= 6) #extended bed
    {
        colnames(d)[4] <- 'id'
        colnames(d)[5] <- 'extra'
        colnames(d)[6] <- 'strand'
    }
    d
}

# ------------------------------------------------------------
# peaks in MACS2 xls format
# ------------------------------------------------------------
read.peaks.xls <- function(filename)
{
    read.table(filename, header = T, check.names = F, stringsAsFactors = F)
}

# ------------------------------------------------------------
# bsf
# ------------------------------------------------------------
bsf <- function(data, cid = NA, sp = FALSE)
{
    names <- if (is.na(cid)) paste('P', rownames(data), sep = '') else data[,cid]
    if (!sp)
    {
        return(data.frame(GeneID = names, Chr = data[,'chr'], 
                      Start = data[,'start'], End = data[,'end'], Strand = 'NA'))
    }
    else
    {
        return(data.frame(GeneID = names, Chr = data[,'chr'], 
                      Start = data[,'start'], End = data[,'end'], Strand = data[,'strand']))
    }
}

# ------------------------------------------------------------
# feature count
# ------------------------------------------------------------
# reads, read.peaks.xls(peaks), read.bed(geneBody), cid=NA, nthreads=8
features.raw.fpkm <- function(bams, features, normFeatures, cid=NA, nthreads=8)
{
    ## /!\ fc parameters added on the 20170530 to include multimapping and feature overlapping reads (fractional counting) ##
    # process annotation
    ft <- bsf(normFeatures, 5, sp = T) #[1:10,]
    # get count of reads by genes
    fn <- fc(bams, annot.ext = ft, useMetaFeatures = T, ignoreDup = T, countMultiMappingReads = T, allowMultiOverlap = T, nthreads = nthreads) #-M -O (see IP-EdgeR for more details)
    # process peaks list
    ft <- bsf(features, cid, sp = F) #[1:10,]
    # get count of reads by peaks
    fr <- fc(bams, annot.ext = ft, useMetaFeatures = F, ignoreDup = T, countMultiMappingReads = T, allowMultiOverlap = T, nthreads = nthreads) #-M -O (see IP-EdgeR for more details)
    # convert raw count to fpkm 
    result <- t(apply(1000 * fr$counts / fr$annotation$Length, 1, '/', apply(fn$counts, 2, sum)) * 1e6)
    # format result
    result <- as.data.frame( result )
    common_string <- getLongestCommonSubstring(colnames( result ))
    colnames( result ) <- paste0( gsub(common_string ,"",colnames( result )), "_fpkm" )
    result <- data.frame( fr$annotation[,c("Chr","Start","End","GeneID","Length", "Strand")], result[fr$annotation$GeneID,])
    rownames( result ) <- result$GeneID
    return( result )
}

fpkm.to.logfpkm <- function(fpkm)
{
    return(apply(fpkm, c(1, 2), function(x) {ifelse(log2(x) < 0, 0, log2(x))}))
}

# ------------------------------------------------------------
# Differential table (first vs second)
# ------------------------------------------------------------
differential.table <- function(fp, s1, s2, i1, i2)
{
    d <- data.frame(row.names = rownames(fp), 
                S1 = fp[,s1],
                S2 = fp[,s2],
                I1 = fp[,i1],
                I2 = fp[,i2],
                LOR =  (fp[,s1] - fp[,i1]) - (fp[,s2] - fp[,i2]), 
                LR = fp[,s1] - fp[,s2], 
                IR = fp[,i1] - fp[,i2]
                )
    colnames(d)[1] <- s1
    colnames(d)[2] <- s2
    colnames(d)[3] <- i1
    colnames(d)[4] <- i2
    d
}

