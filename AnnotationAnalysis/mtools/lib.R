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

