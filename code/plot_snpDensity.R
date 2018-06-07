# Plot the density of snps moving away from the primer site

setwd("~/Work/MobiSeq/mafs")
library(ggplot2)

layout(matrix(1:4, 2,2, byrow = T))
## LINE first as usual
line.snps = read.table("LINE.allSamples.90pct.snps.txt", as.is=T)
line.den  = density(line.snps$V3)
plot(line.den, xlab="Distance from primer site", ylab="Density of SNPs", 
     main="Wolf/LINE", sub=paste("Bandwidth =", round(line.den$bw, 2)))

## SINE
sine.snps = read.table("SINE.allSamples.90pct.snps.txt", as.is=T)
sine.den  = density(sine.snps$V3)
plot(sine.den, xlab="Distance from primer site", ylab="Density of SNPs", 
     main="Wolf/SINE", sub=paste("Bandwidth =", round(sine.den$bw, 2)))

## SINE
deer.snps = read.table("BOV2A.onlyCE.90pct.snps.txt", as.is=T)
deer.den  = density(deer.snps$V3)
plot(deer.den, xlab="Distance from primer site", ylab="Density of SNPs", 
     main="Deer/BOV2A", sub=paste("Bandwidth =", round(deer.den$bw, 2)))

## SINE
rat.snps = read.table("L1.allSamples.90pct.snps.txt", as.is=T)
rat.den  = density(rat.snps$V3)
plot(rat.den, xlab="Distance from primer site", ylab="Density of SNPs", 
     main="Rat/L1", sub=paste("Bandwidth =", round(rat.den$bw, 2)))

layout(matrix(1:4, 2,2, byrow = T))
par(lwd=0.5)
hist(line.snps$V3, xlab="Distance from primer site", ylab="Number of SNPs", 
     main="Wolf/LINE", breaks=50, col = "#beaed4")
hist(sine.snps$V3, xlab="Distance from primer site", ylab="Number of SNPs", 
     main="Wolf/SINE", breaks=50, col="#fdc086")
hist(deer.snps$V3, xlab="Distance from primer site", ylab="Number of SNPs", 
     main="Deer/BOV2A", breaks=50, col="#7fc97f")
hist(rat.snps$V3, xlab="Distance from primer site", ylab="Number of SNPs", 
     main="Rat/L1", breaks=50, col="#ffff99")
par(lwd=1)

