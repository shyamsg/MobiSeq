# Plot the density of snps moving away from the primer site

setwd("~/Work/MobiSeq/mafs")
library(ggplot2)
library(gridExtra)

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


## read the fai files
wolf.chrlengths = read.table("~/Work/MobiSeq/beds/wolf_chrlengths.genome", as.is=T)
deer.chrlengths = read.table("~/Work/MobiSeq/beds/deer_chrlengths.genome", as.is=T)
rats.chrlengths = read.table("~/Work/MobiSeq/beds/rn6_chrlengths.genome", as.is=T)

## LINE
line.primers = read.table("~/Work/MobiSeq/matches/wolf/LINE/LINE.allSamples.90pct.bed", as.is=T)$V1
line.snp.primers = data.frame(chr=wolf.chrlengths$V1)
snps = c()
primers = c()
for (index in 1:nrow(wolf.chrlengths)) {
  snps    = c(snps, sum(line.snps$V1 == line.snp.primers$chr[index]))
  primers = c(primers, sum(line.primers == line.snp.primers$chr[index]))
}
line.snp.primers$snps    = cumsum(as.numeric(snps))/cumsum(as.numeric(wolf.chrlengths$V2))
line.snp.primers$primers = cumsum(as.numeric(primers))/cumsum(as.numeric(wolf.chrlengths$V2))
#line.snp.primers         = line.snp.primers[order(wolf.chrlengths$V2, decreasing=T),]

line.1 = ggplot(line.snp.primers, aes(x=1:nrow(line.snp.primers), y=snps*1e6)) + 
         geom_point(size=3, col="#beaed488", shape=19) +
         xlab("Scaffold") + ylab("Cumulative snps/Mb") + 
         theme_bw()
line.2 = ggplot(line.snp.primers, aes(x=1:nrow(line.snp.primers), y=primers*1e6)) + 
         geom_point(size=2, col="#beaed488", shape=19) +
         xlab("Scaffold") + ylab("Cumulative loci/Mb") +
         theme_bw()

## SINE
sine.primers = read.table("~/Work/MobiSeq/matches/wolf/SINE/SINE.allSamples.90pct.bed", as.is=T)$V1
sine.snp.primers = data.frame(chr=wolf.chrlengths$V1)
snps = c()
primers = c()
for (index in 1:nrow(wolf.chrlengths)) {
  snps    = c(snps, sum(sine.snps$V1 == sine.snp.primers$chr[index]))
  primers = c(primers, sum(sine.primers == sine.snp.primers$chr[index]))
}
sine.snp.primers$snps    = cumsum(as.numeric(snps))/cumsum(as.numeric(wolf.chrlengths$V2))
sine.snp.primers$primers = cumsum(as.numeric(primers))/cumsum(as.numeric(wolf.chrlengths$V2))
#sine.snp.primers         = sine.snp.primers[order(wolf.chrlengths$V2, decreasing=T),]

sine.1 = ggplot(sine.snp.primers, aes(x=1:nrow(sine.snp.primers), y=snps*1e6)) + 
         geom_point(size=3, col="#fdc08666", shape=19) +
         xlab("Scaffold") + ylab("Cumulative snps/Mb") + 
         theme_bw()
sine.2 = ggplot(sine.snp.primers, aes(x=1:nrow(sine.snp.primers), y=primers*1e6)) + 
         geom_point(size=2, col="#fdc08666", shape=19) +
         xlab("Scaffold") + ylab("Cumulative loci/Mb") +
         theme_bw()

## deer
deer.primers = read.table("~/Work/MobiSeq/matches/deer/BOV2A.onlyCE.90pct.bed", as.is=T)$V1
deer.snp.primers = data.frame(chr=deer.chrlengths$V1)
snps = rep(0, nrow(deer.snp.primers))
temp = tapply(deer.snps$V1, deer.snps$V1, length)
tempname = sapply(names(temp), function(x) {
  which(deer.snp.primers$chr == x)
})
snps[tempname] = temp
primers = rep(0, nrow(deer.snp.primers))
temp = tapply(deer.primers, deer.primers, length)
tempname = sapply(names(temp), function(x) {
  which(deer.snp.primers$chr == x)
})
primers[tempname] = temp
deer.snp.primers$snps    = cumsum(as.numeric(snps))/cumsum(as.numeric(deer.chrlengths$V2))
deer.snp.primers$primers = cumsum(as.numeric(primers))/cumsum(as.numeric(deer.chrlengths$V2))
#deer.snp.primers         = deer.snp.primers[order(deer.chrlengths$V2, decreasing=T),]

deer.1 = ggplot(deer.snp.primers, aes(x=1:nrow(deer.snp.primers), y=snps*1e6)) + 
        geom_point(size=3, col="#7fc97f66", shape=19) +
        xlab("Scaffold") + ylab("Cumulative snps/Mb") + 
        theme_bw()
deer.2 = ggplot(deer.snp.primers, aes(x=1:nrow(deer.snp.primers), y=primers*1e6)) + 
        geom_point(size=2, col="#7fc97f66", shape=19) +
        xlab("Scaffold") + ylab("Cumulative loci/Mb") +
        theme_bw()

## rats
rats.primers = read.table("~/Work/MobiSeq/matches/rats/L1.allSamples.90pct.bed", as.is=T)$V1
rats.snp.primers = data.frame(chr=rats.chrlengths$V1)
snps = c()
primers = c()
for (index in 1:nrow(rats.chrlengths)) {
  snps    = c(snps, sum(rat.snps$V1 == rats.snp.primers$chr[index]))
  primers = c(primers, sum(rats.primers == rats.snp.primers$chr[index]))
}
rats.snp.primers$snps    = cumsum(as.numeric(snps))/cumsum(as.numeric(rats.chrlengths$V2))
rats.snp.primers$primers = cumsum(as.numeric(primers))/cumsum(as.numeric(rats.chrlengths$V2))
#rats.snp.primers         = rats.snp.primers[order(rats.chrlengths$V2, decreasing=T),]

rats.1 = ggplot(rats.snp.primers, aes(x=1:nrow(rats.snp.primers), y=snps*1e6)) + 
         geom_point(size=3, col="#1e90ff66", shape=19) +
         xlab("Scaffold") + ylab("Cumulative snps/Mb") + 
         theme_bw()
rats.2 = ggplot(rats.snp.primers, aes(x=1:nrow(rats.snp.primers), y=primers*1e6)) + 
         geom_point(size=2, col="#1e90ff66", shape=19) +
         xlab("Scaffold") + ylab("Cumulative loci/Mb") +
         theme_bw()

#Plot
grid.arrange(line.1, line.2, sine.1, sine.2, deer.1, deer.2, rats.1, rats.2,
             layout_matrix = matrix(1:8, nr=4, nc=2, byrow=TRUE))

