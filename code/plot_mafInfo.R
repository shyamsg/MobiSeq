# Plot maf info, such as number of samples covered, and maf distribution

setwd("~/Work/MobiSeq/mafs")

library(ggplot2)
LINE.inds = read.table("LINE.allSamples.90pct.inds", as.is=T)$V1
LINE.mafs = read.table("LINE.allSamples.90pct.mafs", as.is=T)$V1

SINE.inds = read.table("SINE.allSamples.90pct.inds", as.is=T)$V1
SINE.mafs = read.table("SINE.allSamples.90pct.mafs", as.is=T)$V1

BOV2A.inds = read.table("BOV2A.allSamples.90pct.inds", as.is=T)$V1
BOV2A.mafs = read.table("BOV2A.allSamples.90pct.mafs", as.is=T)$V1

BOV2A.oc.inds = read.table("BOV2A.onlyCE.90pct.inds", as.is=T)$V1
BOV2A.oc.mafs = read.table("BOV2A.onlyCE.90pct.mafs", as.is=T)$V1

L1.inds = read.table("L1.allSamples.90pct.inds", as.is=T)$V1
L1.mafs = read.table("L1.allSamples.90pct.mafs", as.is=T)$V1

layout(matrix(1:10, nr=5, nc=2, byrow=T))
par(oma=c(0,6,4,0), xpd=F)
par(mar=c(2,2,0,0.5)+0.1)
par(mgp=c(3,1,0))

a=hist(LINE.inds, breaks=seq(4.9,10.9), plot=F)
barplot(a$counts, space=0, width=1, col="firebrick", 
     main="", xlab="", ylab="")
axis(1, at=seq(0,5)+.5, labels=seq(5,10), tick=F,mgp=c(3,.2,0))
mtext(text = "LINE", side=2, line = 3, las=2)
mtext(text="Number of individuals with\ngreater than 3 reads at SNP", side=3, line=1, font=2)
plot(density(LINE.mafs, bw=0.02), xlim=c(0,0.5), main="", xlab="", ylab="", axes=F, frame.plot=T)
axis(2)
axis(1, mgp=c(1,0.2,0), tick=F)
mtext(text="Density of minor allele frequencies\n(Bandwidth = 0.02)", side=3, line=1, font=2)

a=hist(SINE.inds, breaks=seq(4.9,10.9), plot=F)
barplot(a$counts, space=0, width=1, col="goldenrod3", 
        main="", xlab="", ylab="")
axis(1, at=seq(0,5)+.5, labels=seq(5,10), tick=F,mgp=c(1,.2,0))
mtext(text = "SINE", side=2, line = 3, las=2)
plot(density(SINE.mafs, bw=0.02), xlim=c(0,0.5), main="", xlab="", ylab="", axes=F, frame.plot=T)
axis(2)
axis(1, mgp=c(1,0.2,0), tick=F)

a=hist(BOV2A.inds, breaks=seq(13.9,28.9), plot=F)
barplot(a$counts, space=0, width=1, col="dodgerblue", 
        main="", xlab="", ylab="")
axis(1, at=seq(0,14)+.5, labels=seq(14,28), tick=F,mgp=c(1,.2,0))
mtext(text = "BOV2A", side=2, line = 3, las=2)
plot(density(BOV2A.mafs, bw=0.02), xlim=c(0,0.5), main="", xlab="", ylab="", axes=F, frame.plot=T)
axis(2)
axis(1, mgp=c(1,0.2,0), tick=F)

a=hist(BOV2A.oc.inds, breaks=seq(12.9,26.9), plot=F)
barplot(a$counts, space=0, width=1, col="dodgerblue", 
        main="", xlab="", ylab="")
axis(1, at=seq(0,13)+.5, labels=seq(13,26), tick=F,mgp=c(1,.2,0))
mtext(text = "BOV2A\nonly\nCervus\nElaphus", side=2, line = 3, las=2)
plot(density(BOV2A.oc.mafs, bw=0.02), xlim=c(0,0.5), main="", xlab="", ylab="", axes=F, frame.plot=T)
axis(2)
axis(1, mgp=c(1,0.2,0), tick=F)

a=hist(L1.inds, breaks=seq(1.9,4.9), plot=F)
barplot(a$counts, space=0, width=1, col="forestgreen", 
        main="", xlab="", ylab="")
axis(1, at=seq(0,2)+.5, labels=seq(2,4), tick=F,mgp=c(1,.2,0))
mtext(text = "L1", side=2, line = 3, las=2)
plot(density(L1.mafs, bw=0.02), xlim=c(0,0.5), main="", xlab="", ylab="", axes=F, frame.plot=T)
axis(2)
axis(1, mgp=c(1,0.2,0), tick=F)
