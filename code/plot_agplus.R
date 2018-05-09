## plot agplus
setwd("~/Work/MobiSeq/aggPlot/")

##LINE
library(ggplot2)
library(tidyr)
library(dplyr)


movingSum = function(counts, window=250) {
  nr = nrow(counts)
  nc = ncol(counts)
  positions = counts[1:(nr-window),1]
  fullcounts = apply(counts[,-1], 2, function(x) {
    temp = rep(0,nr-window)
    for (index in 1:(nr-window)) {
      temp[index] = sum(x[index:(index+window)])
    }
    return(temp)
  })
  return(cbind(positions, fullcounts))
}

inLine = paste0("L", 1:10, ".Wolf_noHets.90pct.nodupsec.agplus.txt")
line.combined = read.table(inLine[1], header=T)[,1]
for (infile in inLine) {
  ff = read.table(infile,header=T)
  line.combined = cbind(line.combined, ff[,2])
}
colnames(line.combined) = c("distance", paste0("L", 1:10))

line.combined = movingSum(line.combined)

line.combined.tib = as_tibble(line.combined)
line.combined.tib = gather(line.combined.tib, sample, coverage, -positions)
pdf(file="LINE.coverage.pdf")
ggplot(line.combined.tib, aes(x=positions, y=coverage, group=sample, col=sample)) + geom_line() +
  xlab("Distance from LINE primer") + ylab ("Coverage per million reads") + xlim(-750,750)
dev.off()

inSine = paste0("S", 1:10, ".Wolf_noHets.90pct.nodupsec.agplus.txt")
sine.combined = read.table(inSine[1], header=T)[,1]
for (infile in inSine) {
  ff = read.table(infile,header=T)
  sine.combined = cbind(sine.combined, ff[,2])
}
colnames(sine.combined) = c("distance", paste0("S", 1:10))

sine.combined = movingSum(sine.combined)

sine.combined.tib = as_tibble(sine.combined)
sine.combined.tib = gather(sine.combined.tib, sample, coverage, -positions)
pdf(file="SINE.coverage.pdf")
ggplot(sine.combined.tib, aes(x=positions, y=coverage, group=sample, col=sample)) + geom_line() +
  xlab("Distance from SINE primer") + ylab ("Coverage per million reads") + xlim(-750,750)
dev.off()

inDeer = paste0(c("CE1", "DD1", "DD2",
                  "MM313", "MM315", "MM335", "MM722",
                  "MM762", "MM763", "MM764", "S10", "S11",
                  "S12", "S13", "S14", "S15", "S16", "S17",
                  "S18", "S1", "S2", "S3", "S4", "S5",
                  "S6", "S7", "S8", "S9"), ".CervusElaphus.allSamples.90pct.nodupsec.agplus.txt")
deer.combined = read.table(inDeer[1], header=T)[,1]
for (infile in inDeer) {
  ff = read.table(infile,header=T)
  deer.combined = cbind(deer.combined, ff[,2])
}
colnames(deer.combined) = c("distance", "CE1", "DD1", "DD2",
                            "MM313", "MM315", "MM335", "MM722",
                            "MM762", "MM763", "MM764", "S10", "S11",
                            "S12", "S13", "S14", "S15", "S16", "S17",
                            "S18", "S1", "S2", "S3", "S4", "S5",
                            "S6", "S7", "S8", "S9")

deer.combined = movingSum(deer.combined)

deer.combined.tib = as_tibble(deer.combined)
deer.combined.tib = gather(deer.combined.tib, sample, coverage, -positions)
pdf(file="BOV2A.allSamples.coverage.pdf")
ggplot(deer.combined.tib, aes(x=positions, y=coverage, group=sample, col=sample)) + geom_line() +
  xlab("Distance from BOV2A primer") + ylab ("Coverage per million reads") + xlim(-750, 750)
dev.off()

inDeer = paste0(c("CE1", "MM313", "MM315", "MM335", "MM722",
                  "MM762", "MM763", "MM764", "S10", "S11",
                  "S12", "S13", "S14", "S15", "S16", "S17",
                  "S18", "S1", "S2", "S3", "S4", "S5",
                  "S6", "S7", "S8", "S9"), ".CervusElaphus.onlyCE.90pct.nodupsec.agplus.txt")
deer.combined = read.table(inDeer[1], header=T)[,1]
for (infile in inDeer) {
  ff = read.table(infile,header=T)
  deer.combined = cbind(deer.combined, ff[,2])
}
colnames(deer.combined) = c("distance", "CE1",
                            "MM313", "MM315", "MM335", "MM722",
                            "MM762", "MM763", "MM764", "S10", "S11",
                            "S12", "S13", "S14", "S15", "S16", "S17",
                            "S18", "S1", "S2", "S3", "S4", "S5",
                            "S6", "S7", "S8", "S9")

deer.combined = movingSum(deer.combined)

deer.combined.tib = as_tibble(deer.combined)
deer.combined.tib = gather(deer.combined.tib, sample, coverage, -positions)
pdf(file="BOV2A.onlyCE.coverage.pdf")
ggplot(deer.combined.tib, aes(x=positions, y=coverage, group=sample, col=sample)) + geom_line() +
  xlab("Distance from BOV2A primer") + ylab ("Coverage per million reads") + xlim(-750, 750)
dev.off()


inRat = paste0("Rat", 1:4, ".rn6.90pct.nodupsec.agplus.txt")
rat.combined = read.table(inRat[1], header=T)[,1]
for (infile in inRat) {
  ff = read.table(infile,header=T)
  rat.combined = cbind(rat.combined, ff[,2])
}
colnames(rat.combined) = c("distance", paste0("Rat", 1:4))

rat.combined = movingSum(rat.combined)

rat.combined = as_tibble(rat.combined)
rat.combined = gather(rat.combined, sample, coverage, -positions)
pdf(file="L1.coverage.pdf")
ggplot(rat.combined, aes(x=positions, y=coverage, group=sample, col=sample)) + geom_line() +
  xlab("Distance from L1 primer") + ylab ("Coverage per million reads") + xlim(-750,750)
dev.off()