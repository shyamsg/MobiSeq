## Preseq plots
setwd("~/Work/MobiSeq/preseq")

library(ggplot2)
library(tidyr)
library(dplyr)

inLine = paste0("L", 1:10, ".Wolf_noHets.allReads.lc")
line.combined = read.table(inLine[1], header=T)[,1]
for (infile in inLine) {
  ff = read.table(infile,header=T)
  line.combined = cbind(line.combined, ff[,2])
}
colnames(line.combined) = c("seqDepth", paste0("L", 1:10))

line.combined = as_tibble(line.combined)
line.combined = gather(line.combined, sample, uniqReads, -seqDepth)
ggplot(line.combined, aes(x=seqDepth, y=uniqReads, group=sample, col=sample)) + geom_line() +
  xlab("Total sequenced reads") + ylab ("Unique reads") + xlim(0,2.5e8)


inSine = paste0("S", 1:10, ".Wolf_noHets.allReads.lc")
sine.combined = read.table(inSine[1], header=T)[,1]
for (infile in inSine) {
  ff = read.table(infile,header=T)
  sine.combined = cbind(sine.combined, ff[,2])
}
colnames(sine.combined) = c("seqDepth", paste0("S", 1:10))

sine.combined = as_tibble(sine.combined)
sine.combined = gather(sine.combined, sample, uniqReads, -seqDepth)
ggplot(sine.combined, aes(x=seqDepth, y=uniqReads, group=sample, col=sample)) + geom_line() +
  xlab("Total sequenced reads") + ylab ("Unique reads") + xlim(0,5e8)

inDeer = paste0(c("CE1", "DD1", "DD2",
                  "MM313", "MM315", "MM335", "MM722",
                  "MM762", "MM763", "MM764", "S10", "S11", 
                  "S12", "S13", "S14", "S15", "S16", "S17",
                  "S18", "S1", "S2", "S3", "S4", "S5",
                  "S6", "S7", "S8", "S9"), ".CervusElaphus.allReads.lc")
deer.combined = read.table(inDeer[1], header=T)[,1]
for (infile in inDeer) {
  ff = read.table(infile,header=T)
  deer.combined = cbind(deer.combined, ff[,2])
}
colnames(deer.combined) = c("seqDepth", "CE1", "DD1", "DD2",
                            "MM313", "MM315", "MM335", "MM722",
                            "MM762", "MM763", "MM764", "S10", "S11", 
                            "S12", "S13", "S14", "S15", "S16", "S17",
                            "S18", "S1", "S2", "S3", "S4", "S5",
                            "S6", "S7", "S8", "S9")

deer.combined = as_tibble(deer.combined)
deer.combined = gather(deer.combined, sample, uniqReads, -seqDepth)
ggplot(deer.combined, aes(x=seqDepth, y=uniqReads, group=sample, col=sample)) + geom_line() +
  xlab("Total sequenced reads") + ylab ("Unique reads")


inRat = c()
for (rat in 1:4) {
  for (tissue in c(1:3,5:14)) {
    inRat = c(inRat, paste0("Rat", rat, ".Tissue", tissue, ".lc"))
  }
}
rat.combined = read.table(inRat[1], header=T)[,1]
for (infile in inRat) {
  ff = read.table(infile,header=T)
  rat.combined = cbind(rat.combined, ff[,2])
}
colnames(rat.combined) = c("seqDepth", sub(".lc", "", inRat))

rat.combined = as_tibble(rat.combined)
rat.combined = gather(rat.combined, sample, uniqReads, -seqDepth)
ggplot(rat.combined, aes(x=seqDepth, y=uniqReads, group=sample, col=sample)) + geom_line() +
  xlab("Total sequenced reads") + ylab ("Unique reads")



