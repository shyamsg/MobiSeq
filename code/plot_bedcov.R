## TESeq bedcov plots
library(ggplot2)
library(GGally)

setwd("~/Work/MobiSeq/bedcovs")

# All samples
## LINE data
line.bedcov = read.table("LINE.allSamples.bedcov", as.is=T)
colnames(line.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("L", c(10,1:9)))
line.covs = line.bedcov[,-c(1:6)]
apply(line.covs, 2, function(x) summary(x) )
apply(line.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## SINE data
sine.bedcov = read.table("SINE.allSamples.bedcov", as.is=T)
colnames(sine.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("S", c(10,1:9)))
sine.covs = sine.bedcov[,-c(1:6)]
apply(sine.covs, 2, function(x) summary(x) )
apply(sine.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## bov2a data
bov2a.bedcov = read.table("BOV2A.allSamples.bedcov", as.is=T)
colnames(bov2a.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", "CE1", "DD1", "DD2",
                          "MM313", "MM315", "MM335", "MM722",
                          "MM762", "MM763", "MM764", "S10", "S11", 
                          "S12", "S13", "S14", "S15", "S16", "S17",
                          "S18", "S1", "S2", "S3", "S4", "S5",
                          "S6", "S7", "S8", "S9")
bov2a.covs = bov2a.bedcov[,-c(1:6)]
apply(bov2a.covs, 2, function(x) summary(x) )
apply(bov2a.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## bov2a data
bov2a.oc.bedcov = read.table("BOV2A.onlyCE.bedcov", as.is=T)
colnames(bov2a.oc.bedcov) = c("Scaffold", "Start", "End", "Name", 
                           "Dummy", "Strand", "CE1",
                           "MM313", "MM315", "MM335", "MM722",
                           "MM762", "MM763", "MM764", "S10", "S11", 
                           "S12", "S13", "S14", "S15", "S16", "S17",
                           "S18", "S1", "S2", "S3", "S4", "S5",
                           "S6", "S7", "S8", "S9")
bov2a.oc.covs = bov2a.oc.bedcov[,-c(1:6)]
apply(bov2a.oc.covs, 2, function(x) summary(x) )
apply(bov2a.oc.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## rat data 
L1.bedcov = read.table("L1.allSamples.bedcov", as.is=T)
colnames(L1.bedcov) = c("Scaffold", "Start", "End", "Name", 
                        "Dummy", "Strand", "Rat1", "Rat2", 
                        "Rat3", "Rat4")
L1.covs = L1.bedcov[,-c(1:6)]
apply(L1.covs, 2, function(x) summary(x) )
apply(L1.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

jpeg(filename="LINE.allSamples.jpg",width = 1200, height = 1200)
ggpairs(line.covs)
dev.off()
jpeg(filename="SINE.allSamples.jpg",width = 1200, height = 1200)
ggpairs(sine.covs)
dev.off()
jpeg(filename="BOV2A.allSamples.jpg",width = 2000, height = 2000)
ggpairs(bov2a.covs)
dev.off()
jpeg(filename="BOV2A.onlyCE.jpg",width = 2000, height = 2000)
ggpairs(bov2a.covs)
dev.off()
jpeg(filename="L1.allSamples.jpg",width = 1200, height = 1200)
ggpairs(L1.covs)
dev.off()


# 90 pct
## LINE data
line.bedcov = read.table("LINE.allSamples.90pct.bedcov", as.is=T)
colnames(line.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("L", c(10,1:9)))
line.covs = line.bedcov[,-c(1:6)]
apply(line.covs, 2, function(x) summary(x) )
apply(line.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## SINE data
sine.bedcov = read.table("SINE.allSamples.90pct.bedcov", as.is=T)
colnames(sine.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("S", c(10,1:9)))
sine.covs = sine.bedcov[,-c(1:6)]
apply(sine.covs, 2, function(x) summary(x) )
apply(sine.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## bov2a data
bov2a.bedcov = read.table("BOV2A.allSamples.90pct.bedcov", as.is=T)
colnames(bov2a.bedcov) = c("Scaffold", "Start", "End", "Name", 
                           "Dummy", "Strand", "CE1", "DD1", "DD2",
                           "MM313", "MM315", "MM335", "MM722",
                           "MM762", "MM763", "MM764", "S10", "S11", 
                           "S12", "S13", "S14", "S15", "S16", "S17",
                           "S18", "S1", "S2", "S3", "S4", "S5",
                           "S6", "S7", "S8", "S9")
bov2a.covs = bov2a.bedcov[,-c(1:6)]
apply(bov2a.covs, 2, function(x) summary(x) )
apply(bov2a.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## bov2a data
bov2a.oc.bedcov = read.table("BOV2A.onlyCE.90pct.bedcov", as.is=T)
colnames(bov2a.oc.bedcov) = c("Scaffold", "Start", "End", "Name", 
                              "Dummy", "Strand", "CE1",
                              "MM313", "MM315", "MM335", "MM722",
                              "MM762", "MM763", "MM764", "S10", "S11", 
                              "S12", "S13", "S14", "S15", "S16", "S17",
                              "S18", "S1", "S2", "S3", "S4", "S5",
                              "S6", "S7", "S8", "S9")
bov2a.oc.covs = bov2a.oc.bedcov[,-c(1:6)]
apply(bov2a.oc.covs, 2, function(x) summary(x) )
apply(bov2a.oc.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## rat data 
L1.bedcov = read.table("L1.allSamples.90pct.bedcov", as.is=T)
colnames(L1.bedcov) = c("Scaffold", "Start", "End", "Name", 
                        "Dummy", "Strand", "Rat1", "Rat2", 
                        "Rat3", "Rat4")
L1.covs = L1.bedcov[,-c(1:6)]
apply(L1.covs, 2, function(x) summary(x) )
apply(L1.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

jpeg(filename="LINE.allSamples.90pct.jpg",width = 1200, height = 1200)
ggpairs(line.covs)
dev.off()
jpeg(filename="SINE.allSamples.90pct.jpg",width = 1200, height = 1200)
ggpairs(sine.covs)
dev.off()
jpeg(filename="BOV2A.allSamples.90pct.jpg",width = 2000, height = 2000)
ggpairs(bov2a.covs)
dev.off()
jpeg(filename="BOV2A.onlyCE.90pct.jpg",width = 2000, height = 2000)
ggpairs(bov2a.covs)
dev.off()
jpeg(filename="L1.allSamples.90pct.jpg",width = 1200, height = 1200)
ggpairs(L1.covs)
dev.off()
