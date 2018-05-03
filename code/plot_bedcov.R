## TESeq bedcov plots
library(ggplot2)
library(GGally)

setwd("~/Work/TESeq/bedcovs")

## LINE data
line.bedcov = read.table("wolf_LINE.bedcov", as.is=T)
colnames(line.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("L", c(10,1:9)))
line.covs = line.bedcov[,-c(1:6)]
apply(line.covs, 2, function(x) summary(x) )
apply(line.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## SINE data
sine.bedcov = read.table("wolf_SINE.bedcov", as.is=T)
colnames(sine.bedcov) = c("Scaffold", "Start", "End", "Name", 
                          "Dummy", "Strand", paste0("S", c(10,1:9)))
sine.covs = sine.bedcov[,-c(1:6)]
apply(sine.covs, 2, function(x) summary(x) )
apply(sine.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

## bov2a data
bov2a.bedcov = read.table("deer_BOV2A.bedcov", as.is=T)
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

## rat data 
L1.bedcov = read.table("rats_L1.bedcov", as.is=T)
colnames(L1.bedcov) = c("Scaffold", "Start", "End", "Name", 
                        "Dummy", "Strand", "Rat1", "Rat2", 
                        "Rat3", "Rat4")
L1.covs = L1.bedcov[,-c(1:6)]
apply(L1.covs, 2, function(x) summary(x) )
apply(L1.covs, 2, function(x) quantile(x,seq(0,1,by=0.1)) )

jpeg(filename="wolf_LINE_covar.jpg",width = 1200, height = 1200)
ggpairs(line.covs)
dev.off()
jpeg(filename="wolf_SINE_covar.jpg",width = 1200, height = 1200)
ggpairs(sine.covs)
dev.off()
jpeg(filename="deer_BOV2A_covar.jpg",width = 2000, height = 2000)
ggpairs(bov2a.covs)
dev.off()
jpeg(filename="rats_L1_covar.jpg",width = 1200, height = 1200)
ggpairs(L1.covs)
dev.off()
