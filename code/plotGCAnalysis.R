## Plot the gc content and distance from primer site

setwd("~/Work/MobiSeq/gcAnalysis")
library(ggplot2)
library(gridExtra)

## LINE details
line.gc           = read.table("LINE.gc_coverage.bed")
colnames(line.gc) = c("scaffold", "start", "end", "name", "gc", "coverage")
## SINE details
sine.gc           = read.table("SINE.gc_coverage.bed")
colnames(sine.gc) = c("scaffold", "start", "end", "name", "gc", "coverage")
## LINE details
deer.gc           = read.table("BOV2A_onlyCE.gc_coverage.bed")
colnames(deer.gc) = c("scaffold", "start", "end", "name", "gc", "coverage")
## LINE details
rats.gc           = read.table("L1.gc_coverage.bed")
colnames(rats.gc) = c("scaffold", "start", "end", "name", "gc", "coverage")

#par(mar=c(2,2,1,1)+0.1, oma=c(5,4,0,0), xpd=FALSE)
line.1 = ggplot(line.gc, aes(x=gc*100, y=coverage)) + geom_point(size=3, col="#beaed488", shape=19) +
  xlab("GC %") + ylab("Coverage") + theme_bw()
line.2 = ggplot(line.gc, aes(x=gc*100)) + geom_histogram(fill="#beaed4", col="black", size=0.25) +
  geom_line(aes(x=x, y=y*50), data=data.frame(loess.smooth(line.gc$gc*100, line.gc$coverage, span=0.05, degree=2))) + 
  scale_y_continuous(sec.axis = sec_axis(~./50, name = "Average coverage")) + xlab("GC %") + ylab("Number of loci") +
  theme_bw()

sine.1 = ggplot(sine.gc, aes(x=gc*100, y=coverage)) + geom_point(size=3, col="#fdc08688", shape=19) +
  xlab("GC %") + ylab("Coverage") + theme_bw()
sine.2 = ggplot(sine.gc, aes(x=gc*100)) + geom_histogram(fill="#fdc086", col="black", size=0.25) +
  geom_line(aes(x=x, y=y*150), data=data.frame(loess.smooth(sine.gc$gc*100, sine.gc$coverage, span=0.05, degree=2))) + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = "Average coverage")) + xlab("GC %") + ylab("Number of loci") +
  theme_bw()

deer.1 = ggplot(deer.gc, aes(x=gc*100, y=coverage)) + geom_point(size=3, col="#7fc97f88", shape=19) +
  xlab("GC %") + ylab("Coverage") + theme_bw()
deer.2 = ggplot(sine.gc, aes(x=gc*100)) + geom_histogram(fill="#7fc97f", col="black", size=0.25) +
  geom_line(aes(x=x, y=y*20), data=data.frame(loess.smooth(deer.gc$gc*100, deer.gc$coverage, span=0.05, degree=2))[1:40,]) + 
  scale_y_continuous(sec.axis = sec_axis(~./20, name = "Average coverage")) + xlab("GC %") + ylab("Number of loci") + 
  theme_bw()

rats.1 = ggplot(rats.gc, aes(x=gc*100, y=coverage)) + geom_point(size=3, col="#1e90ff88", shape=19) +
  xlab("GC %") + ylab("Coverage")+ theme_bw()
rats.2 = ggplot(rats.gc, aes(x=gc*100)) + geom_histogram(fill="#1e90ff", col="black", size=0.25) +
  geom_line(aes(x=x, y=y*200), data=data.frame(loess.smooth(rats.gc$gc*100, rats.gc$coverage, span=0.05, degree=2))) + 
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Average coverage")) + xlab("GC %") + ylab("Number of loci") + 
  theme_bw()

grid.arrange(line.1, line.2, sine.1, sine.2, deer.1, deer.2, rats.1, rats.2, 
             layout_matrix = matrix(1:8, nr=4, nc=2, byrow=TRUE))

plot(sine.gc$gc*100, sine.gc$coverage, pch=19, col="#fdc08666", 
     xlab="GC %", ylab="Average coverage")
smoothScatter(sine.gc$gc*100, sine.gc$coverage, pch=19, col="#fdc086", ylim=c(0,100), 
              colramp = colorRampPalette(c("white", "#fdc086")), cex=0.5, xlab="GC %", ylab="Average coverage")
lines(loess.smooth(sine.gc$gc*100, sine.gc$coverage, span=0.05, degree=2))

plot(deer.gc$gc*100, deer.gc$coverage, pch=19, col="#7fc97f66", 
     xlab="GC %", ylab="Average coverage")
smoothScatter(deer.gc$gc*100, deer.gc$coverage, pch=19, col="#7fc97f", ylim=c(0,100), 
              colramp = colorRampPalette(c("white", "#7fc97f")), cex=0.5, xlab="GC %", ylab="Average coverage")
lines(loess.smooth(deer.gc$gc*100, deer.gc$coverage, span=0.05, degree=2))

plot(rats.gc$gc*100, rats.gc$coverage, pch=19, col="#1e90ff66", 
     xlab="GC %", ylab="Average coverage")
smoothScatter(rats.gc$gc*100, rats.gc$coverage, pch=19, col="#1e90ff", ylim=c(0,100), 
              colramp = colorRampPalette(c("white", "#1e90ff")), cex=0.5, xlab="GC %", ylab="Average coverage")
lines(loess.smooth(rats.gc$gc*100, rats.gc$coverage, span=0.05, degree=2))

test = data.frame(loess.smooth(sine.gc$gc*100, sine.gc$coverage, span=0.05, degree=2))

a = ggplot(sine.gc, aes(x=gc*100)) + geom_histogram(fill="#fdc086", col="black", size=0.25)
a = a + geom_line(aes(x=x, y=y*150), data=test)
a = a + scale_y_continuous(sec.axis = sec_axis(~./150, name = "Average coverage"))
a = a + xlab("GC %") + ylab("Number of loci")
print(a)

