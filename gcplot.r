
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


new_dat <- read.csv(args[1], sep="\t", header=FALSE)
n_dat <- new_dat[,c(1,2,5)]

# change from "0.525000" to "52", yes it ROUNDS DOWN
n_dat$V5 <- as.numeric(as.character(sub("0." , "", n_dat$V5)))
n_dat$V5 <- substr(n_dat$V5, 0, 2)
n_dat$V5 <- formatC(as.numeric(n_dat$V5),width=2,format='f',digits=0,flag='0')

# make as matrix
zzz <- as.matrix(n_dat)

# check for NA and 0's
zzz <- zzz[zzz[, 1] != 0, ]
max.chr <- max(as.numeric(zzz[, 1]), na.rm=TRUE)
if(is.infinite(max.chr))	max.chr <- 0
zzz.xy.index <- which(!as.numeric(zzz[, 1]) %in% c(0 : max.chr))

if(length(zzz.xy.index) != 0){
  chr.xy <- unique(zzz[zzz.xy.index, 1])
  for(i in 1:length(chr.xy)){
    zzz[zzz[, 1] == chr.xy[i], 1] <- max.chr + i
  }
}

zzz <- zzz[order(as.numeric(zzz[, 1]), as.numeric(zzz[, 2])), ]
chr <- as.numeric(zzz[, 1])
pos <- as.numeric(zzz[, 2])
gc <- as.numeric(zzz[,3])
chr.num <- unique(chr)
chorm.maxlen <- max(pos)

# Assign Variables for plotting
plot <- TRUE
band <- 3
main <- "GC Content Percentage"
maxbin.num <- NULL
bin <- args[2]

#Assuming %, this is for GC
legend.max <- 100

#Change colors here
col <- c("grey4", "red", "yellow", "darkgreen", "black")
col.seg <- NULL
width <- 5
legend.len <- 100
legend.cex <- 1
legend.y.intersp <- 1
legend.x.intersp <- 1
legend.pt.cex <- 3
file.output <- TRUE

if(plot)	plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
pos.x <- list()
col.index <- list()

for(i in 1 : length(chr.num)) {
  pos.x[[i]] <- pos[which(chr == chr.num[i])]
  col.index[[i]] <- gc[which(chr == chr.num[i])]
}

# I want this 100 so it is percent by definition
col=colorRampPalette(col)(100)

# THIS IS IT
for(i in 1 : length(chr.num)){
  if(plot)	polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
                   c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
                     width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey", border="grey")
  col.seg <- c(col.seg, col[col.index[[i]]])
  if(plot)	segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
                    col=col[col.index[[i]]], lwd=1)
}

if(length(zzz.xy.index) != 0){
  for(i in 1:length(chr.xy)){
    chr.num[chr.num == max.chr + i] <- chr.xy[i]
  }
}

chr.num <- rev(chr.num)

if(plot)	mtext(at=seq(band, length(chr.num) * band, band),text=paste(chr.num, sep=""), side=2, las=2, font=1, cex=0.6, line=0.2)
if(plot)	axis(3, at=seq(0, chorm.maxlen, length=5), labels=c(NA, paste(round((seq(0, chorm.maxlen, length=5))[-1] / 1e3, 0), "Kb", sep="")),
              font=1, cex.axis=0.8, tck=0.01, lwd=2, padj=1.2)

# Make legend not so large, yet informative
legend.y <- seq(10, 90, 5)
legend.y.col <- legend.y

# pick every 5 colors from palette
legend.col <- col[seq(10,90,5)]

# plot legend
if(plot)	legend("topright", title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=legend.col,
                cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)


