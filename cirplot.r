# args = commandArgs(trailingOnly=TRUE)

new_dat <- read.csv("/Users/daniel/Desktop/final_10000bps.igv", sep="\t", header=FALSE)
n_dat <- new_dat[,c(1,2,5)]

# change from "0.525000" to "52", yes it ROUNDS DOWN
n_dat$V5 <- as.numeric(as.character(sub("0." , "", n_dat$V5)))
n_dat$V5 <- substr(n_dat$V5, 0, 2)
n_dat$V5 <- formatC(as.numeric(n_dat$V5),width=2,format='f',digits=0,flag='0')

circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
{
  curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
  curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
}

# make as matrix
zzz <- as.matrix(n_dat)

# check for NA and 0's
zzz <- zzz[zzz[, 1] != 0, ]
max.chr <- max(as.numeric(zzz[, 1]), na.rm=TRUE)
if(is.infinite(max.chr))  max.chr <- 0
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
chr.labels <- unique(n_dat$V1)

Num <- as.numeric(table(n_dat[,1]))
Nchr <- length(Num)
N <- NULL
pvalueT <- as.matrix(n_dat[,3])
pvalue.pos <- n_dat[, 2]
pvalue.posN <- NULL
ticks <- NULL
pvalue.pos.list <- tapply(pvalue.pos, n_dat[, 1], list)

if (!missing(band)) {
  band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
} else {
  band <- floor((sum(sapply(pvalue.pos.list, max))/100))
}
if (band==0)	band=1

if (LOG10) {
  pvalueT[pvalueT <= 0] <- 1
  pvalueT[pvalueT > 1] <- 1
}

for(i in 0:(Nchr-1)){
  if (i==0){
    #pvalue <- append(pvalue,rep(Inf,band),after=0)
    pvalue.posN <- pvalue.pos.list[[i+1]] + band
    ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
  }else{
    #pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
    pvalue.posN <- c(pvalue.posN, max(pvalue.posN) + band + pvalue.pos.list[[i+1]])
    ticks[i+1] <- max(pvalue.posN)-floor(max(pvalue.pos.list[[i+1]])/2)
  }
} 
pvalue.posN.list <- tapply(pvalue.posN, n_dat[, 1], list)
TotalN <- max(pvalue.posN)

# Assign Variables for plotting
plot <- TRUE
band <- 1
main <- "GC Content Percentage"
maxbin.num <- NULL
bin <- 10000
R <- 1 # this will be changed to 2 when I plot repeats (masking n's)

#Assuming %, this is for GC
legend.max <- 100
#Change colors here
col <- c("black", "red", "orange", "yellow", "darkgreen", "grey4")
col.seg <- NULL
width <- 5
legend.len <- 100
legend.cex <- 4
cex.axis <- 1
legend.y.intersp <- 1
legend.x.intersp <- 1
legend.pt.cex <- 3
file.output <- TRUE
cir.chr=TRUE
cir.chr.h= 8.5
chr.den.col = c("black", "red", "yellow", "darkgreen", "grey4")
cir.legend=TRUE
cir.legend.cex = 8.6
cir.legend.col="black"
LOG10=TRUE
box=FALSE
file="jpg"    #output type -> "jpg", "png", "tiff"
dpi=300
r = 15.5
H = 1
cir.band= 2
Max <- 5

pos.x <- list()
col.index <- list()

for(i in 1 : length(chr.num)) {
  pos.x[[i]] <- pos[which(chr == chr.num[i])]
  col.index[[i]] <- gc[which(chr == chr.num[i])]
}

# I want this 100 so it is percent by definition
col=colorRampPalette(col)(100)

# THIS IS IT
for (i in 1 : length(chr.num)) {
  col.seg <- c(col.seg, col[col.index[[i]]])
}

if (length(zzz.xy.index) != 0) {
  for (i in 1:length(chr.xy)) {
    chr.num[chr.num == max.chr + i] <- chr.xy[i]
  }
}

chr.num <- rev(chr.num)

if (chorm.maxlen > 10000000) {
  top_axis_length <- 10
  kilo_or_mega <- 1e6    # mega
  text_size <- "Mb"
} else {
  top_axis_length <- 5
  kilo_or_mega <- 1e3  # kilo
  text_size <- "Kb"
}

# Make legend not so large, yet informative
legend.y <- seq(10, 90, 5)
legend.y.col <- legend.y

# pick every 5 colors from palette
legend.col <- col[seq(10,90,5)]

if (length(chr.den.col) > 1) {
  cir.density=TRUE
  den.fold <- 20
  density.list <- list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
} else {
  cir.density=FALSE
}
if (!file.output) {
  if (!is.null(dev.list())) dev.new(width=8, height=8)
  par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
}
par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
RR <- 2*r+H*R+cir.band*R
RI <- r-3.5+H*R+cir.band*R
if (cir.density) {
  plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
} else {
  plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
}

if(cir.chr==TRUE) {
  #plot the boundary which represents the chromosomes
  polygon.num <- 1000
  for (k in 1:length(chr)) {
    polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
    X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
    Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
    X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
    Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
    if(is.null(chr.den.col)){
      polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
    } else {
      if (cir.density) {
        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="white",border="white")
      } else {
        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
      }
    }
  }
}

if (cir.density) {
  segments(
    (1.3 * RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    col=density.list$den.col, lwd=0.8
  )
  legend(
    85,
    1,
    title="GC %", legend=density.list$legend.y, pch=15, pt.cex = 3, col = legend.col,
    cex=0.8, bty="n",
    y.intersp=1,
    x.intersp=1,
    yjust=0.5, xjust=0, xpd=TRUE
  )
}

# ticks1=1.4*(RR+cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
# ticks2=1.4*(RR+cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)

# for(i in 1:length(ticks)){
  # angle=360*(1-(ticks-round(band/2))[i]/TotalN)
  # text(ticks1[i],ticks2[i],chr.labels[i],srt=angle,font=2, cex=0.6)
# }

X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN-round(band/2))/TotalN)
Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN-round(band/2))/TotalN)
points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))

if (cir.chr==TRUE) {
  # XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
  # YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
  # lines(XLine,YLine,lwd=1.5)
  
  polygon.num <- 10000
  for (k in 1:length(chr)) {
    polygon.index <- seq(1+round(band/2)+max(pvalue.posN.list[[k-1]]),-round(band/2)+max(pvalue.posN.list[[k]]), length=polygon.num)
    X1chr=(2*cir.band+RI)*sin(4*pi*(polygon.index)/TotalN)
    Y1chr=(2*cir.band+RI)*cos(4*pi*(polygon.index)/TotalN)
    X2chr=(2*cir.band+RI+cir.chr.h)*sin(4*pi*(polygon.index)/TotalN)
    Y2chr=(2*cir.band+RI+cir.chr.h)*cos(4*pi*(polygon.index)/TotalN)
    if (is.null(chr.den.col)) {
      polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k]) 
    }
  }
}

if (cir.density) {
  segments(
    (8*cir.band+RI)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    col=density.list$den.col, lwd=0.8
  )
  
}