args = commandArgs(trailingOnly=TRUE)

new_dat <- read.csv(args[1], sep="\t", header=FALSE)
n_dat <- new_dat[,c(1,2,5)]
masked_dat <- read.csv(args[3], sep="\t", header=FALSE)
m_dat <- masked_dat[,6]

bin_size <- as.numeric(args[2])
# Assign Variables for plotting
plot <- TRUE
band <- 1
main <- "GC Content Percentage"
maxbin.num <- NULL
bin <- bin_size
R <- 1 # this will be changed to 2 when I plot repeats (masking n's)

#Assuming %, this is for GC
legend.max <- 100
#Change colors here
col <- c("black", "red", "orange", "yellow", "darkgreen", "grey4")
col2 <- c("black", "green", "yellow")
col.seg <- NULL
col2.seg <- NULL
width <- 5
legend.len <- 100
legend.cex <- 4
cex.axis <- 1
legend.y.intersp <- 1
legend.x.intersp <- 1
legend.pt.cex <- 3
file.output <- TRUE
cir.chr <- TRUE
cir.chr.h <- 8.5
chr.den.col <- c("black", "red", "yellow", "darkgreen", "grey4")
cir.legend=TRUE
cir.legend.cex <- 8.6
cir.legend.col <- "black"
LOG10 <- TRUE
box <- FALSE
file <- "jpg"    #output type -> "jpg", "png", "tiff"
dpi <- 300
r <- 15.5
H <- 1
cir.band <- 2
Max <- 5

if(file=="jpg") { jpeg(paste("Circular GC and Mask Plot.jpg"), width = 9*dpi,height=7*dpi,res=dpi,quality = 100) }
# if(file=="pdf") { pdf(paste(paste("Circular GC and Mask Plot.pdf"), width = 9,height=7) }
# if(file=="tiff") { tiff(paste("Circular GC and Mask Plot.tiff"), width = 9*dpi,height=7*dpi,res=dpi) }
# if(file=="png") { jpeg(paste("Circular GC and Mask Plot.png")sep=""), width = 9*dpi,height=7*dpi,res=dpi,quality = 100) }
par(xpd=TRUE)


# change from "0.525000" to "52", -- YES IT ROUNDS DOWN
n_dat$V5 <- as.numeric(as.character(sub("0." , "", n_dat$V5)))
n_dat$V5 <- substr(n_dat$V5, 0, 2)
n_dat$V5 <- formatC(as.numeric(n_dat$V5),width=2,format='f',digits=0,flag='0')
m_dat <- as.numeric(m_dat)/bin_size
m_dat <- ceiling(100*m_dat)/100
m_dat <- as.numeric(as.character(sub("0." , "", m_dat)))
m_dat <- substr(m_dat, 0, 2)

# make the curvature for the circular plot
circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
{
  curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
  curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
}

# make as matrix
zzz <- as.matrix(n_dat)
zx <- as.matrix(m_dat)

# check for NA and 0's
zzz <- zzz[zzz[, 1] != 0, ]
max.chr <- max(as.numeric(zzz[, 1]), na.rm=TRUE)
if(is.infinite(max.chr))  max.chr <- 0
zzz.xy.index <- which(!as.numeric(zzz[, 1]) %in% c(0 : max.chr))

#if(length(zzz.xy.index) != 0){
chr.xy <- unique(zzz[zzz.xy.index, 1])
#  for(i in 1:length(chr.xy)){
#    zzz[zzz[, 1] == chr.xy[i], 1] <- max.chr + i
#  }
#}

zzz <- zzz[order(as.numeric(zzz[, 1]), as.numeric(zzz[, 2])), ]
chr <- zzz[, 1]
pos <- as.numeric(zzz[, 2])
gc <- as.numeric(zzz[,3])
mask <- as.numeric(zx[,1])
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
if (band==0)  band=1

if (LOG10) {
  pvalueT[pvalueT <= 0] <- 1
  pvalueT[pvalueT > 1] <- 1
}

for (i in 0:(Nchr-1)) {
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

pos.x <- list()
col.index <- list()
col2.index <- list()

for(i in 1 : length(chr.num)) {
  pos.x[[i]] <- pos[which(chr == chr.num[i])]
  col.index[[i]] <- gc[which(chr == chr.num[i])]
  col2.index[[i]] <- mask[which(chr == chr.num[i])]
}

# I want this 100 so it is percent by definition
col <- colorRampPalette(col)(100)
col2 <- colorRampPalette(col2)(100)

# THIS IS IT
for (i in 1 : length(chr.num)) {
  col.seg <- c(col.seg, col[col.index[[i]]])
}

for (i in 1 : length(chr.num)) {
  col2.seg <- c(col2.seg, col2[col2.index[[i]]])
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
legend2.y <- seq(10,90,5)
legend.y.col <- legend.y
legend2.y.col <- legend2.y

# pick every 5 colors from palette
legend.col <- col[seq(10,90,5)]
legend2.col <- col2[seq(10,90,5)]

if (length(chr.den.col) > 1) {
  cir.density=TRUE
  den.fold <- 20
  #density.list <- list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
  
} else {
  cir.density <- FALSE
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


if (cir.density) {
  segments(
    (1.3 * RR)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (1.3 * RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    col=col.seg, lwd=0.8
  )
  legend(
    75,
    1,
    title="GC %", legend=legend.y, pch=15, pt.cex = 1, col = legend.col,
    cex=0.5, bty="n",
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


if (cir.density) {
  segments(
    (8*cir.band+RI)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2))/TotalN),
    (8*cir.band+RI+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2))/TotalN),
    col=col2.seg, lwd=0.8
  )
  legend(
    86,
    1,
    title="Mask %", legend=legend2.y, pch=15, pt.cex = 1, col = legend2.col,
    cex=0.5, bty="n",
    y.intersp=1,
    x.intersp=1,
    yjust=0.5, xjust=0, xpd=TRUE
  )
}

dev.off()