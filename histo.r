#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
library("stringi")

new_dat <- read.csv(args[1], sep="\t", header=FALSE)
n_dat <- new_dat[,c(1,2,5)]

# change from "0.525000" to "52", yes it ROUNDS DOWN
n_dat$V5 <- as.character(sub("0." , "", n_dat$V5))
n_dat$V5 <- substr(n_dat$V5, 0, 2)
n_dat$V5 <- stri_pad_right(n_dat$V5, 2, 0)
n_dat$V5 <- formatC(as.numeric(n_dat$V5),width=2,format='f',digits=0,flag='0')
zzz <- as.matrix(n_dat)
gc <- as.numeric(zzz[,3])

ggplot(data=n_dat, aes(gc)) + 
    geom_histogram(breaks = seq(0,100,1), 
                   aes(fill=..count..)) + labs(x="GC%", y="Frequency") +
    scale_fill_gradient("Frequency", low = args[2], high = args[3])

gc_tab <- table(gc)
d <- prop.table(gc_tab)

write.table(d)
