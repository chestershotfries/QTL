library(ggplot2)
library(dplyr)
library(tidyverse)
setwd('~/Documents/QTL')

kozakseq <- read.table('all.selected.kozak.uniq.split.fa')
colnames(kozakseq) <- c("gene","seq")
kozakseq$seq <- toupper(kozakseq$seq)
seq <- data.frame(str_split_fixed(kozakseq$seq, "", max(nchar(kozakseq$seq))))
rownames(pwm) <- c("A","C","G","T")
a <- function(x) {
  g <- as.data.frame(c(length(grep("A",x)),length(grep("C",x)),length(grep("G",x)),length(grep("T",x))))
  return(.g)
}

pwm.t <- as.data.frame(apply(seq,2,a), row.names = c("A","C","G","T"), col.names =)
names(pwm.t) <- NULL

pwm <- pwm.t %>%
  mutate(across(where(is.numeric), ~ ./sum(.)))
