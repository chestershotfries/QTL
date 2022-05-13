setwd('~/Documents/QTL')
qtl = read.table('cis.teQTL.mapping.permutation.pass.results.table')
qtlbed = read.table('qtl.bed')
  qtlbed <- qtlbed[,-1]
  qtlbed$V6 <- qtlbed$V2
  qtlbed <- qtlbed[,-1]
  qtlbed$V7 <- qtlbed$V4
  qtlbed$V4 <- qtlbed$V5
  qtlbed$V8 <- qtlbed$V7
  qtlbed$V7 <- 0
  colnames(qtlbed) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
qtlbed$chromStart <- qtlbed$chromStart-10
qtlbed$chromEnd <- qtlbed$chromEnd+10
qtl.kozak10.bed <- qtlbed
write.table(qtl.kozak10.bed, 'kozak10.bed', quote = FALSE, row.names = FALSE, col.names =  FALSE, sep = '\t')
qtlbed$chromStart <- ifelse(qtlbed$strand == '+', qtlbed$chromStart - 2, qtlbed$chromStart - 6)
qtlbed$chromEnd <- ifelse(qtlbed$strand == '+', qtlbed$chromEnd + 6, qtlbed$chromEnd + 2)
write.table(qtlbed, 'kozak.bed', quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')

teqtlbed = read.table('teqtl.selected.start.codons.bed')
colnames(teqtlbed) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
teqtlbed$chromStart <- ifelse(teqtlbed$strand == '+', teqtlbed$chromStart - 6, teqtlbed$chromStart - 4)
teqtlbed$chromEnd <- ifelse(teqtlbed$strand == '+', teqtlbed$chromEnd + 4, teqtlbed$chromEnd + 6)
teqtlbed$chrom <- gsub('^chr', '', teqtlbed$chrom)
write.table(teqtlbed, 'teqtl.kozak.bed', quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')
