#!/usr/bin/env Rscript

tn <- Sys.getenv("TN")
ntn <- Sys.getenv("NTN")
fdr <- Sys.getenv("FDR")

# teQTL
tnfile= paste(get("tn"), "selected.start.codons.bed", sep=".")
tnbed = read.table(tnfile)
colnames(tnbed) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
# Designate Kozak region as -9 to +6 of start codon
# 9 before start position to 6 after start position if strand is +; 6 before start position to 9 after start position if strand is -
tnbed$chromStart <- ifelse(tnbed$strand == '+', tnbed$chromStart - 9, tnbed$chromStart - 6)
tnbed$chromEnd <- ifelse(tnbed$strand == '+', tnbed$chromEnd + 6, tnbed$chromEnd + 9)
# Remove chr from chromosome number column
tnbed$chrom <- gsub('^chr', '', tnbed$chrom)
write.table(tnbed, paste(get("tn"),"kozak.bed", sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')

# Take 200bp upstream
tnupbed <- tnbed
tndownbed <- tnbed
tnupbed$chromStart <- ifelse(tnupbed$strand == '+', tnupbed$chromStart - 200, tnupbed$chromStart + 18)
tnupbed$chromEnd <- ifelse(tnupbed$strand == '+', tnupbed$chromEnd - 18, tnupbed$chromEnd + 200)
# Take 200bp downstream
tndownbed$chromStart <- ifelse(tndownbed$strand == '+', tndownbed$chromStart + 18, tndownbed$chromStart - 200)
tndownbed$chromEnd <- ifelse(tndownbed$strand == '+', tndownbed$chromEnd + 200, tndownbed$chromEnd - 18)
write.table(tnupbed, paste(get("tn"),"kozak.upstream.bed",sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')
write.table(tndownbed, paste(get("tn"),"kozak.downstream.bed",sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')

# Non-teQTL
ntnfile= paste(get("ntn"), "selected.start.codons.bed", sep=".")
ntnbed = read.table(ntnfile)
colnames(ntnbed) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand')
# Designate Kozak region as -9 to +6 of start codon
# 9 before start position to 6 after start position if strand is +; 6 before start position to 9 after start position if strand is -
ntnbed$chromStart <- ifelse(ntnbed$strand == '+', ntnbed$chromStart - 9, ntnbed$chromStart - 6)
ntnbed$chromEnd <- ifelse(ntnbed$strand == '+', ntnbed$chromEnd + 6, ntnbed$chromEnd + 9)
# Remove chr from chromosome number column
ntnbed$chrom <- gsub('^chr', '', ntnbed$chrom)
write.table(ntnbed, paste(get("ntn"),"kozak.bed", sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')

# Take 200bp upstream
ntnupbed <- ntnbed
ntndownbed <- ntnbed
ntnupbed$chromStart <- ifelse(ntnupbed$strand == '+', ntnupbed$chromStart - 200, ntnupbed$chromStart + 18)
ntnupbed$chromEnd <- ifelse(ntnupbed$strand == '+', ntnupbed$chromEnd - 18, ntnupbed$chromEnd + 200)
# Take 200bp downstream
ntndownbed$chromStart <- ifelse(ntndownbed$strand == '+', ntndownbed$chromStart + 18, ntndownbed$chromStart - 200)
ntndownbed$chromEnd <- ifelse(ntndownbed$strand == '+', ntndownbed$chromEnd + 200, ntndownbed$chromEnd - 18)
write.table(ntnupbed, paste(get("ntn"),"kozak.upstream.bed",sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')
write.table(ntndownbed, paste(get("ntn"),"kozak.downstream.bed",sep="."), quote = FALSE, row.names = FALSE, col.names =  FALSE, sep ='\t')