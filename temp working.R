setwd('~/Documents/QTL')
library(GenomicRanges)
genes <- read.table('gencode.v19.annotation.gtf', skip=5, sep='\t', stringsAsFactors = FALSE)
whichCodingTranscripts <- genes[, 3] == "transcript" & grepl("transcript_type protein_coding", genes[, 9], fixed = TRUE)
proteinTranscripts <- genes[whichCodingTranscripts, ]
strands <- proteinTranscripts[, 7]
allFeaturesTranscripts <- gsub("transcript_id ", '', sapply(strsplit(genes[, 9], "; "), '[', 2))
proteinTranscriptsNames <- allFeaturesTranscripts[whichCodingTranscripts]
whichCDS <- genes[, 3] == "CDS" & allFeaturesTranscripts %in% proteinTranscriptsNames
transcriptsCDS <- genes[whichCDS, ]
transcriptsCDS <- split(GRanges(transcriptsCDS[, 1], IRanges(transcriptsCDS[, 4], transcriptsCDS[, 5]), transcriptsCDS[, 7]),
                        factor(allFeaturesTranscripts[whichCDS], levels = proteinTranscriptsNames))
firstCDS <- mapply(function(CDS, strand) {if(strand == '+') {CDS[1]} else {CDS[length(CDS)]}}, transcriptsCDS, strands)
lastCDS <-  mapply(function(CDS, strand) {if(strand == '+') {CDS[length(CDS)]} else {CDS[1]}}, transcriptsCDS, strands)
whichUTR <- genes[, 3] == "UTR" & allFeaturesTranscripts %in% proteinTranscriptsNames
transcriptsUTR <- genes[whichUTR, ]
transcriptsUTR <- split(GRanges(transcriptsUTR[, 1], IRanges(transcriptsUTR[, 4], transcriptsUTR[, 5]), transcriptsUTR[, 7]),
                        factor(allFeaturesTranscripts[whichUTR], levels = names(firstCDS)))

transcriptsUTR5 <- mapply(function(UTR, CDS, strand)
{        
  if(strand == '+') UTR[UTR < CDS[1]] else UTR[UTR > CDS[length(CDS)]]
}, transcriptsUTR, firstCDS, as.list(strands), SIMPLIFY = FALSE)

transcriptsUTR3 <- mapply(function(UTR, CDS, strand)
{        
  if(strand == '+') UTR[UTR > CDS[length(CDS)]] else UTR[UTR < CDS[1]]
}, transcriptsUTR, firstCDS, as.list(strands), SIMPLIFY = FALSE)
