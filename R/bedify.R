bedify <- function(matched = matched, num = c(1, 3)[1], saveBED = F, outname = outname, cores = 8) {
  tID <- matched[[num]]$transcriptID
  eiID <- matched[[num]]$input_id
  toBed <- list()
  toBed <- parallel::mclapply(1:length(tID), mc.cores = cores, function(i) {
    bed <- gtf[gtf$transcriptID == tID[i] & gtf$type == "exon",] %>% dplyr::select(chr, start, stop, transcriptID, geneID, strand)
    bed$eiID <- rep(eiID[i], length(gtf$geneID[gtf$transcriptID == tID[i] & gtf$type == "exon"]))
    bed$score <- 0
    bed$squish <- paste(bed$transcriptID, "#", bed$eiID, sep = "")

    colnames(bed) <- c('chrom','chromStart', 'chromEnd', 'transcriptID', "geneID", "strand", "eiID", "score", "name")
    bed <- bed %>% dplyr::select(chrom, chromStart, chromEnd, name, score, strand)
    if (unique(bed$strand) == "+") {
      bed$chromStart <- as.integer(bed$chromStart) - 1
    } else {
      bed$chromStart <- as.integer(bed$chromEnd) + 1
    }
    bed
  })
  toBed <- do.call(rbind, toBed)
  if (saveBED == T) {
    write.table(toBed, paste("./", outname, ".bed", sep = ""), sep = '\t', col.names = F, quote = F, row.names = F)
  }
  return(toBed)
}
