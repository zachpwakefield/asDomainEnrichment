proteinExtract_pipe <- function(files_dir, background = T, mOverlap = .5, saveOutput = F, inCores = 8, nC = 0, nE = 0, exon_type = "AFE", location = system.file(package="domainEnrichment"), output_location) {

  if (background == T) {
    files <- paste(files_dir, list.files(files_dir)[grep('[.]exon', list.files(files_dir))], sep = "")
    cat(files)
    first_exons <- unique(unlist(lapply(files, function(x) {
      in_file <- read.delim(x)
      in_file <- in_file[in_file$ID == 'first',]
      paste(in_file$gene, ';', in_file$exon, ';',  in_file$strand, sep = "")})))
    redExon <- data.frame(geneR = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 1)), split = '[.]'), "[[", 1)),
                          chr = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 1)),
                          start = unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 1)), split = ":"), "[[", 2)),
                          stop = unlist(lapply(strsplit(unlist(lapply(strsplit(first_exons, split = ";"), "[[", 2)), split = '-'), "[[", 2))
    )
  } else {
    df <- read.delim(files_dir, sep = " ")
    df.l <- lfc(df, numCont = nC, numExp = nE, exon_type = exon_type, cores = inCores)
    lfcPlot <- make_lfcPlot(df.l, location = location)

    redExon <- data.frame(geneR = unlist(lapply(strsplit(df.l$gene, split = "[.]"), "[[", 1)),
                          chr = sapply(strsplit(df.l$exon, split = ":"), "[[", 1),
                          start = sapply(strsplit(sapply(strsplit(df.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 1),
                          stop = sapply(strsplit(sapply(strsplit(df.l$exon, split = ":"), "[[", 2), split = "[-]"), "[[", 2)
    )

  }


  colnames(redExon) <- c("geneR", "chr", "start", "stop")
  redExon$start <- as.numeric(redExon$start)
  redExon$stop <- as.numeric(redExon$stop)
  print("exon loaded...")


  matched <- getTranscript(gtf = gtf, redExon = redExon[1:200,], ex_type = exon_type, minOverlap = mOverlap, swaps = !(background), cores = inCores)
  print("exons matched, bed-ifying...")
  bed <- bedify(matched, saveBED=F, outname = outname, cores = inCores)


  trans <- unlist(lapply(strsplit(unique(bed$name), "#"), "[[", 1))
  possT <- unlist(lapply(strsplit(bed$name, "#"), "[[", 1))


  ## Find annotated proteins for transcripts if possible
  protCode <- unlist(mclapply(trans, mc.cores = 8, function(x) {
    rc <- c_trans[which(c_trans == x)+1]
    if (length(rc) > 0) {
      rc[1]
    } else {"none"}
  }))

  proBed <- data.frame(id = unique(bed$name), strand = unlist(lapply(unique(bed$name), function(x) unique(bed$strand[bed$name == x][1])[1])), prot = protCode) %>% tidyr::separate(id, c("transcript", "id"), "#") %>% tidyr::separate("id", c("gene", "chr"), ";") %>% tidyr::separate('chr', c('chr', 'coords'), ':') %>% tidyr::separate('coords', c('start', 'stop'), '-')

  proFast <- c()
  if (length(proBed[,1]) %% 2 == 0) {
    subVal <- 1
  } else {subVal <- 0}
  for (i in seq(1, (length(proBed[,1])-subVal), by = 2)) {
    proFast <- c(proFast, paste(">", proBed$transcript[i], "#", proBed$gene[i], ";", proBed$chr[i], ":", proBed$start[i], "-", proBed$stop[i], ";", proBed$strand[i], sep = ""),
                 proBed$prot[i], paste(">", proBed$transcript[i+1], "#", proBed$gene[i+1], ";", proBed$chr[i+1], ":", proBed$start[i+1], "-", proBed$stop[i+1], ";", proBed$strand[i+1], sep = ""),
                 proBed$prot[i+1])
  }

  if (background == T) {
    if (saveOutput == T) {
      write_csv(proBed, paste0(output_location, "bgoutBed.csv"))
      write_lines(proFast, paste0(output_location, "bgoutFast.fa"))
      write_csv(matched$out_matched,  paste0(output_location, "bgmatched.csv"))
      write_csv(bed,  paste0(output_location, "bgbed.csv"))

    }
    return(list(matched = matched,
                bed = bed,
                proBed = proBed,
                proFast = proFast))
  } else {
    protAlign <- list()
    protC <- c()
    pMatch <- c()
    alignType <- c()
    cate <- c()
    for (i in seq(from=1,to=(length(protCode)-1), by=2)) {
      if (protCode[i] == "none" | protCode[i+1] == "none") {
        if (protCode[i] == "none" & protCode[i+1] != "none") {
          protC <- c(protC, "nonPC", "PC")
          protAlign[[i]] <- "none"
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else if (protCode[i] != "none" & protCode[i+1] == "none") {
          protC <- c(protC, "PC", "nonPC")
          protAlign[[i]] <- "onePC"
          alignType <- c(alignType, "onePC")
        } else {
          protC <- c(protC, "nonPC", "nonPC")
          protAlign[[i]] <- "none"
          alignType <- c(alignType, "noPC")
        }
        pMatch <- c(pMatch, 0)
      } else if (protCode[i] == protCode[i+1]) {
        protC <- c(protC, "Same", "Same")
        protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])))
        pMatch <- c(pMatch, 1.04)
        alignType <- c(alignType, "Match")
      } else {
        protC <- c(protC, "Different", "Different")
        protAlign[[i]] <- msa::msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE)

        minPc <- min(nchar(protCode[i]), nchar(protCode[i+1]))
        pMatch <- c(pMatch, table(unlist(lapply(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = ""), function(x) x == "?")))[1]/min(nchar(protCode[i]), nchar(protCode[i+1])))
        if (nchar(paste(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]][nchar(strsplit(msa::msaConsensusSequence(protAlign[[i]]), split = "\\?|\\.|!")[[1]]) > (.1*minPc)], collapse = "")) > .2*minPc)  {
          alignType <- c(alignType, "PartialMatch")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_pm_prettyAlignment.pdf", sep = ""), output = "pdf")
        } else {
          alignType <- c(alignType, "FrameShift")
          # msaPrettyPrint(msa(Biostrings::AAStringSet(c(protCode[i], protCode[i+1])), verbose = FALSE), askForOverwrite=FALSE
          # , file = paste(out_dir, "prettyAlignments/", proBed$transcript[i], "_", proBed$transcript[i+1], "_fs_prettyAlignment.pdf", sep = ""), output = "pdf")
        }
      }
    }

    print(table(protC))
    proBed$match <- protC
    proBed$prop <- rep(pMatch, each = 2)

    # Filled Density Plot
    (gdf <- ggplot2::ggplot(data.frame(dens = as.numeric(pMatch), type = alignType), ggplot2::aes(x = dens, fill = type)) +
        ggplot2::geom_histogram(ggplot2::aes(y=after_stat(count)/sum(after_stat(count))), colour = 1,
                       bins = 20) + ggplot2::geom_density(ggplot2::aes(y=.0005*after_stat(count)), color = 'black', fill = "coral2", bw = .1, alpha = .3) +
        ggplot2::scale_fill_manual(values=c('noPC' = "azure4", 'Match' = "#E69F00", 'onePC' = "#56B4E9", 'FrameShift' = "pink", 'PartialMatch' = "deeppink4")) +
        ggplot2::theme_classic() + ggplot2::xlab("Alignment Score") + ggplot2::ylab("Fraction"))


    proBed$matchType <- rep(alignType, each = 2)

    if (saveOutput == T) {
      write_csv(proBed, paste0(output_location, "fgoutBed.csv"))
      write_lines(proFast, paste0(output_location, "fgoutFast.fa"))
      write_csv(matched$out_matched,  paste0(output_location, "fgmatched.csv"))
      write_csv(bed,  paste0(output_location, "fgbed.csv"))
      write_csv(df.l,  paste0(output_location, "fglfc.csv"))

      pdf(file = paste0(output_location, "alignPlot.pdf"))
      gdf
      dev.off()
      pdf(file = paste0(output_location, "volcano.pdf"))
      print(deExons)
      dev.off()
    }
    return(list(matched = matched,
                bed = bed,
                proBed = proBed,
                proFast = proFast,
                df.l = df.l,
                gdf = gdf,
                deExons = lfcPlot))
  }

}
