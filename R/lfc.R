lfc <- function(de_df, numCont, numExp, exon_type, cores = 8) {
  de_df <- de_df[de_df$type == exon_type,]
  samps <- colnames(de_df)

  lfc <- list()
  col <- list()


  lfc <- mclapply(1:length(de_df$gene), mc.cores = cores, function(i) {
    outlier <- which(unlist(lapply(samps, function(x) grepl(x, de_df$outlier[i]))))+10

    cont <- c(11:(11-1+numCont))
    cont <- cont[!(cont %in% outlier)]
    exp <- c((numCont+10+1):(11-1+numCont+numExp))
    exp <- exp[!(exp %in% outlier)]



    if (length(exp) == 0) {
      lfc_t <- log2((.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    }
    if (length(cont) == 0) {
      lfc_t <- log2((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/.01)
    }
    if (length(exp) == 0 & length(cont) == 0) {
      lfc_t <- 0
    }
    if (length(exp) != 0 & length(cont) != 0) {
      lfc_t <- log((as.numeric(rowMeans(de_df[i,exp, drop = FALSE]))+.01)/as.numeric(rowMeans(de_df[i,cont, drop = FALSE])+.01))
    }


    lfc_t
  })

  de_df$lfc <- unlist(lfc)
  de_df <- de_df[de_df$p_value >= 0,]
  col <- list()
  for (i in 1:length(de_df$gene)) {
    col[[i]] <- "#A7A9AC"#natparks.pals("Acadia", 15)[7]
      if (de_df$lfc[i] <= -1.0 & de_df$p_value[i] < .01) {
        col[[i]] <- "#FE9234" #natparks.pals("Acadia", 15)[15]
      }
    if (de_df$lfc[i] >= 1.0 & de_df$p_value[i] < .01) {
      col[[i]] <- "#00A79D"# natparks.pals("Acadia", 15)[1]
    }

  }
  de_df$col <- unlist(col)

  p.scDE <- de_df
  p.scDE$numOutliers <- str_count(p.scDE$outlier, "[.]")
  return (de = p.scDE)
}
