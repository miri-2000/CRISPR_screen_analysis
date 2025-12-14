# Helper file to create correlation and distribution plots

plotDistributionsGeneTypes <- function(
  distr,
  from = NULL,
  to = NULL,
  plotTitle = "",
  xlab = "log2FoldChange",
  yHigh = NULL,
  xBreaks = NULL
) {
  if (is.null(from)) {
    from <- distr$from
  }
  if (is.null(to)) {
    to <- distr$to
  }

  #intersection distribution with the positive contrl
  pd <- distr$p
  nd <- distr$n
  if (length(pd) != 0 & length(nd) != 0) {
    vLine <- distrIntersection(pd, nd)$x[[1]]
  } else {
    vLine <- NULL
  }

  #Annotation
  a <- data.frame(
    breaks = c("p", "n", "x", "a", "o", "o2"),
    labels = c(
      "pos. control",
      "neg. control",
      "experimental",
      "all",
      "other",
      "other2"
    ),
    colors = c("red", "blue", "darkgreen", "darkgrey", "orange", "pink"),
    stringsAsFactors = F
  )
  #subset annotation for the distributions present
  a <- a[a$breaks %in% names(distr), ]

  #for colorValues create named vector
  colorValues <- a$colors
  names(colorValues) <- a$breaks

  p <- plotDistributions(
    distr,
    from = from,
    to = to,
    plotTitle = plotTitle,
    vLine = vLine,
    legendBreaks = a$breaks,
    legendLabels = a$labels,
    colorValues = colorValues,
    xlab = xlab,
    yHigh = yHigh,
    xBreaks = xBreaks
  )

  return(p)
}

distrIntersection <- function(pd, nd) {
  # create True / False comparing nd$y < pd$y
  # create a diff vector in which two consecutive T/F are compared, to see where they change. In the comparison t
  # although the difference vector is one smaller, based on a big enough span (currently using 2048) it is still a good approximation
  x <- pd$x[as.logical(abs(diff(pd$y < nd$y)))]
  y <- pd$y[as.logical(abs(diff(pd$y < nd$y)))]

  # calculate median of negative control
  # http://stats.stackexchange.com/questions/32093/how-to-draw-mean-median-and-mode-lines-in-r-that-end-at-density
  n <- length(nd$y)
  y.cs <- cumsum(nd$y)
  x.med <- nd$x[i.med <- length(y.cs[2 * y.cs <= y.cs[n]])]
  #y.med <- nd$y[i.med]

  ## select the first intersection starting left from the median of the negative control.
  if (length(x[x < x.med]) == 0) {
    o <- NULL
    warning("length(x[x< x.med]==0, no intersection could be determined")
  } else {
    sel <- x == max(x[x < x.med])
    x <- x[sel]
    y <- y[sel]
    o <- list(x = x, y = y)
  }
  return(o)
}

plotDistributions <- function(
  distr,
  from,
  to,
  plotTitle,
  xBreaks = NULL,
  legendBreaks = NULL,
  legendLabels = NULL,
  vLine = NULL,
  colorValues = NULL,
  xlab = "log2FoldChange",
  yHigh = NULL,
  legendLimits = NULL,
  medians = NULL,
  floorCeiling = T,
  yLimits = NULL,
  yBreaks = NULL
) {
  ## Create from each of the distributions and x and y data from, and merge them into a data.frame
  ## with one column x and for every distr a y column
  td <- NULL

  #thought before that values are matched with breaks but it is limits.
  if (is.null(legendLimits) & !is.null(legendBreaks)) {
    legendLimits <- legendBreaks
  }

  for (i in c(1:(length(distr)))) {
    d <- distr[[i]]
    name <- names(distr)[i]
    if (!(name %in% c("from", "to"))) {
      #filter out from and to
      xy <- data.frame(xAxis = d$x, y = d$y)
      colnames(xy)[2] <- name

      #total
      if (is.null(td)) {
        td <- xy
      } else {
        td <- merge(td, xy, by = "xAxis")
      }
    }
  }

  library(reshape2)
  library(ggplot2)

  #all Distr
  td <- melt(td, id.vars = "xAxis")
  colnames(td)[2] <- "type"

  if (floorCeiling) {
    xLow <- floor(from)
    xHigh <- ceiling(to)
  } else {
    xLow <- from
    xHigh <- to
  }

  if (is.null(xBreaks)) {
    xBreaks <- c(xLow:xHigh)
  }

  p <- ggplot(td, aes(x = xAxis, y = value, color = type)) +
    # if td contains xAxis values > 3. Setting the limit to 3, will generate warning: 'Removed 1436 rows containing missing values (geom_path). '
    scale_x_continuous(limits = c(xLow, xHigh), breaks = xBreaks) +
    geom_line(linewidth = 1)

  if (!is.null(medians)) {
    p <- p +
      geom_vline(
        data = medians,
        aes(xintercept = value, color = type),
        linetype = "dashed"
      )
  }

  if (!is.null(yLimits)) {
    if (is.null(yBreaks)) {
      yBreaks <- waiver()
    }
    p <- p + scale_y_continuous(limits = yLimits, breaks = yBreaks)
  } else {
    if (!is.null(yHigh)) {
      p <- p +
        scale_y_continuous(
          limits = c(0, yHigh),
          breaks = c(1:(yHigh * 10) / 10)
        )
    }
  }

  if (!is.null(vLine)) {
    p <- p + geom_vline(xintercept = vLine, linetype = "dotted")
  }

  # draw horizontal and vertical axis
  p <- p + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  p <- p +
    xlab(xlab) +
    ylab("Probability") +
    ggtitle(plotTitle)

  if (!is.null(colorValues)) {
    p <- p +
      scale_colour_manual(
        name = "",
        limits = legendLimits,
        values = colorValues,
        breaks = legendBreaks,
        labels = legendLabels,
        guide = guide_legend(
          title = "type",
          label.theme = element_text(size = 18, color = "black", angle = 0)
        )
      )
  }

  p <- p +
    theme(
      axis.text = element_text(size = 14, face = "bold"), #for ticks
      axis.title = element_text(size = 18, face = "bold", vjust = 0.3), #for the labels
      plot.title = element_text(face = "bold", size = 20, vjust = 1), #
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 18),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
      panel.background = element_rect(fill = "white", colour = "grey"),
      axis.title.y = element_text(vjust = 1)
    ) #legend.title

  return(p)
}


getDistributions <- function(
  x,
  types = c("p", "n", "x", "a", "o", "o2"),
  limit = 10,
  other = NULL
) {
  #In order to be able to set a common range for the distribution. Truncate -Inf and values higher than -10 to -10

  if (!is.null(limit)) {
    x[!is.na(x$value) & x$value < -limit, "value"] <- -limit
    x[!is.na(x$value) & x$value > limit, "value"] <- limit
  }

  pd <- NULL
  nd <- NULL
  xd <- NULL
  ad <- NULL
  od <- NULL
  o2d <- NULL

  from <- min(x[, "value"], na.rm = T)
  to <- max(x[, "value"], na.rm = T)

  if (!is.null(types)) {
    if ("p" %in% types) {
      pv <- x[x$type == "p", "value"]
      if (length(pv) > 0) {
        pd <- density(pv, from = from, to = to, n = 2048, na.rm = T)
      }
    }
    if ("n" %in% types) {
      nv <- x[x$type == "n", "value"]
      if (length(nv) > 0) {
        nd <- density(nv, from = from, to = to, n = 2048, na.rm = T)
      }
    }
    if ("x" %in% types) {
      xv <- x[x$type == "x", "value"]
      if (length(xv) > 0) {
        xd <- density(xv, from = from, to = to, n = 2048, na.rm = T)
      }
    }
    if ("o" %in% types) {
      ov <- x[x$type == "o", "value"]
      if (length(ov) > 0) {
        od <- density(ov, from = from, to = to, n = 2048, na.rm = T)
      }
    }
    if ("o2" %in% types) {
      o2v <- x[x$type == "o2", "value"]
      if (length(o2v) > 0) {
        o2d <- density(o2v, from = from, to = to, n = 2048, na.rm = T)
      }
    }
    if ("a" %in% types) {
      av <- x[, "value"]
      if (length(av) > 0) {
        ad <- density(av, from = from, to = to, n = 2048, na.rm = T)
      }
    }
  }

  if (!is.null(other)) {
    ov <- x[other, "value"]
    if (length(ov) > 0) {
      od <- density(ov, from = from, to = to, n = 2048, na.rm = T)
    }
  }

  res <- list()
  if (!is.null(pd)) {
    res$p <- pd
  }
  if (!is.null(nd)) {
    res$n <- nd
  }
  if (!is.null(xd)) {
    res$x <- xd
  }
  if (!is.null(ad)) {
    res$a <- ad
  }
  if (!is.null(od)) {
    res$o <- od
  }
  if (!is.null(o2d)) {
    res$o2 <- o2d
  }

  res$from = from
  res$to = to
  return(res)
}


fMeasure <- function(p, r, b = 1) {
  fMeasure <- (1 + b^2) * (p * r) / ((b^2 * p) + r)
}

fMeasureScreenHps <- function(x) {
  #get distributions

  distrs <- getDistributions(x)

  # calculate intersection
  int <- distrIntersection(distrs$p, distrs$n)

  #calculate 'recall'
  #=fraction of positive controls selected (=left of intersection) compare tot total
  ps <- x[x$type == "p", "value"]
  psHits <- ps[ps <= int$x]
  recall <- length(psHits) / length(ps)

  #calculate 'precision' , fraction of positives in the "hitlist"
  ns <- x[x$type == "n", "value"]
  nsHits <- ns[ns <= int$x]
  precision <- length(psHits) / (length(psHits) + length(nsHits))
  fm <- fMeasure(p = precision, r = recall)

  return(fm)
}

distrPlot <- function(
  x,
  annotation_columns = NULL,
  exp = NULL,
  ref = NULL,
  to = 3,
  title,
  paired = F,
  other = NULL,
  yAxisMax = NULL,
  drawMedian = F,
  drawMean = F,
  types = c("p", "n", "x", "a", "o", "o2"),
  from = -10
) {
  if (!is.null(exp) | !is.null(ref)) {
    #log2Fcs needs still be calculated
    expData <- x[, grep(exp, colnames(x)), drop = F]
    refData <- x[, grep(ref, colnames(x)), drop = F]

    if (paired) {
      expReps <- sapply(colnames(expData), FUN = function(x) {
        substring(x, nchar(x))
      })
      names(expReps) <- NULL #otherwise checkEquals will give unequal
      refReps <- sapply(colnames(refData), FUN = function(x) {
        substring(x, nchar(x))
      })
      names(refReps) <- NULL
      ## the calcu
      checkEquals(
        expReps[order(expReps)],
        refReps[order(refReps)],
        msg = ". For paired calculation, replicates for two conditions are not identical."
      )
      fAll <- NULL
      for (rep in expReps) {
        e <- expData[, grep(paste(rep, "$", sep = ""), colnames(expData))]
        r <- refData[, grep(paste(rep, "$", sep = ""), colnames(refData))]
        f <- log2(e / r)
        if (is.null(fAll)) {
          fAll <- data.frame(f = f, stringsAsFactors = F)
        } else {
          fAll <- cbind(fAll, f = f, stringsAsFactors = F)
        }
        colnames(fAll)[ncol(fAll)] <- paste("log2Fc.r", rep, sep = "")
      }
      log2FoldChange <- apply(fAll, MARGIN = 1, FUN = mean) #is mean over logValues= geometricMean
    } else {
      expM <- apply(expData, MARGIN = 1, FUN = mean)
      refM <- apply(refData, MARGIN = 1, FUN = mean)
      log2FoldChange <- log2(expM / refM)
    }
    #log2FoldChange <- log2(expData/refData)

    x <- cbind(
      x[, c(1:annotation_columns)],
      log2FoldChange = log2FoldChange,
      stringsAsFactors = F
    )

    #calculate fMeasure
    ## in case geneSymbol is written with capital
    colnames(x) <- gsub("Gene", "gene", colnames(x))
    colnames(x) <- gsub("GRNA", "gRNA", colnames(x))

    # Select all columns except "sequence"
    x = x[, !tolower(colnames(x)) %in% "sequence"]
    colnames(x)[ncol(x)] <- "value"
  }

  #In order to be able to set a common range for the distribution. Truncate -Inf and values higher than -10 to -10
  if (length(x[!is.na(x$value) & x$value < -10, "value"]) > 0) {
    x[!is.na(x$value) & x$value < -10, "value"] <- -10
  }

  if (length(x[!is.na(x$value) & x$value > 10, "value"]) > 0) {
    x[!is.na(x$value) & x$value > 10, "value"] <- 10
  }

  #remove id column
  ##TODO replace with fixed calculation in stead of samples
  res = tryCatch(
    {
      #fMeasureScreen(x[,-1],annotation_columns=2)
      x2 <- x[, c("type", "value")]
      x2 <- x2[!is.na(x2$value), ]
      fMeasureScreenHps(x2)
    },
    error = function(e) {
      print(paste(e, "Could not calculate fMeasure"))
      NULL
    }
  )

  # in case tag present
  if (colnames(x)[1] == "gRNA") {
    x2 <- x[, -1]
  } else {
    x2 <- x
  }

  distr <- getDistributions(x2, other = other, types = types)
  pd <- distr[names(distr) == "p"]
  nd <- distr[names(distr) == "n"]

  if (length(pd) > 0 && length(nd) > 0) {
    ## take the first intersection calculated from 0
    intersect <- distrIntersection(pd[[1]], nd[[1]])$x[[1]]
    if (!is.null(res) && !is.null(intersect)) {
      # could not be determined
      intersect <- signif(intersect, 2)
      title <- paste(
        title,
        " (fM.:",
        signif(res, 2),
        ", int.:",
        intersect,
        ")",
        sep = ""
      )
    } else {
      intersect <- NULL
    }
  }

  #remove geneSymbol
  print(paste("Number of rows with log2Fc=NA:", nrow(x[is.na(x$value), ])))
  #distr <- getDistributions(d)
  plot <- plotDistributionsGeneTypes(
    distr,
    plotTitle = title,
    from = from,
    to = to
  )

  if (!is.null(yAxisMax)) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, yAxisMax),
        breaks = c(0:(10 * yAxisMax)) / 10
      )
  }

  ##!!  CODE FOR DRAWING MEDIAN IS ALREADY PRESENT IN plotDistributionsGeneType used that. Built in option there to do the same for means
  if (drawMedian || drawMean) {
    x2p <- x2[x2$type == "p", "value"]
    x2n <- x2[x2$type == "n", "value"]
    x2x <- x2[x2$type == "x", "value"]
    x2o <- x2[x2$type == "o", "value"]

    if (drawMedian) {
      posMedian <- median(x2p, na.rm = T)
      negMedian <- median(x2n, na.rm = T)
      xMedian <- median(x2x, na.rm = T)
      oMedian <- median(x2o, na.rm = T)
      x3 <- data.frame(
        type = c("p", "n", "x", "o"),
        value = c(posMedian, negMedian, xMedian, oMedian)
      )
      x3 <- x3[x3$type %in% types, ]
      plot <- plot +
        geom_vline(
          data = x3,
          aes(xintercept = value, color = type),
          linetype = "longdash"
        )
    }

    if (drawMean) {
      posMean <- mean(x2p, na.rm = T)
      negMean <- mean(x2n, na.rm = T)
      xMean <- mean(x2x, na.rm = T)
      oMean <- mean(x2o, na.rm = T)
      x3 <- data.frame(
        type = c("p", "n", "x", "o"),
        value = c(posMean, negMean, xMean, oMean)
      )
      x3 <- x3[3$type %in% types, ]
      plot <- plot +
        geom_vline(
          data = x3,
          aes(xintercept = value, color = type),
          linetype = "dashed"
        )
    }
  }

  res <- list(plot = plot, x = x, fm = res, intersect = intersect)
  if (drawMedian) {
    res$posMedian <- posMedian
    res$negMedian <- negMedian
  }
  if (drawMean) {
    res$posMean <- posMean
    res$negMean <- negMean
  }

  return(res)
}

controlDistrPlot <- function(
  x,
  exp,
  ref,
  title,
  paired,
  other,
  from,
  to,
  types,
  drawMedian,
  yAxisMax,
  annotation_columns = NULL
) {
  res <- distrPlot(
    x = x,
    exp = exp,
    ref = ref,
    title = title,
    paired = paired,
    other = other,
    from = from,
    to = to,
    types = types,
    drawMedian = drawMedian,
    annotation_columns = annotation_columns
  )
  plot <- res$plot

  #set to the same y range
  # next lines of codes are within a loop therefore not
  if (!is.null(yAxisMax)) {
    plot <- plot +
      scale_y_continuous(
        limits = c(0, yAxisMax),
        breaks = c(0:(10 * yAxisMax)) / 10
      )
  }

  res <- list(plot = plot, x = res$x)
  return(res)
}


controlDistrs <- function(
  x,
  exps,
  refs,
  paired,
  other = NULL,
  yAxisMax = NULL,
  nCol,
  nRow,
  oFileName,
  widthSizeFactor = 7,
  from = -10,
  to = 3,
  types = c("p", "n", "x", "a", "o", "o2"),
  drawMedian = F,
  writeToFile = T,
  annotation_columns = NULL,
  plotformat = "pdf",
  titleLineBreak = F,
  inclRep = T
) {
  plots <- list()
  pIdx <- 0

  xAll <- NULL # table with log2FoldChange values
  # loop over each treatment condition
  for (i in c(1:length(exps))) {
    print(exps[i])
    print(refs[i])

    #per replicate
    reg <- paste("^", exps[i], sep = "")
    if (inclRep == T) {
      reg <- paste(reg, "_", sep = "")
    }
    expData <- x[, grep(reg, colnames(x)), drop = F]
    expReps <- sapply(colnames(expData), FUN = function(x) {
      substring(x, nchar(x))
    })
    if (paired) {
      for (rep in expReps) {
        if (length(expReps) > 9) {
          stop("Code handles up to 9 replicates")
        }
        exp <- paste("^", exps[i], "_.*", rep, "$", sep = "")
        ref <- paste("^", refs[i], "_.*", rep, "$", sep = "")
        title <- paste(exps[i], " vs ", sep = "")
        if (titleLineBreak) {
          title <- paste(title, "\n", sep = "")
        }
        title <- paste(title, refs[i], " ", "for r", rep, sep = "")
        res <- controlDistrPlot(
          x = x,
          exp = exp,
          ref = ref,
          title = title,
          paired = paired,
          other = other,
          from = from,
          to = to,
          types = types,
          drawMedian = drawMedian,
          yAxisMax = yAxisMax,
          annotation_columns = annotation_columns
        )
        plot <- res$plot
        pIdx <- pIdx + 1
        plots[[pIdx]] <- plot
        if (is.null(xAll)) {
          xAll <- res$x1
        } else {
          xAll <- merge(xAll, res$x1, by = c("tag", "geneSymbol", "type"))
        }
        colnames(xAll)[ncol(xAll)] <- paste(exps[i], rep, sep = "_")
      }
    } else {
      exp <- paste("^", exps[i], sep = "")
      ref <- paste("^", refs[i], sep = "")
      title <- paste(exps[i], " vs ", sep = "")
      if (titleLineBreak) {
        title <- paste(title, "\n", sep = "")
      }
      title <- paste(title, refs[i], sep = "")
      res <- controlDistrPlot(
        x = x,
        exp = exp,
        ref = ref,
        title = title,
        paired = paired,
        other = other,
        from = from,
        to = to,
        types = types,
        drawMedian = drawMedian,
        yAxisMax = yAxisMax,
        annotation_columns = annotation_columns
      )
      plot <- res$plot
      pIdx <- pIdx + 1
      plots[[pIdx]] <- plot

      if (is.null(xAll)) {
        xAll <- res$x
      } else {
        cns <- "tag"
        if (length(grep("geneSymbol", colnames(res$x1))) > 0) {
          cns <- c(cns, "geneSymbol")
        }
        cns <- c(cns, "type")
        xAll <- merge(xAll, res$x1, by = cns)
      }

      colnames(xAll)[ncol(xAll)] <- exps[i]
    }
  }

  if (writeToFile) {
    if (is.null(nRow)) {
      if (paired) {
        nRow <- length(exps)
      } else {
        nRow <- 1
      }
      nCol <- length(expReps)
    }

    if (plotformat == "pdf") {
      width <- nCol * widthSizeFactor
      height <- nRow * 5
    } else if (plotformat == "jpg") {
      width <- nCol * widthSizeFactor * 80
      height <- nRow * widthSizeFactor * 50
    }

    oFileName1 <- paste(oFileName, plotformat, sep = "..")

    # plot <- plot + theme(plot.title = element_text(size = 10))
    plot <- ggarrange(
      plotlist = plots,
      ncol = nCol,
      nrow = nRow,
      common.legend = TRUE,
      legend = "bottom"
    ) %>%
      ggexport(filename = oFileName1, width = width, height = height)
  }
  o <- list()
  o[["xAll"]] <- xAll
  o[["plots"]] <- plots

  return(o)
}

showSampleHeatmap <- function(dds, inclRep = T, margin = c(13, 13)) {
  require(gplots)
  library("RColorBrewer")
  rld <- rlogTransformation(dds, blind = TRUE)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  if (inclRep) {
    names <- with(colData(dds), paste(Conditions, rep, sep = "_"))
  } else {
    names <- colData(dds)$Conditions
  }
  rownames(mat) <- colnames(mat) <- names
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, trace = "none", col = rev(hmcol), margin = margin)
}

nrPlotsPerCond <- function(nrReps) {
  #(full matrix - diagonal) / 2 (to get unique combinations)
  y <- ((nrReps)^2 - nrReps) / 2
  return(y)
}

getCorrelationPlot <- function(
  df,
  point.size = 2,
  slope.1 = T,
  typeColoring = F,
  lim = NULL,
  breaks = NULL,
  doPlot = T,
  inTitle = "rs",
  regressionLineColor = "grey",
  units = NULL,
  label = F
) {
  rm <- 0
  # if (label)
  #   rm <- rm + 1

  if (typeColoring) {
    rm <- rm + 1
  }

  # change the column names to x and y
  if (rm == 1) {
    xyNames <- colnames(df[, -c(1:rm)])
    colnames(df)[-c(1:rm)] <- c("X", "Y")
  } else {
    xyNames <- colnames(df)
    colnames(df) <- c("X", "Y")
  }

  # remove all cell values that are either NaN or infinite
  df1 <- df[!is.na(df$X) & is.finite(df$X) & !is.na(df$Y) & is.finite(df$Y), ]
  # fit a linear model using x and y
  df.lm <- lm(Y ~ X, df1)
  # store the adjusted R^2: 1 – [(1-R^2)*(n-1)/(n-k-1)]
  adj.r.squared <- summary(df.lm)$adj.r.squared
  # calculate the pearson correlation

  pearson.corr = cor(df1$X, df1$Y, use = "pairwise.complete.obs")
  if (doPlot) {
    if (inTitle == "rs") {
      title <- paste("R^2=", signif(adj.r.squared, digits = 3))
    }
    if (inTitle == "pc") {
      title <- paste("Pearson corr.=", signif(pearson.corr, digits = 3))
    }

    if (typeColoring) {
      p <- ggplot(df, aes(x = X, y = Y, color = type))
    } else {
      p <- ggplot(df, aes(x = X, y = Y))
    }

    # if (label) {
    #   labelCol <- colnames(df)[1]
    #   p <- p  + geom_text_repel(aes_string(label=labelCol),color="red")
    # }

    #all black

    xLab <- xyNames[1]
    yLab <- xyNames[2]
    # add the units if they were given
    if (!is.null(units)) {
      xLab <- paste(xLab, " (", units, ")", sep = "")
      yLab <- paste(yLab, " (", units, ")", sep = "")
    }

    p <- p +
      geom_point(size = point.size) +
      geom_abline(
        intercept = df.lm$coeff[1],
        slope = df.lm$coeff[2],
        color = regressionLineColor
      ) +
      xlab(xLab) +
      ylab(yLab) +
      theme_bw() +
      ggtitle(title)

    if (typeColoring) {
      # redo the controls to have color
      colorValues <- c(
        "p" = "red",
        "n" = "blue",
        "x" = "black",
        "o" = "orange",
        "o2" = "pink"
      )
      # store the colors for each type that is present in the dataframe
      colorValues <- colorValues[names(colorValues) %in% unique(df$type)]

      # create a new datframe that stores every type except "x" and "e"
      df1 <- df[!df$type %in% c("x", "e"), ]
      p <- p +
        geom_point(data = df1, mapping = aes(x = X, y = Y, color = type)) +
        scale_colour_manual(
          name = "",
          values = colorValues,
          breaks = names(colorValues),
          labels = names(colorValues), #label o only shown if present
          guide = guide_legend(
            title = "type",
            label.theme = element_text(size = 18, color = "black", angle = 0)
          )
        )
    }

    # create a line with a slope of 1
    if (slope.1 == T) {
      p <- p + geom_abline(slope = 1, color = "black")
    }

    if (!is.null(lim) | !is.null(breaks)) {
      p <- p +
        scale_x_continuous(limits = lim, breaks = breaks) +
        scale_y_continuous(limits = lim, breaks = breaks)
    }

    p <- p +
      theme(
        axis.text = element_text(size = 16, face = "bold"), #for ticks
        axis.title = element_text(size = 16, face = "bold", vjust = 0.3), #for the labels
        plot.title = element_text(face = "bold", size = 20, vjust = 1), #title
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 16),
        panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.title.y = element_text(vjust = 1)
      ) #legend.title
  } else {
    p <- NULL
  }

  o <- list(plot = p, r.squared = adj.r.squared, pearson.corr = pearson.corr)
  return(o)
}

getCorrelationPlots <- function(
  x,
  typeColoring = F,
  lim = NULL,
  breaks = NULL,
  point.size = 2,
  doPlot = T,
  inTitle = "rs",
  units = NULL,
  returnPlots = F,
  label = F
) {
  # if (label) {
  #   x1 <- x[,-1]
  # }else {
  #   x1 <- x
  # }

  # check if the first column is the "type" or not
  if (typeColoring) {
    x1 <- x[, -1]
  } else {
    x1 <- x
  }

  # generate a dataframe for all pairwise comparisons
  print(colnames(x1))
  y <- combn(ncol(x1), 2)
  plots <- c()
  corr <- data.frame(
    col1 = character(0),
    col2 = character(0),
    r.squared = numeric(0),
    pearson.corr = numeric(0)
  )
  # iterate over each comparison
  for (i in c(1:ncol(y))) {
    cols <- y[, i]
    # select the two columns to compare
    df <- x1[, cols]

    if (typeColoring) {
      df <- cbind(x[, "type", drop = F], df)
    }

    # if (label)
    #   df <- cbind(x[,1,drop=F],df)

    res <- getCorrelationPlot(
      df,
      typeColoring = typeColoring,
      lim = lim,
      breaks = breaks,
      point.size = point.size,
      doPlot = doPlot,
      inTitle = inTitle,
      units = units,
      label = label
    )
    plots[[i]] <- res$plot
    corr1 <- data.frame(
      col1 = colnames(x1)[cols[1]],
      col2 = colnames(x1)[cols[2]],
      r.squared = signif(res$r.squared, 3),
      pearson.corr = signif(res$pearson.corr, 3),
      stringsAsFactors = F
    )
    corr <- rbind(corr, corr1)
  }
  o <- list(plots = plots, correlations = corr)
}


## Plots correlation plots for all replicated conditions
#in case of typeColoring there is a extra column with type. The conditions are for the other columns
# x: dataframe with reads (expected first column contains the type, remaining columns contain numeric values)
# 1format: format of the output file: pdf or jpg
# inTitle: rs=rSquared, pc=pearson correlatino
# returnPlot: if True, returns the plots to the function call
correlationPlots <- function(
  x,
  Conditions,
  outputFile,
  typeColoring = F,
  lim = NULL,
  breaks = NULL,
  nrow = NULL,
  ncol = NULL,
  plotformat = "jpg",
  point.size = 2,
  doPlot = T,
  checkReplicates = T,
  inTitle = "rs",
  units = NULL,
  returnPlot = F
) {
  if (!is.null(lim) & is.null(breaks)) {
    stop(
      "If you set a value for 'lim', you also need to set a value for 'breaks'."
    )
  }

  # if (typeColoring) {
  #   # exclude the first column (the type)
  #   x1 <- x[,-1]
  # }else {
  #   x1 <- x
  # }

  # identify unique conditions
  # if (checkReplicates) {
  #   unique_conditions <- names(which(tapply(Conditions,Conditions,length)>1))
  # }else {
  #   unique_conditions <- unique(Conditions)
  # }
  unique_conditions = unique(Conditions)

  nrPlots <- 0
  # calculate the number of pairwise comparisons that can be made for each timepoint
  for (i in c(1:length(unique_conditions))) {
    nrPlots <- nrPlots +
      nrPlotsPerCond(length(grep(
        paste("^", unique_conditions[i], "$", sep = ""),
        Conditions
      )))
  }
  # define the numbe of rows and columns for the final figure
  if (is.null(nrow) | is.null(ncol)) {
    #Per row the Conditions
    nrow = length(unique_conditions)
    #mean nr of plots per conds, take ceiling
    ncol <- ceiling(nrPlots / length(unique_conditions))
  }

  plots <- c()
  corrAll <- data.frame(
    unique_conditions = character(0),
    col1 = character(0),
    col2 = character(0),
    r.squared = numeric(0)
  )
  # column_names=sub("_[^_]+$", "", names(x))
  for (i in c(1:length(unique_conditions))) {
    # select the columns from the current timepoint
    df <- x[, grepl(paste0("^", unique_conditions[i], "_r"), names(x))]
    if (typeColoring) {
      df <- cbind(x[, "type", drop = F], df, stringsAsFactors = F)
    }
    # create all correlation plots for the current timepoint
    res <- getCorrelationPlots(
      df,
      typeColoring = typeColoring,
      lim = lim,
      breaks = breaks,
      point.size = point.size,
      doPlot = doPlot,
      inTitle = inTitle,
      units = units,
      label = label
    )
    plots <- c(plots, res$plots)
    # save the pearson correlation and r^2 value for the created plots in a dataframe
    corr <- data.frame(
      unique_conditions = rep(unique_conditions[i], nrow(res$corr)),
      res$corr,
      stringsAsFactors = F
    )
    corr$col1 <- gsub(paste(unique_conditions[i], "_", sep = ""), "", corr$col1)
    corr$col2 <- gsub(paste(unique_conditions[i], "_", sep = ""), "", corr$col2)
    corrAll <- rbind(corrAll, corr)
  }

  if (doPlot) {
    # ggsave(filename = paste(outputFile,".tiff",sep=""), plot=plots,width=600, height=500, units="mm", dpi=300, compression = "lzw")
    width <- ncol * 300
    height = nrow * 300
    if (typeColoring) {
      width <- 1.5 * width
    }

    jpeg(paste(outputFile, ".jpg", sep = ""), width = width, height = height) #in px

    #Apparently one more than 10 plots, a ex1tra top layer of list is created
    # Therefore check if the first layers are plot, if not remove the first layer of list
    cl <- class(plots[[1]])
    # if (!("ggplot" %in% unlist(cl)))
    #   plots <- unlist(plots,recursive=F)
    args <- c(plots, list(nrow = nrow, ncol = ncol))
    print(length(plots))
    print(nrow)
    print(ncol)
    do.call(grid.arrange, args = args)
    dev.off()
  }

  outputFile <- gsub(paste("..", plotformat, "$", sep = ""), "", outputFile)
  write.table(
    corrAll,
    file = paste(outputFile, ".csv", sep = ""),
    row.names = F,
    sep = ";"
  )

  if (returnPlot) {
    return(plots)
  }
}


qualityControlPlots <- function(
  target_genes,
  annotation_columns = NULL,
  Conditions = NULL,
  basename,
  correlation = T,
  heatmap = T,
  sizeFactors = NULL,
  heatmapFitType = "local",
  ncol = NULL,
  nrow = NULL,
  typeColoring = T,
  distrPlotYaxisMax = NULL,
  plotformat = "pdf",
  other = NULL,
  ctrlDistrNRow = NULL,
  ctrlDistrNCol = NULL,
  ctrlDistrFrom = -10,
  ctrlDistrTo = 3,
  ctrlWidthSizeFactor = 7,
  ctrlTypes = c("p", "n"),
  paired = F,
  checkReplicates = T,
  heatmapWidth = NULL,
  heatmapMargin = c(13, 13),
  inclRep = T,
  corrLim = NULL,
  corrBreaks = NULL
) {
  if (is.null(Conditions)) {
    Conditions <- sub(
      "_[^_]+$",
      "",
      colnames(target_genes)[-c(1:annotation_columns)]
    )
  }

  plots <- list()

  annotations <- target_genes[, c(1:annotation_columns)]
  normCounts <- target_genes[, -c(1:annotation_columns)]

  # create correlation plots of the log10-transformed reads
  if (correlation) {
    # If sizeFactors are not given, it will be generated. sizeFactors based on colSums no longer supported.

    # calculate the log10 of the reads
    lognormCounts <- log10(normCounts)

    if (typeColoring) {
      #add as first column type
      lognormCounts <- data.frame(
        type = annotations$type,
        lognormCounts,
        stringsAsFactors = F
      )
    }

    outputFile <- "log10_correlation_plots"
    correlationPlots(
      x = lognormCounts,
      Conditions,
      outputFile = paste(basename, "correlation_plots", sep = "_"),
      nrow = Nrrows,
      ncol = Nrcols,
      typeColoring = typeColoring,
      plotformat = plotformat,
      units = "log10Counts",
      checkReplicates = checkReplicates,
      returnPlot = F,
      lim = corrLim,
      breaks = corrBreaks
    )
  }

  # create a heatmap
  if (heatmap) {
    bioCounts <- target_genes[, -c(1:annotation_columns)]
    if (inclRep) {
      rep <- sapply(strsplit(colnames(bioCounts), "_"), function(x) {
        paste(x[length(x)], collapse = "_")
      })
      # simple dataset object just for normalization
      colData <- data.frame(
        Conditions = as.factor(Conditions),
        rep = as.factor(rep)
      )
    } else {
      colData <- data.frame(Conditions = as.factor(Conditions))
    }
    # create a new deseq data set with the normalised counts
    dds <- DESeqDataSetFromMatrix(
      countData = bioCounts,
      colData = colData,
      design = ~Conditions
    )

    # set the size factor for each condition to 1
    sizeFactors(dds) <- rep(1, ncol(bioCounts)) #just used the data as presented in target_genes
    # set the width and height of the heatmap
    if (is.null(heatmapWidth)) {
      heatmapWidth <- max(length(Conditions) / 4, 10)
    }

    ## cds <- newCountDataSet(countData=counts,Conditions=Conditions)
    ## cds <- estimateSizeFactors(cds)
    pdf(
      paste(basename, "heatmap.pdf", sep = "_"),
      width = heatmapWidth,
      height = heatmapWidth
    )
    showSampleHeatmap(dds, margin = heatmapMargin, inclRep = inclRep)
    dev.off()
  }
}
