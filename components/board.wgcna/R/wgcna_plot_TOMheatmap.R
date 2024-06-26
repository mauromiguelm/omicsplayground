##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_TOMheatmap_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_TOMheatmap_server <- function(id,
                                         wgcna.compute,
                                         labels2rainbow,
                                         power,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    TOMplot.RENDER <- shiny::reactive({
      out <- wgcna.compute()
      net <- out$net
      datExpr <- out$datExpr
      geneTree <- net$dendrograms[[1]]
      moduleColors <- labels2rainbow(out$net)
      MEs <- out$net$MEs

      ## Calculate topological overlap anew: this could be done
      ## more efficiently by saving the TOM calculated during
      ## module detection, but let us do it again here.


      power <- as.numeric(power())
      dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = power)
      rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)

      nSelect <- 999999
      nSelect <- 400
      ## For reproducibility, we set the random seed
      set.seed(10)
      select <- head(1:ncol(dissTOM), nSelect)
      selectTOM <- dissTOM[select, select]
      ## There’s no simple way of restricting a clustering tree
      ## to a subset of genes, so we must re-cluster.
      #
      selectTree <- hclust(as.dist(selectTOM), method = "average")
      selectColors <- moduleColors[select]
      ## Taking the dissimilarity to a power, say 10, makes the plot
      ## more informative by effectively changing the color palette;
      ## setting the diagonal to NA also improves the clarity of the
      ## plot
      plotDiss <- selectTOM^7
      diag(plotDiss) <- NA
      myheatcol <- gplots::colorpanel(250, "red", "orange", "lemonchiffon")
      myheatcol <- gplots::colorpanel(250, "lemonchiffon", "orange", "red")

      par(oma = c(2, 0, 0, 0))
      plotly::layout(
        matrix(c(
          0, 0, 5, 0,
          0, 0, 2, 0,
          4, 1, 3, 6
        ), nr = 3, byrow = T),
        widths = c(2.3, 0.5, 10, 1.8),
        heights = c(2.3, 0.5, 10)
      )

      WGCNA::TOMplot(
        plotDiss, selectTree, selectColors,
        col = myheatcol,
        setLayout = FALSE,
        main = NULL
      )

      ## add color legend
      frame()
      me.names <- colnames(MEs)
      me.nr <- as.integer(sub("ME", "", me.names))
      ii <- order(me.nr)
      label.colors <- labels2rainbow(net)
      me.colors <- label.colors[!duplicated(names(label.colors))]
      me.colors <- me.colors[as.character(me.nr)]

      legend(-0.1, 1,
        legend = me.names[ii], fill = me.colors[ii],
        cex = 1.2, bty = "n", x.intersp = 0.5
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = TOMplot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 90),
      add.watermark = watermark
    )
  })
}
