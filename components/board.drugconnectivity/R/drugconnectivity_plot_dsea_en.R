drugconnectivity_plot_dsea_en_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)
    info_text = "<strong>Drug connectivity</strong> correlates your signature with known
     drug profiles from the L1000 database, and shows similar and opposite profiles by
      running the GSEA algorithm on the drug profile correlation space."

    opts = shiny::tagList()

    PlotModuleUI(
        ns("pltmod"),
        title = "Drug connectivity",
        label = "a",
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,
        info.text = info_text,
        options = opts,
        download.fmt=c("png","pdf","csv"),
        height = height
    )
}


drugconnectivity_plot_dsea_en_server <- function(id,
                                                pgx,
                                                getActiveDSEA,
                                                dmethod,
                                                dsea_table,
                                                watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        dbg("[dataview_expressionplot_server] created!")

        plot_data <- shiny::reactive({

            ngs <- pgx
            if(is.null(ngs$drugs)) return(NULL)
            shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))

            dbg("[dsea_enplots.RENDER] called!")

      
            dsea <- getActiveDSEA()
            dt <- dsea$table

            ## filter with table selection/search
            ii  <- dsea_table$rows_selected()
            jj  <- dsea_table$rows_all()
            shiny::req(jj)  ## must have non-empty table

            if(length(ii)>0) {
                dt <- dt[ii,,drop=FALSE]
            }
            if(length(ii)==0 && length(jj)>0) {
                dt <- dt[jj,,drop=FALSE]
            }
            if(nrow(dt)==0) return(NULL)

            return(dt)
        })


        plot.RENDER <- function() {

            pd  <- plot_data()
           
            shiny::req(pd)
            ## rank vector for enrichment plots
            
            dmethod <- dmethod()
            dsea <- getActiveDSEA()
            rnk <- dsea$stats
            if(length(rnk)==0) return(NULL)

            ## ENPLOT TYPE
            if(nrow(pd)==1) {
                par(oma=c(1,1,1,1))
                par(mfrow=c(1,1), mar=c(4,4,1.1,2), mgp=c(2.3,0.9,0))
                lab.cex = 1
                xlab="Rank in ordered dataset"
                ylab="Rank metric"
                nc=1
            } else {
                pd <- head(pd, 16)
                lab.cex = 0.75
                xlab=''
                ylab=""
                nc = ceiling(sqrt(nrow(pd)))
                par(oma=c(0,1.6,0,0))
                par(mfrow=c(nc,nc), mar=c(0.3,1.0,1.3,0), mgp=c(1.9,0.6,0))
            }

            i=1
            for(i in 1:nrow(pd)) {
                dx <- rownames(pd)[i]
                dx
                gmtdx <- grep(dx,names(rnk),fixed=TRUE,value=TRUE)  ## L1000 naming
                length(gmtdx)
                ##if(length(gmtdx) < 3) { frame(); next }
                dx1 <- substring(dx,1,26)
                par(cex.axis=0.001)
                if(i%%nc==1) par(cex.axis=0.98)
                suppressWarnings(
                    gsea.enplot( rnk, gmtdx, main=dx1, cex.main=1.2,
                                xlab=xlab, ylab=ylab)
                )
                nes <- round(pd$NES[i],2)
                qv  <- round(pd$padj[i],3)
                tt <- c( paste("NES=",nes), paste("q=",qv) )
                legend("topright", legend=tt, cex=0.8, y.intersp=0.85, bty='n')
                if(i%%nc==1 && nrow(pd)>1) {
                    mtext('rank metric', side=2, line=1.8, cex=lab.cex)
                }
            }



        }

        modal_plot.RENDER <- function() {
            plot.RENDER()
        }

        PlotModuleServer(
            "pltmod",
            plotlib = "base",
            plotlib2 = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            renderFunc = shiny::renderPlot,
            renderFunc2 = shiny::renderPlot,
            ##renderFunc = shiny::renderCachedPlot,
            ##renderFunc2 = shiny::renderCachedPlot,
            res = c(90,170)*1,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}
