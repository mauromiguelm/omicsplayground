drugconnectivity_plot_dsea_moa_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)
    info_text = "<strong>Drug connectivity</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space."
    
    opts = shiny::tagList()
    
    PlotModuleUI(
        ns("pltmod_dsea_moa"),
        title = "Mechanism of action",
        label = "c",
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = info_text,
        options = opts,
        download.fmt=c("png","pdf","csv"),
        height = height
    )
}


drugconnectivity_plot_dsea_moa_server <- function(id,
                                                pgx,
                                                watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        plot_data <- shiny::reactive({

            browser()

            shiny::req(pgx$X)

            ##### new data
            ngs <- pgx
            shiny::req(ngs)
            
            if(is.null(ngs$drugs)) return(NULL)
            shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))    

            dbg("[dsea_moaplot.RENDER] reacted")
            
            res <- getMOA()
            ntop = 16
            jj <- unique(c(head(order(-res$NES),ntop),tail(order(-res$NES),ntop)))
            moa.top <- res$NES[jj]
            names(moa.top) <- res$pathway[jj]
            par(mfrow=c(2,1), mar=c(4,3.5,0.1,0), mgp=c(1.7,0.65,0))
            
            return(moa.top)
        })


        plot.RENDER <- function() {

            browser()
            
            pd  <- plot_data()
            shiny::req(pd)
            browser()

            barplot(pd, horiz=FALSE, las=3,
                ylab="enrichment  (NES)",
                cex.names=0.96 )

        }
        
        modal_plot.RENDER <- function() {
            plot.RENDER()
        }
        
        PlotModuleServer(
            "pltmod_dsea_moa",
            plotlib = "base",
            plotlib2 = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            renderFunc = shiny::renderPlot,
            renderFunc2 = shiny::renderPlot,
            res = c(90,170)*1,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}