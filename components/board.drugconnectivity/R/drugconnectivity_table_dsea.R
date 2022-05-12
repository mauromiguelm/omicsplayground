##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


drugconnectivity_table_dsea_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("tbl"))  
}


drugconnectivity_table_dsea_server <- function(id,
                                                pgx,
                                                getActiveDSEA
                                                )
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns
    
    table_data <- shiny::reactive({

        ngs <- pgx
        shiny::req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        dsea <- getActiveDSEA()
        shiny::req(dsea)
        res <- dsea$table
        res$moa <- shortstring(res$moa,60)
        res$target <- shortstring(res$target,30)
        res$drug   <- shortstring(res$drug,60)
        colnames(res) <- sub("moa","MOA",colnames(res))
    }) 

    
    table.RENDER <- function() {

        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=NULL),
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip', buttons = c('copy','csv','pdf'),
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = '70vh',
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "NES",
                                background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    }

    
    info_text="<b>Enrichment table.</b> Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug."

    ##--------buttons for table
    opts = shiny::tagList(
        withTooltip(shiny::checkboxInput(ns('dseatable_filter'),'only annotated drugs',TRUE),
               "Show only annotated drugs.")
    ) 
    
    shiny::callModule(
      tableModule, "tbl",
      label="b",
      func = table.RENDER,
      csvFunc = enrichment_data,
      options = opts,
      title = "Enrichment table",
      filename = "enrichment.csv",
      info.text = info_text
    )

  })  ## end of moduleServer
} ## end of server

