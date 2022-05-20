##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadBoard <- function(id,
                        pgx_dir,
                        pgx,
                        auth,
                        limits = c(
                          "samples" = 1000,
                          "comparisons" = 20,
                          "genes" = 20000,
                          "genesets" = 10000,
                          "datasets" = 10
                        ),
                        enable_upload = TRUE,
                        enable_save = TRUE,
                        enable_userdir = TRUE
                        )
{
  moduleServer(id, function(input, output, session) 
  {
    ns <- session$ns ## NAMESPACE
    dbg("[UploadBoard] >>> initializing UploadBoard...")

    message("[UploadBoard] pgx_dir = ",pgx_dir)
    message("[UploadBoard] enable_upload = ",enable_upload)    

    if(!enable_upload) return(NULL)

    loadedDataset <- shiny::reactiveVal(0)  ## counts/trigger dataset upload
    height = 720
    
    ##================================================================================
    ##================================= info =========================================
    ##================================================================================
    
    shiny::observeEvent( input$module_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Upload data</strong>"),
            shiny::HTML(module_infotext),
            easyClose = TRUE, size = "l" ))
    })

    module_infotext = paste0(
        'Under the <b>Upload data</b> panel users can upload their transcriptomics and proteomics data to the platform. The platform requires 3 data files as listed below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons as contrasts (contrasts.csv). It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, row names and column names match for all files. On the left side of the panel, users need to provide a unique name and brief description for the dataset while uploading. N.B. Users can now create contrasts from the platform itself, so the contrasts.csv file is optional.

<br><br>
<ol>
<li>counts.csv: Count/expression file with gene on rows, samples as columns.
<li>samples.csv: Samples file with samples on rows, phenotypes as columns.
<li>contrasts.csv: Contrast file with conditions on rows, contrasts as columns.
</ol>

<br><br><br>
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>

')
    
    ##================================================================================
    ##=================================== UI =========================================
    ##================================================================================
    
    observeEvent( input$advanced_mode, {
        message("[upload_server] advanced_mode = ",input$advanced_mode)
        if(input$advanced_mode) {
            showTab("wiz", "page2")
        } else {
            hideTab("wiz", "page2")
        }
        
    })

    
    ##================================================================================
    ##============================== modules =========================================
    ##================================================================================

    WizardDialogServer("wiz", 4)
    
    ## Some 'global' reactive variables used in this file
    files <- shiny::reactiveValues(
      counts    = NULL,
      samples   = NULL,
      contrasts = NULL,
      pgx = NULL,
      meta = NULL
    )
       
    uploaded <- UploadModuleServer(
      id = "upload",
      FILES = FILES,        ## for example.zip
      height = height,
      limits = limits,
      show_tables = reactive(input$show_tables)
    )

    observeEvent( uploaded$counts.csv(), {
      dbg("[upload_server] uploaded$counts.csv() = ",str(uploaded$counts.csv()))
      files$counts <- uploaded$counts.csv()
    })

    observeEvent( uploaded$samples.csv(), {
      dbg("[upload_server] uploaded$samples.csv() = ",str(uploaded$samples.csv()))      
      files$samples <- uploaded$samples.csv()
    })

    observeEvent( uploaded$contrasts.csv(), {
      contrasts <- uploaded$contrasts.csv()
      samples   <- uploaded$samples.csv()      
      contrasts.new <- pgx.convertContrastsExplicit(contrasts, samples)
      files$contrasts <- contrasts.new
    })

    observeEvent( uploaded$meta(), {
      files$meta <- uploaded$meta()
    })
    
    ##correctedX <- shiny::reactive({
    ## normalized_counts <- NormalizeCountsServerRT(
    ##   id = "normalize",
    ##   counts = shiny::reactive(files$counts),
    ##   height = height
    ## )
    
    ##correctedX <- shiny::reactive({
    correctedX <- BatchCorrectServer(
      id = "batchcorrect",
      X = shiny::reactive(files$counts),
      ##X = normalized_counts,  ## NOT YET!!!!
      is.count = TRUE,
      pheno = shiny::reactive(files$samples),
      height = height
    )
    
    ##================================================================================
    ##======================== REACTIVE COUNTS =======================================
    ##================================================================================

    corrected_counts <- shiny::reactive({
      counts <- NULL
      dbg("[UploadModule::corrected_counts] reacted!\n")                
      advanced_mode <- ( length(input$advanced_mode)>0 &&
                           input$advanced_mode[1]==1 )
      if(advanced_mode) {
        message("[UploadModule::corrected_counts] using CORRECTED counts\n")
        out <- correctedX()
        counts <- pmax(2**out$X-1, 0)
      } else {
        message("[UploadModule::corrected_counts] using UNCORRECTED counts\n")
        counts <- files$counts
      }
      counts
    })

    ##mkContrast <- shiny::reactive({
    modified_ct <- MakeContrastServerRT(
      id = "makecontrast",
      phenoRT = shiny::reactive(files$samples),
      contrRT = shiny::reactive(files$contrasts),
      countsRT = shiny::reactive(files$counts),
      ##countsRT = corrected_counts,
      height = height
    )
    
    batch_vectors <- shiny::reactive({
      dbg("batch_vectors reactive")
      correctedX()$B
    })
    
    shiny::observeEvent( modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      ##
      dbg("[observe:modified_ct()] reacted...")                
      modct <- modified_ct()
      dbg("[observe:modified_ct()] dim(modct$contr) = ",dim(modct$contr))
      files$contrasts <- modct$contr
      files$samples   <- modct$pheno
    })

    computed_pgx  <- ComputePgxServer(
      id = "compute", 
      countsRT = shiny::reactive(files$counts),
      ##countsRT = corrected_counts,
      samplesRT = shiny::reactive(files$samples),
      contrastsRT = shiny::reactive(files$contrasts),
      batchRT = batch_vectors, 
      metaRT = shiny::reactive(files$meta),                
      ## enable_button = upload_ok,
      enable_button = reactive(TRUE),      
      alertready = FALSE,
      FILES = FILES,
      ##pgx.dirRT = shiny::reactive(getPGXDIR()),
      pgx.dirRT = getPGXDIR,      
      max.genes = as.integer(limits["genes"]),
      max.genesets = as.integer(limits["genesets"]),
      max.datasets = as.integer(limits["datasets"]),
      height = height
    )          
    
    ##================================================================================
    ##============================== Finish up =======================================
    ##================================================================================

    shiny::observeEvent( computed_pgx(), {
      
      dbg("[observe::computed_pgx] uploaded PGX detected!")
      new_pgx <- computed_pgx()
      
      dbg("[observe::computed_pgx] initializing PGX object")
      new_pgx <- pgx.initialize(new_pgx)
      
      ## update Session PGX
      dbg("[UploadBoard@load_react] **** copying current pgx to session.pgx  ****")        
      for(i in 1:length(new_pgx)) {
        pgx[[names(new_pgx)[i]]] <- new_pgx[[i]]
      }
      
      DT::selectRows(proxy = DT::dataTableProxy(ns("pgxtable")), selected=NULL)            
      
      savedata_button <- NULL
      if(enable_save) {                
        
        dbg("[UploadBoard] observeEvent:savedata reacted")        
        ## -------------- save PGX file/object ---------------
        pgxname <- sub("[.]pgx$","",new_pgx$name)
        pgxname <- gsub("^[./-]*","",pgxname)  ## prevent going to parent folder
        pgxname <- paste0(gsub("[ \\/]","_",pgxname),".pgx")
        pgxname
        
        pgxdir  <- getPGXDIR()
        fn <- file.path(pgxdir,pgxname)
        fn <- iconv(fn, from = '', to = 'ASCII//TRANSLIT')
        
        ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ## Note: Currently we use 'ngs' as object name but want to go
        ## towards 'pgx' as standard name. Actually saving as RDS
        ## should be better.
        ngs=new_pgx
        save(ngs, file=fn)
        
        remove(ngs)
        remove(new_pgx)
                
        message("[UploadBoard::@savedata] updating PGXINFO")                    
        pgx.initDatasetFolder(pgxdir, force=FALSE, verbose=TRUE)
        ## reload_pgxdir(reload_pgxdir()+1)
      }
      
      ## shiny::removeModal()
      msg1 <- "<b>Ready!</b>"
      ##beepr::beep(sample(c(3,4,5,6,8),1))  ## music!!
      beepr::beep(10)  ## short beep
      
      if(enable_save) {
        msg1 <- "<b>Ready!</b><br>Your data is ready and has been saved in your library. You can now start exploring your data."
      } else {
        msg1 <- "<b>Ready!</b><br>Your data is ready. You can now start exploring your data."
      }
      loadedDataset(loadedDataset()+1)  ## notify new data uploaded
      
      showModal(
        modalDialog(
          HTML(msg1),
          title = NULL,
          size = "s",
          footer = tagList(
            ## savedata_button,
            ## shiny::actionButton(ns("sharedata"), "Share with others", icon=icon("share-alt")),
            modalButton("Start!")
          )
        ))

        on.exit({
            session$sendCustomMessage(
                "show-tabs",
                list()
            )
        })
        
        ## updateTabsetPanel(session, "tabs",  selected = "Datasets")
    })
    
    ##================================================================================
    ##======================== HELPER FUNCTIONS ======================================
    ##================================================================================
    
    getPGXDIR <- shiny::reactive({
        ##reload_pgxdir()  ## force reload

        email="../me@company.com"
        email <- auth$email()
        email <- gsub(".*\\/","",email)
        pdir  <- pgx_dir  ## from module input

        ##USERDIR=FALSE
        if(enable_userdir) {
            pdir <- paste0(pdir,"/",email)
            if(!is.null(email) && !is.na(email) && email!="") pdir <- paste0(pdir,'/')
            if(!dir.exists(pdir)) {
                dbg("[LoadingBoard:getPGXDIR] userdir does not exists. creating pdir = ",pdir)
                dir.create(pdir)
                dbg("[LoadingBoard:getPGXDIR] copy example pgx")                
                file.copy(file.path(pgx_dir,"example-data.pgx"),pdir)
            }
        }
        pdir
    })

    ##------------------------------------------------
    ## Board return object
    ##------------------------------------------------
    res <- list(
        loaded = loadedDataset
    )
    return(res)
  })
}
