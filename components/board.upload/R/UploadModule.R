##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##================================================================================
## Upload data Module
##================================================================================

UploadModuleUI <- function(id) {
  ns <- shiny::NS(id)

  csv.info <- "Upload your counts as 'counts.csv' and sample information as 'samples.csv'. The file format must be comma-separated-values (CSV) text. Check that dimensions, rownames and column names match for all files. You can download a zip file with example files here: exampledata.zip. You can upload a maximum of 100 datasets (with each up to 2000 samples and 100 comparisons)."
  
  fluidPage(
    shiny::fluidRow(
      shiny::column(3,
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Upload CSV",
            br(),
            ##h4("Upload from CSV file"),
            csv.info,
            br(), br(),
            shiny::fileInput(ns("upload_files"),"Choose file:"),
            ##shinyWidgets::prettySwitch(ns("load_example"), "load example data", value=FALSE)
            actionButton(ns("load_example"), "load example data", class="btn-outline-success")          
          ),
          shiny::tabPanel(
            "Retrieve GEO",
            br(),
            h4("Retrieve from GEO"),          
            "Please enter the GEO id. We try to download the data automatically. Warning: this does not alwyas work and you may need to download the data manually from GEO.",
            br(), br(),
            shiny::textInput(ns("geoID"), "Enter GEO id:"),
            shiny::actionButton(ns("geoLoad"),"Retrieve from GEO")
          )
        ) ## end of tabsetPanel
      ), 
      column(9,
        CheckDataUI(ns("chkdata"))
      )
    ),      
    shiny::fluidRow(
        shiny::conditionalPanel(
            "input.show_tables == true", ns=ns,
            tabsetPanel(
                tabPanel(
                    "Counts",
                    div(DT::dataTableOutput(ns("countstable"), width="100%"),
                        style = "height:250px; overflow-y:scroll; overflow-x:scroll;")
                ),
                tabPanel(
                    "Samples",
                    div(
                        DT::dataTableOutput(ns("sampletable"), width="100%"),
                        style = "height:250px; overflow-y:scroll; overflow-x:scroll;"
                    )
                ),
                tabPanel(
                    "Contrasts",
                    div(
                        DT::dataTableOutput(ns("contrasttable"), width="100%"),
                        style = "height:250px; overflow-y:scroll; overflow-x:scroll;"        
                    )
                )
            )  ## end of tabsetPanel
        ) ## end of conditionalPanel
    )  ## end of fluidRow

  )  ## end of fluidPage
    
}

UploadModuleServer <- function(id, 
                               FILES,
                               height = 720,
                               limits = c(
                                   samples = 100,
                                   comparisons = 20,
                                   genes = 20000,
                                   genesets = 10000,
                                   datasets = 10 ),
                               show_tables = reactive(TRUE)
                               )
{
    shiny::moduleServer(
        id,
        function(input, output, session) {

            ns <- session$ns
            ## ns <- shiny::NS(id)
                      
          
            output$downloadExampleData <- shiny::downloadHandler(
                filename = "exampledata.zip",
                content = function(file) {
                    zip = file.path(FILES,"exampledata.zip")
                    file.copy(zip,file)
                }
            )
            
            upload_info = "<h4>User file upload</h4><p>Please prepare the data files in CSV format as listed below. It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, rownames and column names match for all files. You can download a zip file with example files here: EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>."
            DLlink = shiny::downloadLink(ns("downloadExampleData"),"exampledata.zip")
            upload_info = sub("EXAMPLEZIP", DLlink, upload_info)
            
            ##upload_filetypes = c("text/csv","text/comma-separated-values,text/plain",".csv")
            upload_filetypes = c(".csv",".pgx")                

            limits0 <- paste(limits["datasets"],"datasets (with each up to",
                             limits["samples"],"samples and",
                             limits["comparisons"],"comparisons)")
            upload_info = sub("LIMITS", limits0, upload_info)
            
            ##========================================================================
            ##================================= UI ===================================
            ##========================================================================

          ## Some 'global' reactive variables used in this file
            uploaded <- shiny::reactiveValues(
              counts.csv = NULL,
              samples.csv = NULL,
              contrasts.csv =  NULL,
              pgx = NULL,
              meta = NULL
            )

          
          CheckDataModule(
            "chkdata",
            countsRT    = shiny::reactive(uploaded$counts.csv), 
            samplesRT   = shiny::reactive(uploaded$samples.csv),
            contrastsRT = shiny::reactive(uploaded$contrasts.csv)
          ) 

            ##=====================================================================
            ##============================== TABS =================================
            ##=====================================================================            

            if(0) {
                ## Hide/show tabpanels upon available data like a wizard dialog
                ##shiny::observeEvent( names(uploaded), {
                shiny::observe({                                    
                    has.upload <- Vectorize(function(f) {
                        ( f %in% names(uploaded) && !is.null(nrow(uploaded[[f]])))
                    })
                    need2 = c("counts.csv","samples.csv")                
                    need3 = c("counts.csv","samples.csv","contrasts.csv")
                    if(all(has.upload(need3))) {
                        shiny::showTab("tabs", "Contrasts")
                        shiny::showTab("tabs", "Compute")
                        if(input$advanced_mode) {
                            shiny::showTab("tabs", "Normalize")                            
                            shiny::showTab("tabs", "BatchCorrect")
                        }
                        ##shinyjs::enable(selector = '.navbar-nav a[data-value="Contrasts"')
                        ##shinyjs::enable(selector = '.navbar-nav a[data-value="Compute"')        
                    } else if(all(has.upload(need2))) {
                        if(input$advanced_mode) {
                            shiny::showTab("tabs", "Normalize")
                            shiny::showTab("tabs", "BatchCorrect")
                        }
                        shiny::showTab("tabs", "Contrasts")
                        shiny::hideTab("tabs", "Compute")
                        ##shinyjs::ensable(selector = '.navbar-nav a[data-value="Contrasts"')
                        ##shinyjs::disable(selector = '.navbar-nav a[data-value="Compute"')                        
                    } else {
                        shiny::hideTab("tabs", "Normalize")
                        shiny::hideTab("tabs", "BatchCorrect")
                        shiny::hideTab("tabs", "Contrasts")
                        shiny::hideTab("tabs", "Compute")
                        ##shinyjs::disable(selector = '.navbar-nav a[data-value="Contrasts"')
                        ##shinyjs::disable(selector = '.navbar-nav a[data-value="Compute"')                        
                    }                    
                })
            }
            
            ##=====================================================================
            ##======================= UI OBSERVERS ================================
            ##=====================================================================            
            
            shiny::observeEvent( input$advanced_mode, {
                if(input$advanced_mode) {
                  shiny::showTab("tabs", "Normalize")   ## NOT YET!!!
                    shiny::showTab("tabs", "BatchCorrect")
                } else {
                    shiny::hideTab("tabs", "Normalize")
                    shiny::hideTab("tabs", "BatchCorrect")
                    message("[UploadModule] hiding tab:",paste0(ns("wiz"),"-page2"))
                    shiny::hideTab("wiz", "page2")                    
                }
            })
            
            ##=====================================================================
            ##================== DATA LOADING OBSERVERS ===========================
            ##=====================================================================            
            
            checkDupRows <- function(F, fn, alert=TRUE) {
                rn <- setdiff(F[[1]],c("","NA",NA))
                sum(duplicated(rn))
                t1 <- sum(duplicated(rn))==0
                if(alert && !t1) {
                    shinyalert::shinyalert(
                                    title= paste(fn,"file has duplicated row names"),
                                    text = "Please correct input file",
                                    type = "error")
                }
                return(t1)
            }
            checkEmptyRows <- function(F, fn, alert=TRUE) {
                t1 <- nrow(F)>0
                if(alert && !t1) {
                    shinyalert::shinyalert(
                                    title= paste(fn,"file is empty!"),
                                    text = "Please correct input file",
                                    type = "error")
                }
                return(t1)
            }

            checkDupCols <- function(F, fn, alert=TRUE) {            
                t1 <- sum(duplicated(colnames(F))) == 0
                if(alert && !t1) {
                    shinyalert::shinyalert(
                                    title= paste(fn,"file has duplicated column names"),
                                    text = "Please correct input file",
                                    type = "error")
                    return(FALSE)
                }
                return(t1)                
            }

            checkMaxSamples <- function(F, alert=TRUE) {                        
                MAXSAMPLES = as.integer(limits["samples"])
                t1 <- (ncol(F)-1) <= MAXSAMPLES
                if(alert && !t1) {
                    shinyalert::shinyalert(
                                    title= "Too many samples",
                                    text = paste("Please decrease number of samples. You can upload maximum",MAXSAMPLES,"samples"),  
                                    type = "error")
                }
                t1
            }

            checkMaxContrasts <- function(F, alert=TRUE) {                        
                MAXCONTRASTS = as.integer(limits["comparisons"])
                t1 <- (ncol(F)-1) <= MAXCONTRASTS
                if(alert && !t1) {
                    shinyalert::shinyalert(
                        title= "Too many contrasts",
                        text = paste("Please decrease number of contrasts. You can upload maximum",MAXCONTRASTS,"contrasts"),  
                        type = "error")
                    return(FALSE)
                }
                t1
            }           

            checkCountsCSV <- function(fn) {
                F  <- data.table::fread(fn)
                t1 <- checkDupRows(F, fn, alert=FALSE)
                t2 <- checkEmptyRows(F, fn)                
                if(!t1) {
                    shinyalert::shinyalert(
                       title= paste(fn,"file has duplicated rows"),
                       text = "Warning. Duplicated intensities will be summed (linear).",
                       type = "warning")
                }
                t3 <- checkDupCols(F, fn)                                
                return(t2 && t3)
            }

            checkSamplesCSV <- function(fn) {
                F <- data.table::fread(fn,header=TRUE)
                t1 <- checkDupRows(F, fn)
                t2 <- checkEmptyRows(F, fn)
                t3 <- checkDupCols(F, fn)
                t4 <- checkMaxSamples(F)
                return(t1 & t2 & t3 & t4)                
            }

            checkContrastsCSV <- function(fn) {
                F  <- data.table::fread(fn)
                t1 <- checkDupRows(F, fn)
                t2 <- checkEmptyRows(F, fn)                
                t3 <- checkDupCols(F, fn)
                checkMaxContrasts                
                vs.names <- colnames(F)[-1]
                t4 <- all(grepl("_vs_",vs.names))
                if(!t4) {
                    shinyalert::shinyalert(
                        title= paste(fn,"file has errors"),
                        text = "Contrast names must include '_vs_'. Must be of the form 'MAIN_vs_REF' or 'VAR:MAIN_vs_REF'. Please correct contrast names.",
                        type = "error")
                    return(FALSE)                                        
                }                
                return(t1 & t2 & t3 & t4)                
            }
            
            ##------------------------------------------------------------------
            ## Observer for uploading data files using fileInput widget.
            ##
            ## Reads in the data files from the file names, checks and
            ## puts in the reactive values object 'uploaded'. Then
            ## uploaded should trigger the computePGX module.
            ## ------------------------------------------------------------------
            shiny::observeEvent( input$upload_files, {
                
                message("[upload_files] >>> reading uploaded files")
                message("[upload_files] upload_files$name=",input$upload_files$name)
                message("[upload_files] upload_files$datapath=",input$upload_files$datapath)
                
                ##for(i in 1:length(uploaded)) uploaded[[i]] <- NULL
                uploaded$pgx <- NULL
                uploaded$last_uploaded <- NULL
                
                ## read uploaded files
                pgx.uploaded <- any(grepl("[.]pgx$",input$upload_files$name))
                matlist <- list()

                if(pgx.uploaded) {
                    
                    message("[upload_files] PGX upload detected")
                    
                    ## If the user uploaded a PGX file, we extract the matrix
                    ## dimensions from the given PGX/NGS object. Really?
                    ##
                    i <- grep("[.]pgx$",input$upload_files$name)
                    load(input$upload_files$datapath[i])  ## load NGS/PGX                    
                    ##matlist[["counts.csv"]] <- ngs$counts
                    ##matlist[["samples.csv"]] <- type.convert(ngs$samples)
                    ##matlist[["contrasts.csv"]] <- ngs$model.parameters$exp.matrix
                    uploaded$pgx <- ngs
                    
                } else {

                    ## If the user uploaded CSV files, we read in the data
                    ## from the files.
                    ##
                    message("[upload_files] getting matrices from CSV")

                    ii <- grep("csv$",input$upload_files$name)
                    ii <- grep("sample|count|contrast|expression",
                               input$upload_files$name, ignore.case=TRUE)
                    if(length(ii)==0) return(NULL)
                    
                    inputnames  <- input$upload_files$name[ii]
                    uploadnames <- input$upload_files$datapath[ii]
                    
                    if(length(uploadnames)>0) {
                        i=1
                        for(i in 1:length(uploadnames)) {
                            fn1 <- inputnames[i]
                            fn2 <- uploadnames[i]
                            matname <- NULL
                            df <- NULL
                            if(grepl("count",fn1, ignore.case=TRUE)) {
                                ## allows duplicated rownames
                                df0 <- read.as_matrix(fn2)
                                if(TRUE && any(duplicated(rownames(df0)))) {
                                  ndup <- sum(duplicated(rownames(df0)))
                                  shinyWidgets::sendSweetAlert(
                                    session=session,
                                    title = "Duplicated gene names",
                                    text = paste("Your counts matrix has",ndup,"duplicated gene names.\nCounts of those genes will be merged."),
                                    type = "warning",
                                    btn_labels = "OK",
                                    closeOnClickOutside = FALSE,
                                  )
                                }
                                
                                if(nrow(df0)>1 && NCOL(df0)>1) {
                                  df <- as.matrix(df0)
                                  matname <- "counts.csv"
                                }
                                
                            } else if(grepl("expression",fn1,ignore.case=TRUE)) {
                                ## allows duplicated rownames
                                df0 <- read.as_matrix(fn2)
                                if(TRUE && any(duplicated(rownames(df0)))) {
                                  ndup <- sum(duplicated(rownames(df0)))                                    
                                  shinyWidgets::sendSweetAlert(
                                    session=session,
                                    title = "Duplicated gene names",
                                    text = paste("Your counts matrix has",ndup,"duplicated gene names.\nCounts of those genes will be merged."),
                                    type = "warning",
                                    btn_labels = "OK",
                                    closeOnClickOutside = FALSE,
                                  )
                                }
                                if(nrow(df0)>1 && NCOL(df0)>1) {
                                    df <- as.matrix(df0)
                                    message("[UploadModule::upload_files] converting expression to counts...")
                                    df <- 2**df
                                    matname <- "counts.csv"
                                }
                                
                            } else if(grepl("sample",fn1,ignore.case=TRUE)) {
                                df0 <- read.as_matrix(fn2)
                                if(any(duplicated(rownames(df0)))) {
                                  dup.rows <- rownames(df0)[which(duplicated(rownames(df0)))]
                                  msg <- paste("Your samples file has duplicated entries: ", 
                                               dup.rows, ". This is not allowed, please correct.")
                                  shinyWidgets::sendSweetAlert(
                                    session=session,
                                    title = "Duplicated sample name",
                                    text = msg,
                                    type = "error",
                                    btn_labels = "OK",
                                    ##btn_colors = "red",
                                    closeOnClickOutside = FALSE,
                                  )
                                  
                                } else if(nrow(df0)>1 && NCOL(df0)>=1) {
                                  df <- as.data.frame(df0)
                                  matname <- "samples.csv"
                               }
                              
                            } else if(grepl("contrast",fn1,ignore.case=TRUE)) {
                                df0 <- read.as_matrix(fn2)
                                if(any(duplicated(rownames(df0)))) {
                                  dup.rows <- rownames(df0)[which(duplicated(rownames(df0)))]
                                  msg <- paste("Your contrasts file has duplicated entries: ", 
                                               dup.rows, ". This is not allowed, please correct.")
                                  shinyWidgets::sendSweetAlert(
                                    session=session,
                                    title = "Duplicated contrast name",
                                    text = msg,
                                    type = "error",
                                    btn_labels = "OK",
                                    ## btn_colors = "red",
                                    closeOnClickOutside = FALSE,
                                  )
                                } else if(nrow(df0)>1 && NCOL(df0)>=1) {
                                    df <- as.matrix(df0)
                                    
                                    matname <- "contrasts.csv"
                                }
                          }
                            if(!is.null(matname)) {
                                matlist[[matname]] <- df
                            }
                        }
                    }            
                }
                
                if("counts.csv" %in% names(matlist)) {
                    ## Convert to gene names (need for biological effects)
                    dbg("[upload_files] converting probe names to symbols")
                    X0 <- matlist[['counts.csv']]
                    pp <- rownames(X0)
                    rownames(X0) <- probe2symbol(pp)
                    sel <- !(rownames(X0) %in% c(NA,'','NA'))
                    X0 <- X0[sel,]
                    xx <- tapply(1:nrow(X0), rownames(X0), function(i) colSums(X0[i,,drop=FALSE]))
                    X0 <- do.call(rbind, xx)
                    matlist[['counts.csv']] <- X0
                }
                
                ## put the matrices in the reactive values 'uploaded'        
                files.needed = c("counts.csv","samples.csv","contrasts.csv")
                if(length(matlist)>0) {
                    matlist = matlist[ which(names(matlist) %in% files.needed) ]                
                    for(i in 1:length(matlist)) {
                        colnames(matlist[[i]]) <- gsub("[\n\t ]","_",colnames(matlist[[i]]))
                        rownames(matlist[[i]]) <- gsub("[\n\t ]","_",rownames(matlist[[i]]))
                        if(names(matlist)[i] %in% c("counts.csv","contrasts.csv")) {
                            matlist[[i]] <- as.matrix(matlist[[i]])
                        } else {
                            matlist[[i]] <- type.convert(matlist[[i]],as.is=FALSE)
                        }
                        m1 <- names(matlist)[i]
                        message("[upload_files] updating matrix ",m1)                        
                        uploaded[[m1]] <- matlist[[i]]
                    }
                    uploaded$last_uploaded <- names(matlist)
                }
                
                message("[upload_files] done!\n")
            })


            ##------------------------------------------------------------------
            ## Observer for loading from local exampledata.zip file
            ##
            ## Reads in the data files from zip and puts in the
            ## reactive values object 'uploaded'. Then uploaded should
            ## trigger the computePGX module.
            ## ------------------------------------------------------------------
            shiny::observeEvent( input$load_example, {

              dbg("[UploadModule] loading example data")
              dbg("[UploadModule] input$load_example = ",input$load_example)              

              if(input$load_example) {
                zipfile = file.path(FILES,"exampledata.zip")
                readfromzip1 <- function(file) {
                  read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE,
                    row.names=1)
                }
                readfromzip2 <- function(file) {
                  ## allows for duplicated names
                  df0 <- read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE)
                  mat <- as.matrix(df0[,-1])
                  rownames(mat) <- as.character(df0[,1])
                  mat
                }
                uploaded$counts.csv <- readfromzip2("exampledata/counts.csv")
                uploaded$samples.csv <- readfromzip1("exampledata/samples.csv")
                uploaded$contrasts.csv <- readfromzip1("exampledata/contrasts.csv")
                
              } else {
                dbg("[UploadModule] resetting upload")
                uploaded$counts.csv <- NULL
                uploaded$samples.csv <- NULL
                uploaded$contrasts.csv <- NULL
              }
            })

            ##------------------------------------------------------------------
            ## Observer for loading CSV from local folder on
            ## host/server using URL. Reads the CSV files from folder
            ## and puts in the reactive values object 'uploaded'.
            ## ------------------------------------------------------------------

            if(ALLOW_URL_QUERYSTRING) {

                shiny::observeEvent(session$clientData$url_search, {
                    ##-------------------------------------------------------------
                    ## Parse URL query string
                    ##-------------------------------------------------------------
                    query <- parseQueryString(session$clientData$url_search)
                    if(length(query)>0) {
                        dbg("[UploadModule:parseQueryString] names.query =",names(query))
                        for(i in 1:length(query)) {
                            dbg("[UploadModule:parseQueryString]",names(query)[i],"=>",query[[i]])
                        }
                    } else {
                        dbg("[UploadModule:parseQueryString] no queryString!")
                    }
                    
                    if(!is.null(query[['csv']])) {
                        qdir <- query[['csv']]
                        dbg("[UploadModule:parseQueryString] *** parseQueryString ***")
                        dbg("[UploadModule:parseQueryString] qdir = ",qdir)                         

                        counts_file = file.path(qdir,"counts.csv")
                        samples_file = file.path(qdir,"samples.csv")
                        if(!file.exists(counts_file)) {
                            dbg("[SERVER:parseQueryString] ***ERROR*** missing counts.csv in dir = ",qdir)
                        }
                        if(!file.exists(samples_file) ) {
                            dbg("[SERVER:parseQueryString] ***ERROR*** missing samples.csv in dir = ",qdir)
                        }
                        if(!file.exists(counts_file) || !file.exists(samples_file)) {
                            return(NULL)
                        }                        

                        FUN.readfromdir <- function() {
                            dbg("[UploadModule:parseQueryString] *** loading CSV from dir = ",qdir,"***")
                            
                            readfromdir1 <- function(file) {
                                read.csv(file, check.names=FALSE, stringsAsFactors=FALSE,
                                         row.names=1)
                            }
                            readfromdir2 <- function(file) {
                                ## allows for duplicated names
                                df0 <- read.csv(file, check.names=FALSE, stringsAsFactors=FALSE)
                                mat <- as.matrix(df0[,-1])
                                rownames(mat) <- as.character(df0[,1])
                                mat
                            }

                            dbg("[UploadModule:parseQueryString] reading samples_csv = ",samples_file)
                            uploaded$samples.csv <- readfromdir1(samples_file)

                            dbg("[UploadModule:parseQueryString] reading samples_csv = ",samples_file)
                            uploaded$counts.csv  <- readfromdir2(counts_file)
                            uploaded$contrasts.csv <- NULL

                            meta_file = file.path(qdir,"meta.txt")
                            uploaded$meta <- NULL
                            if(file.exists(meta_file)) {
                                dbg("[UploadModule:parseQueryString] reading meta file = ",meta_file)
                                ##meta <- read.table(meta_file,sep='\t',header=TRUE,row.names=1)
                                meta <- read.table(meta_file,sep='',header=TRUE,row.names=1)
                                meta <- as.list(array(meta[,1],dimnames=list(rownames(meta))))
                                uploaded$meta <- meta
                            }
                        }
                        
                        shinyalert::shinyalert(
                                        title= "Load CSV data from folder?",
                                        text = paste0("folder = ",qdir),
                                        callbackR = FUN.readfromdir,
                                        confirmButtonText = "Load!",
                                        type = "info")
                        
                        dbg("[UploadModule:parseQueryString] dim(samples) = ",dim(uploaded$samples.csv))
                        dbg("[UploadModule:parseQueryString] dim(counts) = ",dim(uploaded$counts.csv))

                        ## focus on this tab
                        updateTabsetPanel(session, "tabs", selected = "Upload data")
                    }
                    
                    if(0 && !is.null(query[['pgx']])) {

                        qdir <- query[['pgx']]
                        dbg("[UploadModule:parseQueryString] pgx =>",qdir)

                        pgx_file <- "../data/example-data.pgx"
                        pgx_file <- query[['pgx']]
                        pgx_file <- paste0(sub("[.]pgx$","",pgx_file),".pgx")
                        dbg("[UploadModule:parseQueryString] pgx_file = ",pgx_file)                         

                        if(!file.exists(pgx_file)) {
                            dbg("[SERVER:parseQueryString] ***ERROR*** missing pgx_file",pgx_file)
                            return(NULL)
                        }                        

                        dbg("[UploadModule:parseQueryString] 1:")
                        
                        FUN.readPGX <- function() {
                            dbg("[UploadModule:parseQueryString] *** loading PGX file = ",pgx_file,"***")

                            load(pgx_file)  ## load NGS/PGX                    
                            uploaded$pgx <- ngs
                            remove(ngs)
                            
                            uploaded$meta <- NULL
                        }
                        
                        dbg("[UploadModule:parseQueryString] 2:")

                        shinyalert::shinyalert(
                                        title= "Load PGX data from folder?",
                                        text = paste0("folder = ",qdir),
                                        callbackR = FUN.readPGX,
                                        confirmButtonText = "Load!",
                                        type = "info")

                        dbg("[UploadModule:parseQueryString] 3:")
                        
                        ## focus on this tab
                        updateTabsetPanel(session, "tabs", selected = "Upload data")

                        dbg("[UploadModule:parseQueryString] 4:")
                        
                    }

                })
    

                shiny::observeEvent( input$load_example, {

                    if(input$load_example) {
                        zipfile = file.path(FILES,"exampledata.zip")
                        readfromzip1 <- function(file) {
                            read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE,
                                     row.names=1)
                        }
                        readfromzip2 <- function(file) {
                            ## allows for duplicated names
                            df0 <- read.csv(unz(zipfile, file), check.names=FALSE, stringsAsFactors=FALSE)
                            mat <- as.matrix(df0[,-1])
                            rownames(mat) <- as.character(df0[,1])
                            mat
                        }
                        uploaded$counts.csv <- readfromzip2("exampledata/counts.csv")
                        uploaded$samples.csv <- readfromzip1("exampledata/samples.csv")
                        uploaded$contrasts.csv <- readfromzip1("exampledata/contrasts.csv")
                    } else {
                        ## Remove files                        
                        uploaded$counts.csv <- NULL
                        uploaded$samples.csv <- NULL
                        uploaded$contrasts.csv <- NULL
                    }
                })

            }
            
            
            ##=====================================================================
            ##========================= SUBMODULES/SERVERS ========================
            ##=====================================================================            

            upload_ok <- shiny::reactive({
                dbg("[UploadModule] upload_ok reactive")
                check <- checkUpload()
                all(check[,"status"]=="OK")
                all(grepl("ERROR",check[,"status"])==FALSE)
            })

            ##=====================================================================
            ##===================== PLOTS AND TABLES ==============================
            ##=====================================================================            
                
            checkUpload <- shiny::reactive({        
                ##
                ##
                ##
                
                ## check dimensions
                status = rep("please upload",3)
                files.needed = c("counts.csv","samples.csv","contrasts.csv")        
                names(status) = files.needed
                files.nrow = rep(NA,3)
                files.ncol = rep(NA,3)
                
                for(i in 1:3) {
                    fn = files.needed[i]
                    upfile <- uploaded[[fn]]
                    if(fn %in% names(uploaded) && !is.null(upfile)) {
                        status[i] = "OK"
                        files.nrow[i] = nrow(upfile)
                        files.ncol[i] = ncol(upfile)
                    }
                }

                has.pgx <- ("pgx" %in% names(uploaded))
                if(has.pgx) has.pgx <- has.pgx && !is.null(uploaded$pgx)
                if(has.pgx==TRUE) {
                    
                    ## Nothing to check. Always OK.            
                    
                } else if(!has.pgx) {
                    
                    ## check rownames of samples.csv
                    if(status["samples.csv"]=="OK" && status["counts.csv"]=="OK") {
                        
                        samples1 <- uploaded$samples.csv
                        counts1  <- uploaded$counts.csv
                        a1 <- mean(rownames(samples1) %in% colnames(counts1))
                        a2 <- mean(samples1[,1] %in% colnames(counts1))
                        
                        if(a2 > a1 && NCOL(samples1)>1 ) {
                            message("[UploadModuleServer] getting sample names from first column\n")
                            rownames(samples1) <- samples1[,1]
                            uploaded$samples.csv <- samples1[,-1,drop=FALSE]
                        }                        
                    }
                    
                    ## check files: matching dimensions
                    if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        nsamples   <- max( ncol(uploaded$counts.csv), nrow(uploaded$samples.csv) )
                        ok.samples <- intersect(rownames(uploaded$samples.csv),
                                                colnames(uploaded$counts.csv))
                        n.ok <- length(ok.samples)
                        message("[UploadModule::checkTables] n.ok = ",n.ok)
                        if(n.ok > 0 && n.ok < nsamples) {
                            ## status["counts.csv"]  = "WARNING: some samples with missing annotation)"
                        }

                        if(n.ok > 0) {
                            message("[UploadModule::checkTables] conforming samples/counts...")
                            uploaded$samples.csv <- uploaded$samples.csv[ok.samples,,drop=FALSE]
                            uploaded$counts.csv  <- uploaded$counts.csv[,ok.samples,drop=FALSE]
                        }

                        if(n.ok == 0) {
                            status["counts.csv"]  = "ERROR: colnames do not match (with samples)"
                            status["samples.csv"] = "ERROR: rownames do not match (with counts)"
                        }

                    }
                                                     
                    if(status["contrasts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        samples1   <- uploaded$samples.csv
                        contrasts1 <- uploaded$contrasts.csv                        
                        contrasts.new <- pgx.convertContrastsExplicit(contrasts1, samples1)
                        ok.contrast <- all(rownames(samples1) == rownames(contrasts.new))
                        ok.contrast                        
                        if(ok.contrast && NCOL(contrasts.new)>0) {
                            uploaded$contrasts.csv <- contrasts.new
                            status["contrasts.csv"] <- "OK"
                        } else {
                            uploaded$contrasts.csv <- NULL
                            status["contrasts.csv"] <- "ERROR: dimension mismatch"
                        }
                    }

                    MAXSAMPLES   = 25
                    MAXCONTRASTS = 5
                    MAXSAMPLES   = as.integer(limits["samples"])
                    MAXCONTRASTS = as.integer(limits["comparisons"])

                    ## check files: maximum contrasts allowed
                    if(status["contrasts.csv"]=="OK") {
                        if( ncol(uploaded$contrasts.csv) > MAXCONTRASTS ) {
                            status["contrasts.csv"] = paste("ERROR: max",MAXCONTRASTS,"contrasts allowed")
                        }
                    }
                    
                    ## check files: maximum samples allowed
                    if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
                        if( ncol(uploaded$counts.csv) > MAXSAMPLES ) {
                            status["counts.csv"]  = paste("ERROR: max",MAXSAMPLES," samples allowed")
                        }
                        if( nrow(uploaded$samples.csv) > MAXSAMPLES ) {
                            status["samples.csv"] = paste("ERROR: max",MAXSAMPLES,"samples allowed")
                        }
                    }
                    
                    ## check samples.csv: must have group column defined
                    if(status["samples.csv"]=="OK" && status["contrasts.csv"]=="OK") {
                        samples1   = uploaded$samples.csv
                        contrasts1 = uploaded$contrasts.csv
                        if(!all(rownames(contrasts1) %in% rownames(samples1))) {
                            status["contrasts.csv"] = "ERROR: contrasts do not match samples"
                        }
                    }                    
                    
                } ## end-if-from-pgx
                
                e1 <- grepl("ERROR",status["samples.csv"])
                e2 <- grepl("ERROR",status["contrasts.csv"]) 
                e3 <- grepl("ERROR",status["counts.csv"])
                s1 <- "samples.csv" %in% uploaded$last_uploaded
                s2 <- "contrasts.csv" %in% uploaded$last_uploaded
                s3 <- "counts.csv" %in% uploaded$last_uploaded
                
                if( e1 || e2 || e3 ) {
                    message("[checkTables] ERROR in samples table : e1 = ",e1)
                    message("[checkTables] ERROR in contrasts table : e2 = ",e2)
                    message("[checkTables] ERROR in counts table : e2 = ",e3)

                    if(e1 && !s1) {
                        uploaded$samples.csv <- NULL
                        status["samples.csv"] = "please upload"
                    }
                    if(e2 && !s2) {
                        uploaded$contrasts.csv <- NULL
                        status["contrasts.csv"] = "please upload"
                    }
                    if(e3 && !s3) {
                        uploaded$counts.csv <- NULL
                        status["counts.csv"] = "please upload"
                    }
                }


                if( !is.null(uploaded$contrasts.csv) &&
                    (is.null(uploaded$counts.csv) ||
                     is.null(uploaded$samples.csv)) )
                {
                    uploaded$contrasts.csv <- NULL
                    status["contrasts.csv"] = "please upload"
                }
                
                
                ## check files
                description = c(
                    "Count/expression file with gene on rows, samples as columns",
                    "Samples file with samples on rows, phenotypes as columns",
                    ## "Gene information file with genes on rows, gene info as columns.",
                    "Contrast file with conditions on rows, contrasts as columns"        
                )
                df <- data.frame(
                    filename = files.needed,
                    description = description,
                    nrow = files.nrow,
                    ncol = files.ncol,
                    status = status
                )
                rownames(df) <- files.needed
                
                ## deselect
                ## DT::selectRows(proxy = DT::dataTableProxy("pgxtable"), selected=NULL)
                return(df)    
            })
            
            output$checkTablesOutput <- DT::renderDataTable({
                ## Render the upload status table
                ##
                if(!input$advanced_mode) return(NULL)
                df <- checkTables()
                dt <- DT::datatable(
                    df,
                    rownames=FALSE,
                    selection = 'none',
                    class="compact cell-border",
                    options = list(
                        dom = 't'
                    )                                
                ) %>%
                    DT::formatStyle(0, target='row', fontSize='12px', lineHeight='100%')                
            })               


          output$countstable <- DT::renderDataTable({
            df <- uploaded$counts.csv
            if(!is.null(df) && nrow(df)) {
              df <- round(df, digits=2)
              df <- data.frame("(rownames)"=rownames(df),df,check.names=FALSE)
            } else {
              mx1 <- matrix(HTML(" "), 10, 10)
              genes <- paste0("gene_",1:nrow(mx1))
              colnames(mx1) <- paste0("sample_",1:10)
              df <- data.frame("(rownames)"=genes, mx1, check.names=FALSE)
            }
            dt <- DT::datatable(
              df,
              rownames = FALSE,
              selection = 'none',
              class="compact cell-border",
              options = list(
                dom = 't'
              )                                
            ) %>%
              DT::formatStyle(0, target='row', fontSize='12px', lineHeight='80%')                
          })

          output$sampletable <- DT::renderDataTable({
            df <- uploaded$samples.csv
            if(!is.null(df) && nrow(df)) {
              df <- data.frame("(rownames)"=rownames(df),df,check.names=FALSE)
            } else {
              mx1 <- matrix(HTML(" "), 10, 8)
              colnames(mx1) <- paste0("phenotype_",1:8)
              samples <- paste0("sample_",1:nrow(mx1))
              df <- data.frame("(rownames)"=samples, mx1, check.names=FALSE)
            }

            dt <- DT::datatable(
              df,
              rownames = FALSE,
              selection = 'none',
              class="compact cell-border",
              options = list(
                dom = 't'
              )                                
            ) %>%
              DT::formatStyle(0, target='row', fontSize='12px', lineHeight='80%')                
          })

          output$contrasttable <- DT::renderDataTable({
            df <- uploaded$contrasts.csv
            if(!is.null(df) && nrow(df)) {
              df <- data.frame("(rownames)"=rownames(df),df,check.names=FALSE)            
            } else {
              mx1 <- matrix(HTML(" "), 10, 8)
              colnames(mx1) <- paste0(LETTERS[1:8],"_vs_ctrl")              
              samples <- paste0("sample_",1:nrow(mx1))              
              df <- data.frame("(rownames)"=samples, mx1, check.names=FALSE)
            }
            dt <- DT::datatable(
              df,
              rownames = FALSE,
              selection = 'none',
              class="compact cell-border",
              options = list(
                dom = 't'
              )                                
            ) %>%
              DT::formatStyle(0, target='row', fontSize='12px', lineHeight='80%')                
          })

          
          ##========================================================================
          ## return results as reactive object
          ##========================================================================            
          ## return(shiny::reactive(uploaded$pgx))  ## pointing to reactive results object
                    
          return(list(
            counts.csv    = reactive(uploaded$counts.csv),
            samples.csv   = reactive(uploaded$samples.csv),
            contrasts.csv = reactive(uploaded$contrasts.csv),
##          pgx           = reactive(uploaded$pgx),
            meta          = reactive(uploaded$meta)
          ))
            
        })  ## end moduleServer

} ## end UploadModuleServer
