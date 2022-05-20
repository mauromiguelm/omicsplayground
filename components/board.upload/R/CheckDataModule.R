##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##================================================================================
## Upload data Module
##================================================================================

# From https://stackoverflow.com/questions/36995142/get-the-size-of-the-window-in-shiny
# https://stackoverflow.com/a/37060206


CheckDataUI <- function(id) {
  ns <- shiny::NS(id)
  fluidRow(
    column(4,
      shiny::plotOutput(ns("countStats")) %>% shinycssloaders::withSpinner()
    ),
    column(4,
      shiny::plotOutput(ns("phenoStats")) %>% shinycssloaders::withSpinner()
    ),
    column(4,
      plotWidget(ns("pcaplot")) %>% shinycssloaders::withSpinner()
    ) #,
#    column(3,
#      shiny::plotOutput(ns("contrastStats")) %>% shinycssloaders::withSpinner()
#    )
  )
}


CheckDataModule <- function(id, 
                            countsRT = reactive(matrix(0)), 
                            samplesRT = reactive(matrix(0)), 
                            contrastsRT = reactive(matrix(0))
                            )
{
  shiny::moduleServer(
    id,
    function(input, output, session) {

      ns <- session$ns
      
      dbg("[CheckDataModuleServer] called!")            
      
      files <- reactiveValues(
        counts.csv    = NULL,
        samples.csv   = NULL,
        contrasts.csv = NULL,
        pgx = NULL
      )
        
      observeEvent( countsRT(), {
        dbg("[CheckDataModule] countsRT reacted! dim.countsRT = ", dim(countsRT()) )                    
        if(is.null(countsRT) || is.null(countsRT())) {
          files$counts.csv  <- NULL          
        } else {
          files$counts.csv <- countsRT()
        }
      })
      
      observeEvent( samplesRT(), {
        dbg("[CheckDataModule] samplesRT reacted! dim.samplesRT = ", dim(samplesRT()) )
        if(is.null(samplesRT) || is.null(samplesRT())) {
          files$samples.csv  <- NULL
        } else {
          files$samples.csv <- samplesRT()
        }
      })
      
      observeEvent( contrastsRT(), {
        dbg("[CheckDataModule] samplesRT reacted! dim.samplesRT = ", dim(samplesRT()) )        
        if(is.null(contrastsRT) || is.null(contrastsRT())) {
          files$contrasts.csv  <- NULL
        } else {
          contrasts.new <- pgx.convertContrastsExplicit(contrastsRT(), samplesRT())          
          files$contrasts.csv  <- contrasts.new
        }
      })
      
      ##------------------------------------------------------------------------------
      ##------------------------------------------------------------------------------
      ##------------------------------------------------------------------------------
      
      checkTables <- shiny::reactive({        
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
          upfile <- files[[fn]]
          if(fn %in% names(uploaded) && !is.null(upfile)) {
            status[i] = "OK"
            files.nrow[i] = nrow(upfile)
            files.ncol[i] = ncol(upfile)
          }
        }
        
        has.pgx <- ("pgx" %in% names(uploaded))
        if(has.pgx) has.pgx <- has.pgx && !is.null(files[["pgx"]])
        if(has.pgx==TRUE) {
          
          ## Nothing to check. Always OK.            
          
        } else if(!has.pgx) {
          
          ## check rownames of samples.csv
          if(status["samples.csv"]=="OK" && status["counts.csv"]=="OK") {
            
            samples1 <- files$samples.csv
            counts1  <- files$counts.csv
            a1 <- mean(rownames(samples1) %in% colnames(counts1))
            a2 <- mean(samples1[,1] %in% colnames(counts1))
            
            if(a2 > a1 && NCOL(samples1)>1 ) {
              message("[CheckDataModuleServer] getting sample names from first column\n")
              rownames(samples1) <- samples1[,1]
              files$samples.csv <- samples1[,-1,drop=FALSE]
            }                        
          }
          
          ## check files: matching dimensions
          if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
            nsamples   <- max( ncol(files$counts.csv), nrow(files$samples.csv) )
            ok.samples <- intersect(rownames(files$samples.csv),
              colnames(files$counts.csv))
            n.ok <- length(ok.samples)
            message("[CheckDataModule::checkTables] n.ok = ",n.ok)
            if(n.ok > 0 && n.ok < nsamples) {
              ## status["counts.csv"]  = "WARNING: some samples with missing annotation)"
            }
            
            if(n.ok > 0) {
              message("[CheckDataModule::checkTables] conforming samples/counts...")
              files$samples.csv <- files$samples.csv[ok.samples,,drop=FALSE]
              files$counts.csv  <- files$counts.csv[,ok.samples,drop=FALSE]
            }
            
            if(n.ok == 0) {
              status["counts.csv"]  = "ERROR: colnames do not match (with samples)"
              status["samples.csv"] = "ERROR: rownames do not match (with counts)"
            }
                        
          }
          
          if(status["contrasts.csv"]=="OK" && status["samples.csv"]=="OK") {
            samples1   <- files$samples.csv
            contrasts1 <- files$contrasts.csv
            group.col <- grep("group", tolower(colnames(samples1)))
            old1 = (length(group.col)>0 &&
                      nrow(contrasts1) < nrow(samples1) &&
                      all(rownames(contrasts1) %in% samples1[,group.col])
            )
            old2 = all(rownames(contrasts1)==rownames(samples1)) &&
              all(unique(as.vector(contrasts1)) %in% c(-1,0,1,NA))
            
            old.style <- (old1 || old2)
            if(old.style && old1) {
              
              message("[CheckDataModule] WARNING: converting old1 style contrast to new format")
              new.contrasts <- samples1[,0]
              if(NCOL(contrasts1)>0) {
                new.contrasts <- contrastAsLabels(contrasts1)
                grp = as.character(samples1[,group.col])
                new.contrasts <- new.contrasts[grp,,drop=FALSE]
                rownames(new.contrasts) <- rownames(samples1)
              }
              
              contrasts1 <- new.contrasts
            }
            if(old.style && old2 ) {
              message("[CheckDataModule] WARNING: converting old2 style contrast to new format")
              new.contrasts <- samples1[,0]
              if(NCOL(contrasts1)>0) {
                new.contrasts <- contrastAsLabels(contrasts1)
                rownames(new.contrasts) <- rownames(samples1)
              }
              contrasts1 <- new.contrasts
            }
                        
            ok.contrast <- length(intersect(rownames(samples1),rownames(contrasts1)))>0
            if(ok.contrast && NCOL(contrasts1)>0) {
              ## always clean up
              contrasts1 <- apply(contrasts1,2,as.character)
              rownames(contrasts1) <- rownames(samples1)
              for(i in 1:ncol(contrasts1)) {
                isz = (contrasts1[,i] %in% c(NA,"NA","NA ",""," ","  ","   "," NA"))
                if(length(isz)) contrasts1[isz,i] <- NA
              }
              files$contrasts.csv <- contrasts1
              status$ontrasts.csv <- "OK"
            } else {
              files$contrasts.csv <- NULL
              status$contrasts.csv <- "ERROR: dimension mismatch"
            }
          }
          
          MAXSAMPLES   = 25
          MAXCONTRASTS = 5
          ##MAXSAMPLES   = as.integer(limits["samples"])
          ##MAXCONTRASTS = as.integer(limits["comparisons"])
          
          ## check files: maximum contrasts allowed
          if(status["contrasts.csv"]=="OK") {
            if( ncol(files$contrasts.csv) > MAXCONTRASTS ) {
              status["contrasts.csv"] = paste("ERROR: max",MAXCONTRASTS,"contrasts allowed")
            }
          }
          
          ## check files: maximum samples allowed
          if(status["counts.csv"]=="OK" && status["samples.csv"]=="OK") {
            if( ncol(files$counts.csv) > MAXSAMPLES ) {
              status["counts.csv"]  = paste("ERROR: max",MAXSAMPLES," samples allowed")
            }
            if( nrow(files$samples.csv) > MAXSAMPLES ) {
              status["samples.csv"] = paste("ERROR: max",MAXSAMPLES,"samples allowed")
            }
          }
          
          ## check samples.csv: must have group column defined
          if(status["samples.csv"]=="OK" && status["contrasts.csv"]=="OK") {
            samples1   = files$samples.csv
            contrasts1 = files$contrasts.csv
            if(!all(rownames(contrasts1) %in% rownames(samples1))) {
              status["contrasts.csv"] = "ERROR: contrasts do not match samples"
            }
          }                    
          
        } ## end-if-from-pgx
        
        e1 <- grepl("ERROR",status["samples.csv"]) || is.null(files$samples.csv)
        e2 <- grepl("ERROR",status["contrasts.csv"]) || is.null(files$contrasts.csv)
        e3 <- grepl("ERROR",status["counts.csv"]) || is.null(files$counts.csv)
        
        if( e1 || e2 || e3 ) {
          message("[checkTables] ERROR in samples table : e1 = ",e1)
          message("[checkTables] ERROR in contrasts table : e2 = ",e2)
          message("[checkTables] ERROR in counts table : e2 = ",e3)
          
          if(e1) {
            files$samples.csv <- NULL
            status["samples.csv"] = "please upload"
          }
          if(e2) {
            files$contrasts.csv <- NULL
            status["contrasts.csv"] = "please upload"
          }
          if(e3) {
            files$counts.csv <- NULL
            status["counts.csv"] = "please upload"
          }
        }
                
        if( !is.null(files$contrasts.csv) &&
              (is.null(files$counts.csv) ||
                 is.null(files$samples.csv)) )
        {
          files$contrasts.csv <- NULL
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
      
      ##=====================================================================
      ##===================== PLOTS AND TABLES ==============================
      ##=====================================================================            

      ##--------------------------------------------
      ##--------------------------------------------

      output$countStats <- shiny::renderPlot({

        dbg("[CheckDataModule:countStats] renderPlot called")
        ##req(files$counts.csv)                
        has.counts <- !is.null(files$counts.csv) && NCOL(files$counts.csv)>0
        
        if(!has.counts) {
          status.ok <- "MISSING FILE"
          frame()
          status.ds <- "(Required) Please upload a 'counts.csv' file containing your counts/expression values with genes as rows and samples as columns."          
          ##status.ds <- check["counts.csv","description"]
          msg <- paste(toupper(status.ok),"\n\n",status.ds)
          graphics::text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
          graphics::box(lty=2, col="grey60")
          return(NULL)
        }
        
        counts <- files$counts.csv
        xx <- log2(1 + counts)
        if(nrow(xx)>1000) xx <- xx[sample(1:nrow(xx),1000),,drop=FALSE]
        ##dc <- reshape::melt(xx)
        suppressWarnings( dc <- data.table::melt(xx) )
        dc$value[dc$value==0] <- NA
        tt2 <- paste(nrow(counts),"genes x",ncol(counts),"samples")
        ggplot2::ggplot(dc, ggplot2::aes(x=value, color=Var2)) +
          ggplot2::geom_density() + ggplot2::xlab("log2(1+counts)") +
          ggplot2::theme( legend.position = "none") +
          ggplot2::ggtitle("Counts", subtitle=tt2)
      })

      ##--------------------------------------------
      ##--------------------------------------------

      output$phenoStats <- shiny::renderPlot({

        dbg("[phenoStats] renderPlot called \n")
        ##req(files$samples.csv)                
        has.samples <- !is.null(files$samples.csv) && NCOL(files$samples.csv)>0
        
        if(!has.samples) {
          ##check <- checkTables()
          ##status.ok <- check["samples.csv","status"]                
          ##if(status.ok!="OK") {
          status.ok = "MISSING FILE"
          frame()
          status.ds <- "(Required) Please upload a 'samples.csv' file specifying your sample information with samples as rows and phenotype conditions as columns."
          ##status.ds <- check["samples.csv","description"]
          msg <- paste(toupper(status.ok),"\n\n",status.ds)
          graphics::text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
          graphics::box(lty=2, col="grey60")
          return(NULL)
        }
        
        pheno <- files$samples.csv             
        px <- head(colnames(pheno),20)  ## show maximum??

        df <- type.convert(pheno[,px,drop=FALSE],as.is=FALSE)
        vt <- df %>% inspectdf::inspect_types()
        vt
        
        ## discretized continuous variable into 10 bins
        ii <- unlist(vt$col_name[c("numeric","integer")])
        ii
        if(!is.null(ii) && length(ii)) {
          cat("[CheckDataModule::phenoStats] discretizing variables:",ii,"\n")
          df[,ii] <- apply(df[,ii,drop=FALSE], 2, function(x) {
            if(any(is.infinite(x))) x[which(is.infinite(x))] <- NA
            cut(x, breaks=10)
          })
        }
        
        p1 <- df %>% inspectdf::inspect_cat() %>% inspectdf::show_plot()
        tt2 <- paste(nrow(pheno),"samples x",ncol(pheno),"phenotypes")
        ## tt2 <- paste(ncol(pheno),"phenotypes")
        p1 <- p1 + ggplot2::ggtitle("Phenotypes", subtitle=tt2) +
          ggplot2::theme(
            ##axis.text.x = ggplot2::element_text(size=8, vjust=+5),
            axis.text.y = ggplot2::element_text(
              size = 12,
              margin = ggplot2::margin(0,0,0,25),
              hjust = 1)
          )

        dbg("[CheckDataModule::phenoStats] done!")
        
        p1
      })
      
      ##--------------------------------------------
      ##--------------------------------------------

      output$contrastStats <- shiny::renderPlot({
        
        contrasts <- files$contrasts.csv
        has.contrasts <- !is.null(contrasts) && NCOL(contrasts)>0

        status.ok <- "OK"
        ##check <- checkTables()
        ##status.ok <- check["contrasts.csv","status"]
        
        if( status.ok!="OK" || !has.contrasts) {
          frame()
          status.ds <- "(Optional) You can optionally upload a 'contrasts.csv' file specifying your contrasts. Or you can create your contrasts in the next step."
          ##status.ds <- check["contrasts.csv","description"]
          msg <- paste(toupper(status.ok),"\n\n",status.ds)
          ##text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
          graphics::text(0.5,0.5,paste(strwrap(msg,30),collapse="\n"),col="grey25")
          graphics::box(lty=2, col="grey60")
          return(NULL)
        }

        ##contrasts <- sign(contrasts)
        ##df <- contrastAsLabels(contrasts)
        df <- contrasts
        px <- head(colnames(df),20)  ## maximum to show??
        df <- data.frame(df[,px,drop=FALSE],check.names=FALSE)
        tt2 <- paste(nrow(contrasts),"samples x",ncol(contrasts),"contrasts")
        p1 <- df %>% inspectdf::inspect_cat() %>% inspectdf::show_plot()                    
        p1 <- p1 + ggplot2::ggtitle("Contrasts", subtitle=tt2) +
          ggplot2::theme(
            ##axis.text.x = ggplot2::element_text(size=8, vjust=+5),
            axis.text.y = ggplot2::element_text(size = 12,
              margin = ggplot2::margin(0,0,0,25),
              hjust = 1)
          )
        
        p1
      })
      
      ##------------------------------------------------------------------------------
      ##------------------------------------------------------------------------------
      ##------------------------------------------------------------------------------



      getClust <- reactive({
        
        counts <- countsRT()
        shiny::req(counts)
        
        method  <- input$pcaplot.method
        
        X <- log2(1 + counts)
        clust <- pgx.clusterMatrix(
          X,
          dims = 2,
          method = method,
          ntop = 400,
          npca = 20,
          find.clusters = FALSE
        )
        clust
      })
      
      pcaplot.RENDER <- shiny::reactive({

        dbg("[CheckDataModule] pcaplot.RENDER : reacted")
        ##ngs <- inputData()
        ##X <- ngs$X
        pheno  <- samplesRT()
        counts <- countsRT()
        shiny::req(pheno)
        shiny::req(counts)

        dbg("[CheckDataModule] pcaplot.RENDER : 1")
        
        clust <- getClust()
        names(clust)

        dbg("[MakeContrastServer] pcaplot.RENDER : 2")
        colorby <- input$pcaplot.colorby
        if(!is.null(colorby) && colorby %in% colnames(samplesRT())) {
          cond <- samplesRT()[,colorby]
        } else {
          cond <- NULL
        }
        ##par(mar=c(4,1,1,1))
        pgx.scatterPlotXY(
          clust$pos2d, var=cond, plotlib="plotly",
          legend = FALSE ##, labels=TRUE
        )
        
      })
      

      observeEvent(samplesRT(), {
        dbg("[CheckDataModule] updating colorby : dim.samples = ",dim(samplesRT()))
        if(NCOL(samplesRT())>0) {
          vars <- colnames(samplesRT())
          dbg("[CheckDataModule] updating colorby: ",vars)          
          freezeReactiveValue(input, ns("pcaplot.colorby"))          
          updateSelectInput(session, "pcaplot.colorby", choices=vars )
        }
      })

      pcaplot.opts = shiny::tagList(
        withTooltip( shiny::selectInput( ns("pcaplot.method"),
          "Method:", choices = c("pca","tsne","umap"),
          width = '100%'),"Choose clustering method.",
          placement="right", options = list(container = "body")),
        withTooltip( shiny::selectInput( ns("pcaplot.colorby"),
          "Color by:", choices = NULL,
          width = '100%'),"Choose color variable.",
          placement="right", options = list(container = "body"))        
      )
      
      shiny::callModule(
        plotModule, 
        id = "pcaplot",
        func = pcaplot.RENDER, ## ns=ns,
        plotlib = "plotly", 
        options = pcaplot.opts,
        height = c(320,700), width=c("auto",800),
        pdf.width=8, pdf.height=8,
        title="Samples"
      )

      ##outputOptions(output, "pcaplot.colorby", suspendWhenHidden = FALSE)
      ##outputOptions(output, "pcaplot-widget", suspendWhenHidden = FALSE)            


      ##========================================================================
      ## return results as reactive object
      ##========================================================================            
      ## return(shiny::reactive(files$pgx))  ## pointing to reactive results object
      ##return(computed_pgx)

      ok <- reactive({
        !is.null(files$counts.csv) && !is.null(files$samples.csv)
      })
      return(ok)
      
    })  ## end moduleServer

} ## end CheckDataModuleServer
