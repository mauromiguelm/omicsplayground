##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadBoardInputs <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  bigdash::tabSettings(
    ## shiny::actionLink(ns("module_info"), "Tutorial", icon = shiny::icon("youtube"))
      shinyWidgets::prettySwitch(ns("advanced_mode"),"batch correction (beta)"),
      shinyWidgets::prettySwitch(ns("upload-show_tables"),"show tables")
  )
}

UploadBoardUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace
 
  shiny::tagList(
    ## header ---------------------------
    shiny::div(
      id="navheader-current-section",
      HTML("Upload data &nbsp;"), 
      shiny::actionLink(
        ns("module_info"), "",
        icon = shiny::icon("info-circle"),
        style = "color: #ccc; padding: 0px 10px;"
      )
    ),        
    WizardDialogUI(
      id = ns("wiz"),
      pages = list(
        "Upload data"    = UploadModuleUI(ns("upload")),
        ## "Normalize" = shiny::uiOutput(ns("normalize_UI")),
        "BatchCorrect" = BatchCorrectUI(ns("batchcorrect")),        
        "Create comparisons" =  MakeContrastUI(ns("makecontrast")),
        "Compute!"   = ComputePgxUI(ns("compute"))
      ),
      showButtons = FALSE,
      doneButton = NULL
    )
  )
  
}
