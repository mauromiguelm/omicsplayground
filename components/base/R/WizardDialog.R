##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

WizardDialogUI <- function(id, pages, showButtons=TRUE, doneButton = NULL) {

  require(shiny)

  stopifnot(is.list(pages))
  n <- length(pages)
  ns <- NS(id)  
  nextPage <- function(i) {
    actionButton( ns(paste0("go_", i, "_", i + 1)), "next \u203A",
      class = "btn btn-primary action-button", style="margin: 0px 0px;")
  }
  prevPage <- function(i) {
    actionButton( ns(paste0("go_", i, "_", i - 1)), "\u2039 back",
      class = "btn btn-seconday action-button", style="margin: 0px 3px;")
  }

  arrow <- "clip-path: polygon(0% 0%, 93% 0%, 93.5% 0.5%, 94% 2%, 99.95% 48%, 99.95% 52%, 94% 98%, 93.5% 99.5%, 93% 100%, 0% 100%);"    
  arrow.style <- paste0(arrow,";position:relative;height:44px;width:",60/n,"vw;padding:10px 25px;border:0px;background-color:#eee;border-radius:8px;font-size:1.3em;")
  
  wrapPage <- function(i, page, button_left = NULL, button_right = NULL) {
    title <- names(pages)[i]
    id <- ns(paste0("page",i))
    nav.buttons <- NULL
    if(showButtons) {
      nav.buttons <- tagList(
        hr(),
        fillRow(
          flex = c(NA,NA,1),
          button_left,
          button_right,
          br()
        )
      )
    }

    tabPanel(
      title = div(title,id=id,style=arrow.style,class="wiz-tabs"),
      value = id,
      flex = c(NA,NA,NA,NA,NA,1),
      hr(),      
      fluidRow(
        column(12, page)
      ),
      nav.buttons,
      br()
    )
  }

  wrapped <- vector("list", n)
  for (i in seq_along(pages)) {
    # First page only has next; last page only prev + done
    lhs <- if (i > 1) prevPage(i)
    rhs <- if (i < n) nextPage(i) else doneButton
    wrapped[[i]] <- wrapPage(i, pages[[i]], lhs, rhs)
  }
  
  # Create tabsetPanel
  # https://github.com/rstudio/shiny/issues/2927
  wrapped$id <- ns("wizpanel")
##  wrapped$type <- "hidden"
  
  inline.css <- tags$head(tags$style(HTML("
      .wizpanel { margin-top:20px;}
      .wizpanel>.tabbable>.nav-tabs { border:0px;}
      .wizpanel>.tabbable>.nav-tabs>li>a { border:0px;}
      .wizpanel>.tabbable>.nav-tabs>li.active>a { border:0px;}
      .wizpanel>.tabbable>.nav>li>a:hover { border: 0px; background-color: white;}
      .wizpanel .active>a .wiz-tabs {color: white !important; font-weight:900;background-color: crimson !important;}
    ")))
  
  tagList(
    shinyjs::useShinyjs(),  # Include shinyjs
    inline.css,
    div( do.call("tabsetPanel", wrapped), class="wizpanel")
  )
}


WizardDialogServer <- function(id, n) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    changePage <- function(from, to) {
      observeEvent(input[[(paste0("go_", from, "_", to))]], {
        ##message("reacted on button: ",ns(paste0("go_", from, "_", to)))
        page.from <- paste0("page", from)
        page.to <- paste0("page", to)        
        updateTabsetPanel(session, "wizpanel", selected = ns(page.to))
      })  
    }

    ids <- seq_len(n)
    lapply(ids[-1], function(i) changePage(i, i-1))
    lapply(ids[-n], function(i) changePage(i, i+1))
  })
}

if(1 && interactive()) {

  ##
  ## Example 
  ##
  
  page1 <- fluidRow(
    column(6, textInput("name", "What's your name?")),
    column(6,
      tabsetPanel(
        tabPanel("tab1","Hello this Tab number 1"),
        tabPanel("tab2","Hello this Tab number 2")        
      )
    )
  )
  page2 <- tagList(
    numericInput("age", "How old are you?", 20)
  )
  page3 <- tagList(
    tags$b("Is this data correct?"),
    verbatimTextOutput("info")
  )
  pages = list("Enter name"=page1, "Enter age"=page2, "Check"=page3)
  
  ui <- fluidPage(
    WizardDialogUI(
      id = "wizdemo", 
      pages = pages,
##      showButtons = FALSE,
      doneButton = actionButton("done", "Submit", class = "btn btn-primary action-button")
    )
  )
  
  server <- function(input, output, session) {
    WizardDialogServer("wizdemo", 3)
    
    observeEvent(input$done, showModal(
      modalDialog("Thank you!", footer = NULL)
    ))
    
    output$info <- renderText(paste0(
      "Age: ", input$age, "\n",
      "Name: ", input$name, "\n"
    ))
  }
  
  shinyApp(ui, server)

}
