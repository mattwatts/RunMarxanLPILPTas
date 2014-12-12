# Author: Matt Watts
# Date: 12 Dec 2014
# Purpose: RunMarxanLPILPTas web app ui.R

require(shiny)

shinyUI(pageWithSidebar(

    headerPanel("Tas Activity: compare Marxan, LP, and ILP algorithms"),

    sidebarPanel(
        actionButton("mrun","Run"), 
        br(),
        br(),
        textOutput("textfeedback"),
        br(),
        br(),
        selectInput("feature", "Choose a species to edit:",
                    choices = c("10","11","12","13","14","15","16","17","18","19",
                                "20","21","22","23","24","25","26")),
        numericInput("target", "Target:",0.1,min=0),
        numericInput("spf", "SPF:",0.1,min=0),
        actionButton("savetargetspf","Save Target and SPF"),
        conditionalPanel(condition = "input.tabs == 'Map'",
            br(),
            br(),
            radioButtons("map", "Map to display:",
                         list("ILP" = "ilpmap",
                              "LP" = "lpmap",     	                  
                              "Best solution" = "bestmap",
                              "Solution M" = "runMmap",
                              "Selection frequency" = "ssolnNmap"))
        ),
        conditionalPanel(condition = "input.tabs == 'Map' & input.map == 'runMmap'",
            br(),
            br(),
            sliderInput("m", "Solution M:",
                        value = 1,
                        min = 1,
                        max = iNUMREPS, step = 1)

        ),
        conditionalPanel(condition = "input.tabs == 'Table'",
            br(),
            br(),
            radioButtons("table", "Table to display:",
                         list("Summary" = "sumtable",
                              "Conservation Features" = "spec"))
        ),
        conditionalPanel(condition = "input.prop == -1",
                         numericInput("refreshinput", "Refresh Input", 0))
    ),

    mainPanel(
        tabsetPanel(id="tabs",
            tabPanel("Map", plotOutput('marxanmap')),
            tabPanel("Table", tableOutput('marxantable')),
            tabPanel("Cluster", plotOutput("plot2ds"),
                                plotOutput("plotdendogram"))
                  )
             )
))
