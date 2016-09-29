library(shiny)

shinyUI(fluidPage(
    titlePanel("Normalization of FCS Data"),
    hr(),

    fluidRow(
        column(3,
               h4('Load Data:'),
               wellPanel(
                   fileInput(inputId = 'fcsFiles',
                             label = "FCS files",
                             multiple = TRUE,
                             accept = c('text/fcs', '.fcs')),
                   fileInput(inputId = 'compensationFile',
                             label = "Compensation File(csv)",
                             multiple = FALSE,
                             accept = c(
                                 "text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
                   selectInput('transMethod', 'Transformation Method:',
                               choices = c("autoLgcl", "cytofAsinh", "logicle",
                                           "arcsinh", "none"),
                               selected = "autoLgcl",
                               width = "100%"),
                   actionButton("loadButton", "Load Data", icon = icon("hand-o-right"))
               ),

               hr(),
               uiOutput("normDataCheck"),

               hr(),
               h4("Sample Filter:"),
               wellPanel(uiOutput("sampleFilter")),

               hr(),
               actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ",
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
               tabsetPanel(type = "pills",
                           tabPanel("Data Visualization", fluidPage(
                               hr(),
                               tabsetPanel(id="D_visualization",
                                           tabPanel(title="Histogram", value="V_panel1",
                                                    br(),
                                                    uiOutput("markerFilter"),
                                                    hr(),
                                                    actionButton("updateHistogram", "Update Plot", icon = icon("hand-pointer-o")),
                                                    plotOutput("HistogramPlot", width = "100%")
                                                    ),
                                           tabPanel(title="Dot Plot", value="V_panel2",
                                                    br(),
                                                    fluidRow(
                                                        column(4,
                                                               uiOutput("xlabMarker")
                                                        ),
                                                        column(4,
                                                               uiOutput("ylabMarker")
                                                        ),
                                                        column(2,
                                                               numericInput("dotSize", "Dot Size:",
                                                                            value = 0.1,
                                                                            min = 0.1,
                                                                            max = 2,
                                                                            step = 0.1)
                                                        ),
                                                        column(2,
                                                               checkboxInput("freeScales",
                                                                             label = "Free Scales",
                                                                             value = FALSE)
                                                        )),
                                                    hr(),
                                                    actionButton("updateDotPlot", "Update Plot", icon = icon("hand-pointer-o")),
                                                    plotOutput("DotPlot", width = "100%")
                                                    )
                                           )
                           )),
                           tabPanel("Data Normalization", fluidPage(
                               hr(),
                               tabsetPanel(id="D_normalization",
                                           tabPanel(title="GaussNorm", value="N_panel1",
                                                    br(),
                                                    h4("Type in the Number of Landmarks for Each Marker:"),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('MarkerLandmark', i))
                                                    }),
                                                    hr(),
                                                    actionButton("doGaussNorm", "Start Normalization", icon = icon("hand-pointer-o"))
                                                    ),
                                           tabPanel(title="fdaNorm", value="N_panel2")
                                           )
                           ))
                           )
               )
        )
    ))
