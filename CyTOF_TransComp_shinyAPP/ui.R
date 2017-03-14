library(shiny)

shinyUI(fluidPage(
    titlePanel("Comparation of Transformation Methods for Flow Cytometry Data"),
    hr(),

    fluidRow(
        column(3,
               fileInput('fcsFile', strong('Choose FCS files:'), multiple = FALSE,
                         accept = c('text/fcs', '.fcs')),
               actionButton("goButton", "Load Data", icon = icon("hand-o-right")),
               hr(),

               selectInput('plotType', 'Visualization Method:',
                           choices = c("DensityPlot", "boxplot", "DataSummary"),
                           selected = "Density",
                           width = "100%"),
               uiOutput("markerChoose"),

               checkboxInput("scaleHeight", label = "Scale Density Height", value = TRUE),
               checkboxInput("scaleWidth", label = "Scale Density Width", value = TRUE),
               actionButton("goButton2", "Creat Plot", icon = icon("hand-o-right")),

               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ",
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
               fluidRow(
                   column(6,
                          selectInput('transformation1', 'Transformation Method 1:',
                                      choices = c("autoLgcl", "logicle", "cytofAsinh", "arcsinh", "autoAsinh"),
                                      selected = "autoLgcl",
                                      width = "100%"),

                          conditionalPanel("input.transformation1 == 'autoLgcl'",
                                           fluidRow(
                                               column(6,
                                                      numericInput("t1_autoLgcl_m", "m:", value = 4.5,
                                                                   min = 4, max = 10, step = 0.1)
                                               ),
                                               column(6,
                                                      numericInput("t1_autoLgcl_q", "q:", value = 0.05,
                                                                   min = 0, max = 0.5, step = 0.01)
                                               )
                                           )
                                           ),
                          conditionalPanel("input.transformation1 == 'logicle'",
                                           fluidRow(
                                               column(3,
                                                      numericInput("t1_logicle_w", "w:", value = 0.5)
                                               ),
                                               column(3,
                                                      numericInput("t1_logicle_t", "t:", value = 10000)
                                               ),
                                               column(3,
                                                      numericInput("t1_logicle_m", "m:", value = 4.5)
                                               ),
                                               column(3,
                                                      numericInput("t1_logicle_a", "a:", value = 0)
                                               )
                                           )
                                           ),
                          conditionalPanel("input.transformation1 == 'cytofAsinh'",
                                           fluidRow(
                                               column(4,
                                                      numericInput("t1_cytofAsinh_a", "a:", value = 0)
                                               ),
                                               column(4,
                                                      numericInput("t1_cytofAsinh_b", "b:", value = 0.2)
                                               ),
                                               column(4,
                                                      numericInput("t1_cytofAsinh_c", "c:", value = 0)
                                               )
                                           )
                                           ),
                          conditionalPanel("input.transformation1 == 'arcsinh'",
                                           fluidRow(
                                               column(4,
                                                      numericInput("t1_arcsinh_a", "a:", value = 1)
                                               ),
                                               column(4,
                                                      numericInput("t1_arcsinh_b", "b:", value = 1)
                                               ),
                                               column(4,
                                                      numericInput("t1_arcsinh_c", "c:", value = 1)
                                               )
                                           )
                                           )
                          ),
                   column(6,
                          selectInput('transformation2', 'Transformation Method 2:',
                                      choices = c("cytofAsinh", "arcsinh", "autoLgcl", "logicle", "autoAsinh"),
                                      selected = "cytofAsinh",
                                      width = "100%"),

                          conditionalPanel("input.transformation2 == 'autoLgcl'",
                                           fluidRow(
                                               column(6,
                                                      numericInput("t2_autoLgcl_m", "m:", value = 4.5,
                                                                   min = 4, max = 10, step = 0.1)
                                               ),
                                               column(6,
                                                      numericInput("t2_autoLgcl_q", "q:", value = 0.05,
                                                                   min = 0, max = 0.5, step = 0.01)
                                               )
                                           )
                          ),
                          conditionalPanel("input.transformation2 == 'logicle'",
                                           fluidRow(
                                               column(3,
                                                      numericInput("t2_logicle_w", "w:", value = 0.5)
                                               ),
                                               column(3,
                                                      numericInput("t2_logicle_t", "t:", value = 10000)
                                               ),
                                               column(3,
                                                      numericInput("t2_logicle_m", "m:", value = 4.5)
                                               ),
                                               column(3,
                                                      numericInput("t2_logicle_a", "a:", value = 0)
                                               )
                                           )
                          ),
                          conditionalPanel("input.transformation2 == 'cytofAsinh'",
                                           fluidRow(
                                               column(4,
                                                      numericInput("t2_cytofAsinh_a", "a:", value = 0)
                                               ),
                                               column(4,
                                                      numericInput("t2_cytofAsinh_b", "b:", value = 0.2)
                                               ),
                                               column(4,
                                                      numericInput("t2_cytofAsinh_c", "c:", value = 0)
                                               )
                                           )
                          ),
                          conditionalPanel("input.transformation2 == 'arcsinh'",
                                           fluidRow(
                                               column(4,
                                                      numericInput("t2_arcsinh_a", "a:", value = 1)
                                               ),
                                               column(4,
                                                      numericInput("t2_arcsinh_b", "b:", value = 1)
                                               ),
                                               column(4,
                                                      numericInput("t2_arcsinh_c", "c:", value = 1)
                                               )
                                           )
                          )
                          )

        ),
        hr(),
        conditionalPanel("input.plotType == 'DensityPlot' || input.plotType == 'boxplot'",
                         plotOutput("comparationPlot", width = "90%", height = "700px")

        ),
        conditionalPanel("input.plotType == 'DataSummary'",
                         dataTableOutput("summaryTable")
        )
    ))
))
