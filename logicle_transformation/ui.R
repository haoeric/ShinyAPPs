
library(shiny)

shinyUI(pageWithSidebar(
        
        # Application title
        headerPanel("Optimizing parameters for logicle transformation"),
        
        sidebarPanel(
                fileInput('fcsFile', 'Choose FCS file:', accept=c('.fcs')),
                
                uiOutput("marker1_select"),
                
                uiOutput("marker2_select"),
                        
#                 radioButtons("plotMethos", "Plot Method select:",
#                         c("Contour" = "con",
#                         "densityPlot" = "dens",
#                         "Exponential" = "exp")),
#                 checkboxGroupInput('plot_option', 'Plot Options:',
#                                    c("Original Data","sel-defined","auto estimation","MCL estimation"), 
#                                    selected = c("Original Data","sel-defined","auto estimation","MCL estimation")),
                
                helpText("Parameters selection for sel-defined logicle transform"),
                
                numericInput("t", 
                             label = "t, top of the scale data value",
                             value = 262144),
                
                sliderInput("m", 
                            "m, the full width(decades):", 
                            value = 4.5, min = 1, max = 10, step = 0.1),
                
                sliderInput("w", 
                            "w, the linearization width(decades):", 
                            value = 0.5, min = 0, max = 5, step = 0.1),
                
                sliderInput("a", 
                            "a, additional negative range(decades):", 
                            value = 0, min = 0, max = 5, step = 0.1),
                
                submitButton("Update View")
                        
                        
                ),
                
                mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel("Original Plot", plotOutput("original_plot")),
                                    tabPanel("SelDefined Plot", plotOutput("selDefined_plot", height="auto", width = "100%")),
                                    tabPanel("Auto Plot", plotOutput("auto_plot", height="auto", width = "100%")),
                                    tabPanel("MCL Plot", plotOutput("mcl_plot", height="auto", width = "100%")),
                                    tabPanel("Density plot", plotOutput("density_plot")),
                                    tabPanel("Logicle function", plotOutput("lgcl_function_plot"))
                        )
                )
        )
)