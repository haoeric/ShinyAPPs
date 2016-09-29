
shinyUI(fluidPage(
        titlePanel("Data Transformation Functions"),
        
        sidebarLayout(
                sidebarPanel(
                        helpText("Plot the function of a transformation method."),
                        
                        selectInput("funcName", 
                                    label = "Choose a transformation method to display",
                                    choices = c("truncateTransform", 
                                                "scaleTransform",
                                                "linearTransform", 
                                                "quadraticTransform",
                                                "lnTransform",
                                                "logTransform",
                                                "arcsinhTransform",
                                                "biexponentialTransform",
                                                "logicleTransform"),
                                    selected = "truncateTransform"
                        ),
                        
                        selectInput("lineColor", 
                                    label = "The color of the function line",
                                    choices = c("red", 
                                                "black",
                                                "blue"),
                                    selected = "red"
                        ),
                        
                        numericInput("lineWidth", 
                                     label = "The width of the function line", 
                                     value = 2),
                
                        numericInput("xfrom", 
                                     label = "The x range starts from:", 
                                     value = -100),
                        
                        numericInput("xto", 
                                     label = "The x range ends to:", 
                                     value = 100),
                        
                        numericInput("yfrom", 
                                     label = "The y range starts from:", 
                                     value = -100),
                        
                        numericInput("yto", 
                                     label = "The y range ends to:", 
                                     value = 100)
                ),
                
                mainPanel(
                        textOutput("funcSelec"),
                        plotOutput("funcPlot"),
                        conditionalPanel("input.funcName == 'truncateTransform'", 
                                         numericInput("trun_a", 
                                                      label = "a, the value at which to truncate", 
                                                      value = 1)),
                        
                        conditionalPanel("input.funcName == 'scaleTransform'", 
                                         numericInput("scale_a", 
                                                      label = "a, value that will be transformed to 0", 
                                                      value = 0),
                                         
                                         numericInput("scale_b", 
                                                      label = "b, value that will be transformed to 1",
                                                      value = 100)),
                        
                        conditionalPanel("input.funcName == 'linearTransform'", 
                                         numericInput("linear_a", 
                                                      label = "a, additive factor in the equation", 
                                                      value = 0),
                                         
                                         numericInput("linear_b", 
                                                      label = "b, the multiplicative factor in the equation",
                                                      value = 5)),
                        
                        conditionalPanel("input.funcName == 'quadraticTransform'", 
                                         numericInput("quadratic_a", 
                                                      label = "a, the quadratic coefficient in the equation", 
                                                      value = 1),
                                         
                                         numericInput("quadratic_b", 
                                                      label = "b, the linear coefficient in the equation",
                                                      value = 2),
                                         
                                         numericInput("quadratic_c", 
                                                      label = "c, the intercept in the equation",
                                                      value = 5)),
                        
                        conditionalPanel("input.funcName == 'lnTransform'", 
                                         numericInput("ln_r", 
                                                      label = "r, corresponds to a scale factor",
                                                      value = 2),
                                         numericInput("ln_d", 
                                                      label = "d, corresponds to a scale factor",
                                                      value = 5)),
                        
                        conditionalPanel("input.funcName == 'logTransform'", 
                                         numericInput("log_base", 
                                                      label = "base, corresponds to a scale factor",
                                                      value = 10),
                                         numericInput("log_r", 
                                                      label = "r, corresponds to a scale factor",
                                                      value = 2),
                                         numericInput("log_d", 
                                                      label = "d, corresponds to a scale factor",
                                                      value = 5)),
                        
                        conditionalPanel("input.funcName == 'arcsinhTransform'", 
                                         numericInput("acs_a", 
                                                      label = "a, corresponds to the base of the logarithm", 
                                                      value = 0),
                                         
                                         numericInput("acs_b", 
                                                      label = "b, corresponds to a scale factor",
                                                      value = 1),
                                         
                                         numericInput("acs_c", 
                                                      label = "c, corresponds to a scale factor",
                                                      value = 0)),
                        
                        conditionalPanel("input.funcName == 'biexponentialTransform'",    
                                         fluidRow(
                                                 column(3,
                                                        numericInput("bie_a", 
                                                                     label = "a", 
                                                                     value = 0.5),
                                                        numericInput("bie_d", 
                                                                     label = "d",
                                                                     value = 1),
                                                        numericInput("bie_tol", 
                                                                     label = "tol, a tolerance to pass to the inversion routine",
                                                                     value = .Machine$double.eps^0.25)),
                                                 column(4, offset = 1,
                                                        numericInput("bie_b", 
                                                                     label = "b",
                                                                     value = 1),
                                                        numericInput("bie_f", 
                                                                     label = "f, a constant bias for the intercept",
                                                                     value = 0),
                                                        numericInput("bie_maxit", 
                                                                     label = "maxit, a maximum number of iterations to use",
                                                                     value = as.integer(5000))),
                                                 column(4,
                                                        numericInput("bie_c", 
                                                                     label = "c",
                                                                     value = 0.5),
                                                        numericInput("bie_w", 
                                                                     label = "w, a constant bias for the 0 point of the data",
                                                                     value = 0))
                                                 )
                                         ),
                         
                        conditionalPanel("input.funcName == 'logicleTransform'", 
                                         fluidRow(
                                                 column(3,
                                                        numericInput("lgcl_w", 
                                                                     label = "w, the linearization width in asymptotic decades", 
                                                                     value = 0.5),
                                                        numericInput("lgcl_m", 
                                                                     label = "m, the full width of the transformed display in asymptotic decades",
                                                                     value = 4.5)),
                                                 column(4,offset = 2,
                                                        numericInput("lgcl_t", 
                                                                     label = "t, top of the scale data value",
                                                                     value = 262144),
                                                        numericInput("lgcl_a", 
                                                                     label = "a, additional negative range to be included in the display in asymptotic decades",
                                                                     value = 0))
                                         )
                       
                ))
        )
))