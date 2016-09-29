library(flowCore)
library(flowTrans)
library(flowDensity)
library(shiny)
options(shiny.maxRequestSize=100*1024^2) 

plotDens1 <- function(flow.frame, channels, col, trans, pch = ".", ...){
        
        f.exprs <- exprs(flow.frame)
        f.data <- pData(parameters(flow.frame))
        f.col.names <- colnames(flow.frame)
        
        colPalette <- colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
        col <- densCols(f.exprs[,channels], colramp = colPalette)
        
        if(is.numeric(channels))
                channels <- f.col.names[channels]
        xlab <- paste("<", channels[1], ">:", f.data$desc[which(f.col.names==channels[1])], sep = "")
        ylab <- paste("<", channels[2], ">:", f.data$desc[which(f.col.names==channels[2])], sep = "")
        main <- "All Events"
        plot(apply(f.exprs[,channels], 2, trans), col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ...)
}

shinyServer(function(input, output) {

        ## load reactive data from fcs Files
        fcsData <- reactive({
                inFile <- input$fcsFile
                if (is.null(inFile))
                        return(NULL)
                read.FCS(inFile$datapath, transformation = FALSE)   
        })
        
        ## extract markers
        paraData <- reactive({
                if(is.null(fcsData()))
                        return(NULL)
                channels <- colnames(fcsData())
                desname <- pData(fcsData()@parameters)$desc
                markers <- vector()
                for(i in 1:length(channels)){
                        markers[i] <- paste("<", channels[i], ">:", desname[i], sep = "")
                } 
                markers
        })
        
        ## markers selection widget 
        output$marker1_select <- renderUI({
                if(is.null(paraData())){
                        selectInput("marker1", "Choose marker1:", as.list("NULL"))
                }else{
                        selectInput("marker1", "Choose marker1:", as.list(paraData())) 
                }   
        })
                
        output$marker2_select <- renderUI({
                if(is.null(paraData())){
                        selectInput("marker2", "Choose marker2:", as.list("NULL"))
                }else{
                        selectInput("marker2", "Choose marker2:", as.list(paraData())) 
                }   
        })
        
        ## original plot
        output$original_plot <- renderPlot({
                channels <- colnames(fcsData())       
                m1 <- channels[grep(input$marker1, paraData())]
                m2 <- channels[grep(input$marker2, paraData())]
                plotDens(fcsData(), channels=c(m1, m2), asp = 1) 
        })
        
        ## self-defined plot
        output$selDefined_plot <- renderPlot({
                channels <- colnames(fcsData())       
                m1 <- channels[grep(input$marker1, paraData())]
                m2 <- channels[grep(input$marker2, paraData())]
                
                lgcl1 <- logicleTransform( w = input$w, t = input$t, m = input$m, a = input$a)
                plotDens1(fcsData(), channels=c(m1, m2), trans = lgcl1, asp = 1)
        }, height = 600, width = 600)
        
        ## Auto-estimate plot
        output$auto_plot <- renderPlot({
                channels <- colnames(fcsData())       
                m1 <- channels[grep(input$marker1, paraData())]
                m2 <- channels[grep(input$marker2, paraData())]
                
                lgcl2 <- estimateLogicle(fcsData(), channels = c(m1, m2))
                after2 <- transform(fcsData(), lgcl2)
                plotDens(after2, channels=c(m1, m2), asp = 1)         
        }, height = 600, width = 600)
        
        ## MCL optimization plot
        output$mcl_plot <- renderPlot({
                channels <- colnames(fcsData())       
                m1 <- channels[grep(input$marker1, paraData())]
                m2 <- channels[grep(input$marker2, paraData())]
                transformed <- flowTrans(dat=fcsData(),fun="mclMultivBiexp",dims=c(m1, m2),n2f=FALSE,parameters.only=FALSE)
                plotDens(transformed$result, channels=c(m1, m2), asp = 1)    
        }, height = 600, width = 600)
        
        ## MCL optimization plot
        output$density_plot <- renderPlot({
                channels <- colnames(fcsData())       
                m1 <- channels[grep(input$marker1, paraData())]
                m2 <- channels[grep(input$marker2, paraData())]
                par(mfrow = c(1, 2), mar=c(2,2,2,2))
                plot(density(fcsData()@exprs[,m1]))
                plot(density(fcsData()@exprs[,m2]))
        })
        
        ## logicle function curve
        output$lgcl_function_plot <- renderPlot({ 
                lgclTrans <- function(x, w = input$w, t = input$t, m = input$m, a = input$a){
                        x <- .Call("logicle_transform", as.double(x),
                                   as.double(t), as.double(w), 
                                   as.double(m), as.double(a))}
                curve(lgclTrans, from = -10000, to = 500000, 
                      ylim = c(-5, 5),
                      xlab = "original value(x)",
                      ylab = "transformed value(y)",
                      col = "red", lwd = 2)
                abline(v = 0, h = 0, col = "green")             
        })  


})