library(flowCore)

shinyServer(function(input, output) {
        output$funcSelec = renderText({input$funcName}) 
        
        output$funcPlot <- renderPlot({ 
                if (input$funcName == "truncateTransform"){
                        truncateTrans <- function(x, a = input$trun_a){
                                x[x <= a] <- a
                                x}
                        curve(truncateTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth) 
                        abline(v = 0, h = 0, col = "green")
                }else if(input$funcName == "scaleTransform"){
                        scaleTrans <- function(x, a = input$scale_a, b = input$scale_b){
                               x <- (x - a)/(b - a)}
                        curve(scaleTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")
                }else if(input$funcName == "linearTransform"){
                        linearTrans <- function(x, a = input$linear_a, b = input$linear_b) {
                                x <- a + b * x}
                        curve(linearTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")
                }else if(input$funcName == "quadraticTransform"){
                        quadraticTrans <- function(x, a = input$quadratic_a, b = input$quadratic_b, c = input$quadratic_c){
                                x <- a * x^2 + b * x + c}
                        curve(quadraticTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")
                }else if(input$funcName == "lnTransform"){
                        lnTrans <- function(x, r = input$ln_r, d = input$ln_d){
                                x <- log(x) * (r/d)}
                        curve(lnTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")    
                }else if(input$funcName == "logTransform"){
                        logTrans <- function(x, logbase = input$log_base, 
                                             r = input$log_r, d = input$log_d){
                                x <- log(x, logbase) * (r/d)}
                        curve(logTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")      
                }else if(input$funcName == "biexponentialTransform"){
                        biexTrans <- function(x, a = input$bie_a, b = input$bie_b, c = input$bie_c, d = input$bie_d, 
                                              f = input$bie_f, w = input$bie_w, tol = input$bie_tol, 
                                              maxit = input$bie_maxit){
                                x <- .Call("biexponential_transform", x, a, b, c, d, f, w, tol, maxit)
                        }
                        curve(biexTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")      
                }else if(input$funcName == "arcsinhTransform"){
                        arcsTrans <- function(x, a = input$acs_a, b = input$acs_b, c = input$acs_c){
                                x <- asinh(a + b * x) + c
                        } 
                        curve(arcsTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green") 
                }else if(input$funcName == "logicleTransform"){
                        lgclTrans <- function(x, w = input$lgcl_w, t = input$lgcl_t, m = input$lgcl_m, a = input$lgcl_a){
                                x <- .Call("logicle_transform", as.double(x),
                                    as.double(t), as.double(w), 
                                    as.double(m), as.double(a))}
                        curve(lgclTrans, from = input$xfrom, to = input$xto, 
                              ylim = c(input$yfrom, input$yto),
                              xlab = "original value(x)",
                              ylab = "transformed value(y)",
                              col = input$lineColor, lwd = input$lineWidth)
                        abline(v = 0, h = 0, col = "green")       
                }
                })  
})