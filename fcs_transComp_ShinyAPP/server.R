## max data size
options(shiny.maxRequestSize=1024^3)

shinyServer(function(input, output, session) {

    ##-------------------------------------------------------------------------

    ## load flowframe data
    fcs <- reactive({
        if (input$goButton == 0)
            return()
        isolate({fcsFile <- input$fcsFile
                 if (is.null(fcsFile))
                     return(NULL)
                 fcs <- read.FCS(fcsFile$datapath, transformation = FALSE)
                 fcs@description$FILENAME <- fcsFile$name
                 })
        return(fcs)
    })

    ## channel names and marker names
    markers <- reactive({
        if(is.null(fcs()))
            return(NULL)
        fcs <- fcs()
        iChannel <- grep("Time|Event|Length", pData(fcs@parameters)$name, ignore.case = TRUE)
        markers <- paste0(pData(fcs@parameters)$name,
                         "<",
                         pData(fcs@parameters)$desc,
                         ">")[setdiff(1:ncol(fcs), iChannel)]
        return(markers)
    })


    channels <- reactive({
        if(is.null(fcs()))
            return(NULL)
        fcs <- fcs()
        iChannel <- grep("Time|Event|Length", pData(fcs@parameters)$name, ignore.case = TRUE)
        channels <- pData(fcs@parameters)$name[setdiff(1:ncol(fcs), iChannel)]
        names(channels) <- NULL
        return(channels)
    })

    output$markerChoose <- renderUI({
        if(is.null(fcs()) || is.null(markers())){
            return(NULL)
        }else{
            selectizeInput('transMarkers', 'Choose Markers:',
                           choices = markers(), selected = markers(),
                           multiple = TRUE, width = "100%")
        }
    })


    gp <- eventReactive(input$goButton2, {
        if(is.null(fcs()) || is.null(markers) || is.null(channels()))
            return(NULL)

        transChannels <- channels()[match(input$transMarkers, markers())]

        withProgress(message="Creating Plot in Progress...", value=0, {
            ## Transformation 1
            switch(input$transformation1,
                   autoLgcl={
                       trans1 <- autoLgcl(fcs(),
                                          channels = transChannels,
                                          m = input$t1_autoLgcl_m,
                                          q = input$t1_autoLgcl_q)
                       fcs1 <- transform(fcs(), trans1)
                   },
                   logicle={
                       lgcl <- logicleTransform(w = input$t1_logicle_w,
                                                t = input$t1_logicle_t,
                                                m = input$t1_logicle_m,
                                                a = input$t1_logicle_a)
                       trans1 <- transformList(transChannels, lgcl)
                       fcs1 <- transform(fcs(), trans1)
                   },
                   cytofAsinh={
                       fcs1 <- fcs()
                       data <- fcs1@exprs
                       data[ ,transChannels] <- apply(data[ ,transChannels, drop=FALSE], 2,
                                                      function(x,
                                                               a=input$t1_cytofAsinh_a,
                                                               b=input$t1_cytofAsinh_b,
                                                               c=input$t1_cytofAsinh_c) {
                                                          x <- x-1
                                                          loID <- which(x < 0)
                                                          if(length(loID) > 0)
                                                              x[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)

                                                          x <- asinh(a + b * x) + c  # x <- log(x + sqrt(x^2 + 1))
                                                          return(x)
                                                          }
                                                      )
                       fcs1@exprs <- data
                   },
                   arcsinh={
                       asinhTrans <- arcsinhTransform(a = input$t1_arcsinh_a,
                                                      b = input$t1_arcsinh_b,
                                                      c = input$t1_arcsinh_c)
                       trans1 <- transformList(transChannels, asinhTrans)
                       fcs1 <- transform(fcs(), trans1)
                   })

            ## Transformation 2
            switch(input$transformation2,
                   autoLgcl={
                       trans2 <- autoLgcl(fcs(),
                                          channels = transChannels,
                                          m = input$t2_autoLgcl_m,
                                          q = input$t2_autoLgcl_q)
                       fcs2 <- transform(fcs(), trans2)
                   },
                   logicle={
                       lgcl <- logicleTransform(w = input$t2_logicle_w,
                                                t = input$t2_logicle_t,
                                                m = input$t2_logicle_m,
                                                a = input$t2_logicle_a)
                       trans2 <- transformList(transChannels, lgcl)
                       fcs2 <- transform(fcs(), trans2)
                   },
                   cytofAsinh={
                       fcs2 <- fcs()
                       data <- fcs2@exprs
                       data[ ,transChannels] <- apply(data[ ,transChannels, drop=FALSE], 2,
                                                      function(x,
                                                               a=input$t2_cytofAsinh_a,
                                                               b=input$t2_cytofAsinh_b,
                                                               c=input$t2_cytofAsinh_c) {
                                                          x <- x-1
                                                          loID <- which(x < 0)
                                                          if(length(loID) > 0)
                                                              x[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)

                                                          x <- asinh(a + b * x) + c  # x <- log(x + sqrt(x^2 + 1))
                                                          return(x)
                                                      }
                       )
                       fcs2@exprs <- data
                   },
                   arcsinh={
                       asinhTrans <- arcsinhTransform(a = input$t2_arcsinh_a,
                                                      b = input$t2_arcsinh_b,
                                                      c = input$t2_arcsinh_c)
                       trans2 <- transformList(transChannels, asinhTrans)
                       fcs2 <- transform(fcs(), trans2)
                   })

            flowSet <- list()
            transName1 <- input$transformation1
            transName2 <- input$transformation2

            if(transName1 == transName2){
                if(!(identical(fcs1@exprs, fcs2@exprs))){
                    transName1 <- paste0(transName1, "_1")
                    transName2 <- paste0(transName1, "_2")
                }
            }
            flowSet[[transName1]] <- fcs1
            flowSet[[transName2]] <- fcs2
            flowSet <- as(flowSet, "flowSet")

            ## creat plot
            switch(input$plotType,
                   DensityPlot = {
                       gp <- flowSet_stackDenistyPlot(flowSet,
                                                      channels = transChannels,
                                                      kernel = "gaussian",
                                                      bw = "nrd0",
                                                      adjust = 1,
                                                      stackRotation = 0,
                                                      stackSeperation = 0.8,
                                                      scaleHeight = input$scaleHeight,
                                                      scaleWidth = input$scaleWidth,
                                                      x_text_size = 1,
                                                      strip_text_size = 8,
                                                      legend_text_size = 1.5,
                                                      legend_title = "")
                   },
                   boxplot = {
                       gp <- flowSet_wrapBWplot(flowSet,
                                                channels = transChannels,
                                                box_ymin = 0.05,
                                                box_lower = 0.25,
                                                box_middle = 0.5,
                                                box_upper = 0.75,
                                                box_ymax = 0.95,
                                                x_text_size = 2,
                                                strip_text_size = 8,
                                                legend_text_size = 1.5,
                                                legend_title = "")
                   }
            )
        })

        return(gp)
    })


    output$comparationPlot <- renderPlot({
        print(gp())
    })


    output$summaryTable <- renderDataTable({
        if(is.null(fcs()))
            return(NULL)
        para <- pData(fcs()@parameters)
        para$minExprs <- apply(fcs()@exprs, 2, min)
        para$maxExprs <- apply(fcs()@exprs, 2, max)

        para
    })

})

