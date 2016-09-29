## max data size
options(shiny.maxRequestSize=1024^3)

shinyServer(function(input, output, session) {
    v <- reactiveValues(exprsData = NULL, normExprsData = NULL)

    ##-------------------Side Panel-------------------

    observeEvent(input$loadButton, {
        fcsFiles <- input$fcsFiles
        compFile <- input$compensationFile
        if (is.null(fcsFiles)){
            v$exprsData <- NULL
        }else{
            withProgress(message="Loading and Processing Data...", value=0, {
                #print(fcsFiles$datapath)
                print(fcsFiles$name)

                if(is.null(compFile)){
                    compMatrix <- FALSE
                }else{
                    compMatrix <- read.csv(compFile$datapath,
                                           check.names = FALSE,
                                           header = TRUE, row.names = 1)
                    tidynames <- sub("\\s\\:\\:\\s.+$", "", colnames(compMatrix), perl = TRUE)
                    row.names(compMatrix) <- tidynames
                    colnames(compMatrix) <- tidynames
                    compMatrix <- as.matrix(compMatrix)
                }
                data <- cytof_exprsMerge(fcsFiles = fcsFiles$datapath,
                                         markers = NULL,
                                         comp = compMatrix,
                                         transformMethod = input$transMethod,
                                         mergeMethod = "all")

                ## update rownames with real original file names
                suffix <- do.call(c, regmatches(row.names(data), gregexpr("_[0-9]*$", row.names(data), perl=T)))
                prefix <- sub(".fcs$", "", fcsFiles$name)[as.numeric(sub("_[0-9]*$", "", row.names(data)))+1]
                row.names(data) <- paste0(prefix, suffix)
                v$exprsData <- data
            })
        }
    })

    output$normDataCheck <- renderUI({
        checkboxInput("useNormData",
                      label = "Normalized Data",
                      value = !is.null(v$normExprsData))
    })

    output$sampleFilter <- renderUI({
        if(is.null(v$exprsData)){
            return(NULL)
        }else{
            sampleNames <- unique(as.character(sub("_[0-9]*$", "", row.names(v$exprsData))))
            checkboxGroupInput('samples', NULL,
                               sampleNames,
                               selected = sampleNames)
        }
    })

    observeEvent(input$saveButton, {
        if(!is.null(input$fcsFiles) && !is.null(v$normExprsData)){
            withProgress(message="Saving Normalized FCS files...", value=0, {
                print(getwd())
                dir.create("Normalized_FCS_files")
                resultDir <- paste0(getwd(), .Platform$file.sep, "Normalized_FCS_files")
                fcsPaths <- input$fcsFiles$datapath
                cellFromFileName <- sub("_[0-9]*$", "", row.names(v$normExprsData))

                for(i in seq_along(fcsPaths)){
                    fi <- read.FCS(fcsPaths[i], transformation = FALSE)
                    fi_name <- sub(".fcs$", "", input$fcsFiles$name[i])
                    fi_exprs <- v$normExprsData[cellFromFileName == fi_name, ]
                    markerMatchID <- match(sub("\\<.+\\>$", "",colnames(fi_exprs), perl = TRUE),
                                           colnames(fi@exprs))
                    print(markerMatchID)
                    fi@exprs[ ,markerMatchID] <- fi_exprs

                    filename <- paste0(resultDir, .Platform$file.sep,
                                       fi_name, "_normalized.fcs")
                    write.FCS(fi, filename = filename)
                }
            })

            ## open the results directory
            opendir(resultDir)
        }
    })

    ##---------------Data Visualization tabset-------------------

    ## Histogram

    output$markerFilter <- renderUI({
        if(is.null(v$exprsData)){
            return(NULL)
        }else{
            markers <- colnames(v$exprsData)
            selectizeInput('markers', 'Choose Markers:',
                           choices = markers, selected = markers,
                           multiple = TRUE, width = "100%")
        }
    })

    observeEvent(input$updateHistogram, {
        markers <- isolate(input$markers)
        samples <- isolate(input$samples)
        output$HistogramPlot <- renderPlot({
            if(is.null(v$exprsData) || is.null(markers) || is.null(samples)){
                return(NULL)
            }else{
                if(input$useNormData && !is.null(v$normExprsData)){
                    data <- v$normExprsData
                }else{
                    data <- v$exprsData
                }
                withProgress(message="Generating Histogram Plot...", value=0, {
                    cellSample <- sub("_[0-9]*$", "", row.names(data))
                    sampleCheck <- cellSample %in% samples
                    gp <- stackDenistyPlot(data = data[sampleCheck,],
                                           densityCols = markers,
                                           stackFactor = cellSample[sampleCheck],
                                           x_text_size = 2,
                                           legend_text_size = 0.5,
                                           legendRow = ceiling(length(samples)/4),
                                           legend_title = "sample")
                    plot(gp)
                })
            }
        }, height = 800, width = 850)
    })

    ## Dot plot

    output$xlabMarker <- renderUI({
        if(is.null(v$exprsData)){
            return(NULL)
        }else{
            markers <- colnames(v$exprsData)
            selectInput('marker1', 'Choose X Marker:',
                        choices = markers,
                        selected = markers[1],
                        width = "100%")
        }
    })

    output$ylabMarker <- renderUI({
        if(is.null(v$exprsData)){
            return(NULL)
        }else{
            markers <- colnames(v$exprsData)
            selectInput('marker2', 'Choose Y Marker:',
                        choices = markers,
                        selected = markers[2],
                        width = "100%")
        }
    })

    observeEvent(input$updateDotPlot, {
        marker1 <- isolate(input$marker1)
        marker2 <- isolate(input$marker2)
        samples <- isolate(input$samples)
        output$DotPlot <- renderPlot({
            if(is.null(v$exprsData) || is.null(marker1) ||
               is.null(marker2) || is.null(samples)){
                return(NULL)
            }else{
                if(input$useNormData && !is.null(v$normExprsData)){
                    data <- v$normExprsData
                }else{
                    data <- v$exprsData
                }
                withProgress(message="Generating Dot Plot...", value=0, {
                    cellSample <- sub("_[0-9]*$", "", row.names(data))
                    sampleCheck <- cellSample %in% samples
                    gp <- dotPlot(data = data[sampleCheck, ],
                                  xlab = marker1,
                                  ylab = marker2,
                                  pointSize = input$dotSize,
                                  freeScales = input$freeScales)
                    plot(gp)
                })
            }
        }, height = 800, width = 850)
    })


    ##---------------Data Normalization tabset-------------------

    lapply(1:100, function(i) {
        output[[paste0('MarkerLandmark', i)]] <- renderUI({
            if(is.null(v$exprsData)){
                return(NULL)
            }

            markers <- colnames(v$exprsData)
            if (i <= length(markers)){
                x <- markers[i]
                textInput(paste0('marker_lms_', i), x,
                          value = "", width = "30%",
                          placeholder = "Type in the number of landmarks")
            }
        })
    })

    observeEvent(input$doGaussNorm,{
        if(!is.null(v$exprsData)){
            withProgress(message="Performing Gaussian Normalization...", value=0, {
                landmarks <- NULL
                markers <- colnames(v$exprsData)
                for (i in 1:length(markers)){
                    ilms <- input[[paste0('marker_lms_', i)]]
                    ilms <- ifelse(ilms=="", NA, as.numeric(ilms))
                    landmarks <- c(landmarks, ilms)}
                ifNormMarker <- !is.na(landmarks)
                cellSample <- sub("_[0-9]*$", "", row.names(v$exprsData))
                sampleCheck <- cellSample %in% input$samples

                data_norm <- gaussNorm(exprs=v$exprsData[sampleCheck, ],
                                       samplesInfo=cellSample[sampleCheck],
                                       channel.names=markers[ifNormMarker],
                                       max.lms=landmarks[ifNormMarker])
                v$normExprsData <- data_norm

                updateTabsetPanel(session, "D_visualization", selected = "V_panel1")
            })
        }
    })
})
