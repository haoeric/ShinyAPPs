
library(flowCore)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)


cytofAsinh <- function(x, a=0, b=0.2, c=0) {
    x <- x-1
    loID <- which(x < 0)
    if(length(loID) > 0)
        x[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)

    x <- asinh(a + b * x) + c  # x <- log(x + sqrt(x^2 + 1))
    return(x)
}


autoLgcl <- function(x, channels, m = 4.5, q = 0.05) {
    if (!is(x, "flowFrame"))
        stop("x has to be an object of class \"flowFrame\"")
    if (missing(channels))
        stop("Please specify the channels to be logicle transformed")
    indx <- channels %in% colnames(x@exprs)
    if (!all(indx))
        stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ",
                   sep = " "))

    trans <- lapply(channels, function(p) {
        data <- x@exprs[, p]
        w <- 0
        t <- max(data)
        ndata <- data[data < 0]
        ## use 1.5 * IQR to filter outliers in negative values
        nThres <- quantile(ndata, 0.25) - 1.5 * IQR(ndata)
        ndata <- ndata[ndata >= nThres]
        transId <- paste(p, "autolgclTransform", sep = "_")

        if (length(ndata)) {
            r <- .Machine$double.eps + quantile(ndata, q)
            ## Check to avoid failure of negative w
            if (10^m * abs(r) <= t) {
                w <- 0
            } else {
                w <- (m - log10(t/abs(r)))/2
                if(w>2) {
                    w <- 0.5
                }
            }
        }
        logicleTransform(transformationId = transId,
                         w = w, t = t, m = m, a = 0)
    })
    transformList(channels, trans)
}


flowSet_stackDenistyPlot <- function(flowSet,
                                     channels,
                                     kernel = c("gaussian", "epanechnikov", "rectangular",
                                                "triangular", "biweight",
                                                "cosine", "optcosine"),
                                     bw = "nrd0", adjust = 1,
                                     stackRotation = 0, stackSeperation = "auto",
                                     scaleHeight = FALSE, scaleWidth = FALSE,
                                     x_text_size = 2, strip_text_size = 7,
                                     legend_text_size = 0.5, legend_title = "stackName"){

    marker_id <- match(channels, colnames(flowSet))

    if(any(is.na(marker_id))){
        stop(paste0("Can not find the specified channels:",
                    paste(channels[is.na(marker_id)], collapse = " "),
                    " in the flowSet data!"))
    }

    data <- NULL
    stackFactor <- NULL
    for(sample in names(flowSet@frames)){
        fcs <- flowSet@frames[[sample]]
        pd <- fcs@parameters@data
        exprs <- fcs@exprs[ ,marker_id, drop=FALSE]
        colnames(exprs) <- paste0(pd$name, "<", pd$desc,">")[marker_id]
        data <- rbind(data, exprs)
        stackFactor <- c(stackFactor, rep(sample, nrow(exprs)))
    }

    stackDenistyPlot(data = data,
                     stackFactor = stackFactor,
                     kernel = kernel,
                     bw = bw,
                     adjust = adjust,
                     stackRotation = stackRotation,
                     stackSeperation = stackSeperation,
                     scaleHeight = scaleHeight,
                     scaleWidth = scaleWidth,
                     x_text_size = x_text_size,
                     strip_text_size = strip_text_size,
                     legend_text_size = legend_text_size,
                     legend_title = legend_title)
}


#' @param densData Data frame.
#' @param stackRotation Rotation degree of density plot to the right side, range (0-90).
#' @param stackSeperation Control factor for stack seperation interval, numeric value from 0-1, or auto.
#'
#' @import plyr ldply
#' @import reshape2 melt
stackDenistyPlot <- function(data, densityCols, stackFactor,
                             kernel = c("gaussian", "epanechnikov", "rectangular",
                                        "triangular", "biweight",
                                        "cosine", "optcosine"),
                             bw = "nrd0", adjust = 1,
                             stackRotation = 0, stackSeperation = "auto",
                             scaleHeight = FALSE, scaleWidth = FALSE,
                             x_text_size = 2, strip_text_size = 7,
                             legend_text_size = 0.5, legend_title = "stackName"){

    if(!is.numeric(stackRotation)){
        stop("stackRotation must be a numeric number")
    }else if(stackRotation < 0 || stackRotation > 90){
        stop("stackRotation must be a numeric number in range 0-90")
    }

    if(missing(densityCols)){
        densityCols <- colnames(data)
    }else if(any(!(densityCols %in% colnames(data)))){
        stop("Unmatch densityCols found:", paste(densityCols[!(densityCols %in% colnames(data))], collapse = " "))
    }

    if(missing(stackFactor)){
        warning("no stackFactor was provided!")
        stackFactor <- rep("stack", length = nrow(data))
    }else if(length(stackFactor) != nrow(data)){
        stop("Length of stackFactor unequal row number of input data")
    }
    kernel <- match.arg(kernel)

    stackCount <- length(unique(stackFactor))
    densityCount <- length(densityCols)

    data <- data.frame(data[ ,densityCols, drop=FALSE], stackFactor = stackFactor, check.names = FALSE)
    densData <- .densityCal(data, kernel = kernel, bw = bw, adjust = adjust, scaleHeight = scaleHeight, scaleWidth = scaleWidth)
    ## dataframe densData contains {stackName, x , y , densityName}
    xStat <- aggregate(x ~ stackName + densityName, densData, max)
    yStat <- aggregate(y ~ stackName + densityName, densData, max)

    if(stackSeperation == "auto"){
        stackIntervals <- aggregate(y ~ densityName, yStat, function(x){0.8*median(x) * (1-(stackRotation/90)^0.2)^2})
    }else if(stackSeperation < 0 || stackSeperation > 1){
        stop("stackSeperation must be value in range 0-1")
    }else{
        stackIntervals <- aggregate(y ~ densityName, yStat, function(x){median(x)*stackSeperation})
    }

    stackShifts <- aggregate(x ~ densityName, xStat, function(x){max(x) * (stackRotation/90)})

    densData$stack_x <- densData$x + (as.numeric(densData$stackName)-1) * stackShifts$x[match(densData$densityName, stackShifts$densityName)]
    densData$stack_y <- densData$y + (as.numeric(densData$stackName)-1) * stackIntervals$y[match(densData$densityName, stackIntervals$densityName)]

    ## segment lines, x tick, x label
    alignSegments <- ldply(split(densData$x, densData$densityName),
                           function(x){seq(min(x), max(x), length.out=5)},
                           .id = "densityName")
    alignSegments <- melt(alignSegments, id.vars="densityName", variable.name="x_tick", value.name = "x")
    alignSegments$y <- min(densData$y)
    alignSegments$xend <- alignSegments$x + (length(unique(densData$stackName))-1) * stackShifts$x[match(alignSegments$densityName, stackShifts$densityName)]
    alignSegments$yend <- min(densData$y) + (length(unique(densData$stackName))-1) * stackIntervals$y[match(alignSegments$densityName, stackIntervals$densityName)]

    densityHeights <- aggregate(y ~ densityName, yStat, max)
    alignSegments$tickXend <- alignSegments$x
    alignSegments$tickYend <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.01
    alignSegments$tickText <- format(alignSegments$x,scientific=TRUE, digits=3)
    alignSegments$textY <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.03

    cat(" Plotting ...\n")
    stackDensityPlot_theme <- theme(legend.position = "top",
                                    legend.title = element_text(size = rel(0.6)),
                                    legend.text = element_text(size = rel(legend_text_size)),
                                    legend.background = element_rect(fill=NULL, size=0.8, linetype="solid"),
                                    strip.text = element_text(size=strip_text_size, lineheight=0.1, hjust = 0.5, vjust = 0.5),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank(),
                                    strip.background=element_rect(fill = "grey90", colour = NA))

    gp <- ggplot(densData, aes(x=stack_x, y=stack_y)) +
        geom_segment(data = alignSegments,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     color = "grey80", size=0.3) +
        geom_segment(data = alignSegments,
                     aes(x = x, y = y, xend = tickXend, yend = tickYend),
                     color = "grey20", size=0.3) +
        geom_text(data = alignSegments, aes(x = x, y = textY, label = tickText),
                  hjust = 0.3, vjust = 1.1, size = x_text_size) +
        geom_polygon(aes(fill=stackName, color=stackName), alpha = 0.15) +
        facet_wrap(~densityName, scale = "free") +
        xlab("") + ylab("") +
        guides(col = guide_legend(title.hjust = 0.5, title = legend_title),
               fill = guide_legend(keywidth = 0.5, keyheight = 0.5,
                                   title.position = "top", title = legend_title)) +
        theme_bw() + stackDensityPlot_theme

    gp
}


#' Internal density calculation function serves for \code{stackDenistyPlot}
#'
#' Output data frame with columns: stackName, x , y , densityName
.densityCal <- function(data, kernel, bw, adjust, scaleHeight = FALSE, scaleWidth = FALSE){
    cat("  Calculating Density for each stack column...\n")
    dataBystackFactor <- split(subset(data, select = -stackFactor), data$stackFactor)
    densityWrap <- function(d, ...){
        resOut <- NULL
        for(i in colnames(d)){
            x <- d[,i]
            dens <- density(x, ...)
            dx <- dens$x
            dy <- dens$y
            if(scaleWidth){
                dx <- rescale(dx, c(0,5))
            }
            if(scaleHeight){
                dy <- rescale(dy, c(0,1))
            }
            densOut <- data.frame(x=dx, y=dy, densityName = i)
            resOut <- rbind(resOut, densOut)
        }
        return(resOut)
    }

    r <- ldply(dataBystackFactor, densityWrap,
               kernel = kernel, bw = bw, adjust = adjust,
               .progress = "text",
               .id = "stackName")
    return(r)
}



flowSet_wrapBWplot <- function(flowSet,
                               channels,
                               box_ymin = 0.05,
                               box_lower = 0.25,
                               box_middle = 0.5,
                               box_upper = 0.75,
                               box_ymax = 0.95,
                               x_text_size = 2,
                               strip_text_size = 7,
                               legend_text_size = 0.5,
                               legend_title = "facetName"){

    marker_id <- match(channels, colnames(flowSet))

    if(any(is.na(marker_id))){
        stop(paste0("Can not find the specified channels:",
                    paste(channels[is.na(marker_id)], collapse = " "),
                    " in the flowSet data!"))
    }

    data <- NULL
    boxClassFactor <- NULL
    for(sample in names(flowSet@frames)){
        fcs <- flowSet@frames[[sample]]
        pd <- fcs@parameters@data
        exprs <- fcs@exprs[ ,marker_id, drop=FALSE]
        colnames(exprs) <- paste0(pd$name, "<", pd$desc,">")[marker_id]
        data <- rbind(data, exprs)
        boxClassFactor <- c(boxClassFactor, rep(sample, nrow(exprs)))
    }

    wrapBWplot(data = data,
               boxClassFactor = boxClassFactor,
               box_ymin = box_ymin,
               box_lower = box_lower,
               box_middle = box_middle,
               box_upper = box_upper,
               box_ymax = box_ymax,
               x_text_size = x_text_size,
               strip_text_size = strip_text_size,
               legend_text_size = legend_text_size,
               legend_title = legend_title)
}

#' @importFrom dplyr group_by summarize
#' @importFrom reshape melt
#' @import ggplot2
wrapBWplot <- function(data,
                       boxplotCols,
                       boxClassFactor,
                       box_ymin = 0.05,
                       box_lower = 0.25,
                       box_middle = 0.5,
                       box_upper = 0.75,
                       box_ymax = 0.95,
                       x_text_size = 2,
                       strip_text_size = 7,
                       legend_text_size = 0.5,
                       legend_title = "facetName"
){

    if(missing(boxplotCols)){
        boxplotCols <- colnames(data)
    }else if(any(!(boxplotCols %in% colnames(data)))){
        stop("Unmatch boxplotCols found:", paste(boxplotCols[!(boxplotCols %in% colnames(data))], collapse = " "))
    }

    if(missing(boxClassFactor)){
        warning("no boxClassFactor was provided!")
        boxClassFactor <- rep("stack", length = nrow(data))
    }else if(length(boxClassFactor) != nrow(data)){
        stop("Length of boxClassFactor unequal row number of input data")
    }

    cat(" Reshaping data ...\n")
    data <- data.frame(data[ ,boxplotCols, drop=FALSE], boxClassFactor = boxClassFactor, check.names = FALSE)
    ldata <- melt(data, id.vars = "boxClassFactor", variable.name = "boxName", value.name = "value")

    cat(" Pre-stat(box range) calculating ...\n")
    ldata_group <- dplyr::group_by(ldata, boxClassFactor, boxName) %>%
        dplyr::summarise(ymin=quantile(value, box_ymin),
                  lower=quantile(value, box_lower),
                  middle=quantile(value, box_middle),
                  upper=quantile(value, box_upper),
                  ymax=quantile(value, box_ymax))

    cat(" Plotting ...\n")
    wrapBWplot_theme <- theme(legend.position = "top",
                              legend.title = element_text(size = rel(0.6)),
                              legend.text = element_text(size = rel(legend_text_size)),
                              legend.background = element_rect(fill=NULL, size=0.8, linetype="solid"),
                              strip.text = element_text(size=strip_text_size, lineheight=0.1, hjust = 0.5, vjust = 0.5),
                              axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              strip.background=element_rect(fill = "grey90", colour = NA))

    gp <- ggplot(ldata_group, aes(x=boxClassFactor, ymin=ymin, lower=lower, middle=middle,
                                  upper=upper, ymax=ymax, fill = boxClassFactor)) +
        geom_boxplot(stat="identity") + facet_wrap(~ boxName, scales = "free") +
        guides(fill = guide_legend(title.hjust = 0.5, keywidth = 0.5, keyheight = 0.5,
                                   title.position = "top", title = legend_title)) +
        theme_bw() + wrapBWplot_theme

    gp
}
