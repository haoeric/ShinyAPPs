library(ggplot2)
library(cytofkit)
library(reshape2)
library(dplyr)
library(flowCore)

stackDenistyPlot <- function(data, densityCols, stackFactor,
                             kernel = c("gaussian", "epanechnikov", "rectangular",
                                        "triangular", "biweight",
                                        "cosine", "optcosine"),
                             bw = "nrd0", adjust = 1,
                             reomoveOutliers = FALSE,
                             stackRotation = 0,
                             stackSeperation = "auto",
                             x_text_size = 2,
                             strip_text_size = 7,
                             legend_text_size = 0.5,
                             legendRow = 1,
                             legend_title = "stackName"){

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

    densData <- .densityCal(data, kernel = kernel, bw = bw, adjust = adjust, reomoveOutliers = reomoveOutliers)
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
                                    legend.title = element_text(size = rel(1)),
                                    legend.text = element_text(size = rel(legend_text_size)),
                                    strip.text = element_text(size=strip_text_size, lineheight=1, hjust = 0.5, vjust = 0.5),
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
        guides(col = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE),
               fill = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE)) +
        theme_bw() + stackDensityPlot_theme

    gp
}



#' Internal density calculation function serves for \code{stackDenistyPlot}
#'
#' Output data frame with columns: stackName, x , y , densityName
.densityCal <- function(data, kernel, bw, adjust, reomoveOutliers = FALSE){
    cat("  Calculating Density for each stack column...\n")
    print(table(data$stackFactor))
    dataBystackFactor <- split(subset(data, select = -stackFactor), data$stackFactor)
    densityWrap <- function(d, ...){
        resOut <- NULL
        for(i in colnames(d)){
            x <- d[,i]
            if(reomoveOutliers){
                cat("  Remove outliers...\n")
                x_IQR <- IQR(x)
                x_lowLimit <- quantile(x, 0.25) - 1.5*x_IQR
                x_highLimit <- quantile(x, 0.75) + 1.5*x_IQR
                x <- x[x >= x_lowLimit && x <= x_highLimit]
            }
            dens <- density(x, ...)
            densOut <- data.frame(x=dens$x, y=dens$y, densityName = i)
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




dotPlot <- function(data, xlab, ylab, pointSize=1, freeScales = FALSE){
    colPalette <- colorRampPalette(c("blue", "turquoise", "green",
                                     "yellow", "orange", "red"))
    data <- as.data.frame(data)
    mnames <- c(xlab, ylab)
    data <- data[ ,mnames]
    colnames(data) <- c("col1", "col2")
    densCol <- densCols(data, colramp = colPalette)
    data$sample <- sub("_[0-9]*$", "", row.names(data))

    gp <- ggplot(data, aes_string(x="col1", y="col2")) +
        geom_point(colour=densCol, size = pointSize)
    if(freeScales){
        gp <- gp + facet_wrap(~sample, scales = "free")
    }else{
        gp <- gp + facet_wrap(~sample, scales = "fixed")
    }

    gp <- gp +
        xlab(mnames[1]) + ylab(mnames[2]) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))

    return(gp)
}



gaussNorm <- function(exprs,
                      samplesInfo,
                      channel.names,
                      channelRanges,
                      max.lms=2,
                      base.lms=NULL,
                      peak.density.thr=0.05,
                      peak.distance.thr=0.05){

    if(length(max.lms)==1){
        max.lms=rep(max.lms, times=length(channel.names))
    }
    if(length(max.lms)!=length(channel.names)){
        cat('Error: length of max.lms and channel.names doesn\'t match\n')
        return(NULL)
    }
    names(max.lms)=channel.names
    expr.list <- unique(samplesInfo)

    ## Extract landmarks
    cat("  Identify landmarks ...\n")
    lms.list=list()
    lms.list$original=matrix(vector("list"), length(channel.names), length(expr.list))
    lms.list$score=matrix(vector("list"), length(channel.names), length(expr.list))
    lms.list$filter=matrix(vector("list"), length(channel.names), length(expr.list))
    row.names(lms.list$original)=channel.names
    row.names(lms.list$score)=channel.names
    row.names(lms.list$filter)=channel.names
    colnames(lms.list$original)=expr.list
    colnames(lms.list$score)=expr.list
    colnames(lms.list$filter)=expr.list
    max.lms.num=0
    c=1

    for( p in channel.names){
        cat(paste0("Normalize channel - ", p, "...\n"))
        if(max.lms[c]!=0){
            j=1
            for(i in expr.list){
                data=exprs[samplesInfo == i, ]
                ## finding the landmarks of channel p of sample i.
                #data[indices$index[[i]],]<-NA
                lms=landmarker(data, p, max.lms[c])
                lms.list$original[p, j]=list(lms)
                ## returns the max.lms[c] top score landmarks.
                filtered=filter.lms(lms, data[,p], max.lms[c],  peak.density.thr, peak.distance.thr)
                lms.list$filter[p, j]=list(filtered$lms)
                lms.list$score[p, j]=list(filtered$score)
                j=j+1
            }
        }
        c=c+1
    }

    if(is.null(base.lms))
        base.lms=extract.base.landmarks(lms.list$filter, channel.names, max.lms)

    ## finds the best matching between landmarks and base landmarks.
    cat("  Register landmarkers ...\n")
    matched.lms=match.all.lms(lms.list, base.lms, channel.names, max.lms)
    colnames(matched.lms)=expr.list
    confidence=compute.confidence(matched.lms, base.lms)
    cat('  Aligning based on landmarks ...\n')
    for(i in expr.list){
        exprs[samplesInfo == i, ]=normalize.one.expr(exprs[samplesInfo == i, ],
                                                     base.lms, lms$filter[,i],
                                                     lms$original[,i], matched.lms[,i],
                                                     channel.names, max.lms)
    }

    return(exprs)
}


## shifts the data in such a way that the peak at matched.lms[i] is moved to matched.lms[i+1] for each i.
## base.lms, lms.lis and lms.original are passed for the debug purposes only.
normalize.one.expr <- function(data, base.lms, lms.list, lms.original, matched.lms, channel.names, max.lms, fname='', debug=FALSE){
    norm.data=data
    if(debug==TRUE){
        if(fname=='')
            x11()
        else
            png(filename=fname, bg="white",width=1800,height=1000)
        par(mfrow=c(1, length(which(max.lms!=0))))
    }
    cn=1
    i=0
    ## normalizing each channel.
    for (p in channel.names){
        i=i+1
        if(max.lms[i]==0){
            norm.data[,p]=data[,p]
        }
        else{
            ## registering the data of channel p given the matching landmarks.
            norm.data[,p]=register.channel(data[,p], matched.lms[[p]])
            ## drawing the pre- and post- normalization density curves.
            if(debug==TRUE){
                if(length(lms.list[[p]])==0){
                    A1=density(na.omit(data[,p]))
                    xlim=c(min(A1$x, na.rm=T), max(A1$x, na.rm=T))
                    ylim=c(min(A1$y, na.rm=T), max(A1$y, na.rm=T))
                    par(mfg=c(1, cn))
                    plot(A1, type="l", col= "black", xlim=xlim, ylim=ylim, xlab=p, main='No lms is found')
                }
                else{
                    A1=density(na.omit(data[,p]))
                    A2=density(na.omit(norm.data[,p]))
                    xlim=c(min(min(A1$x, na.rm=T), min(A2$x, na.rm=T)), max(max(A1$x, na.rm=T), max(A2$x, na.rm=T)))
                    ylim=c(min(min(A1$y, na.rm=T), min(A2$y, na.rm=T)), max(max(A1$y, na.rm=T), max(A2$y, na.rm=T)))
                    par(mfg=c(1, cn))
                    plot(A1, type="l", col= "blue", xlim=xlim, ylim=ylim, xlab=p, main='')
                    Lms=matched.lms[[p]][seq(1, length(matched.lms[[p]]), by=2)]
                    M.Lms=matched.lms[[p]][seq(2, length(matched.lms[[p]]), by=2)]
                    M.Lms.index=1
                    for(k in 1:(length(M.Lms)))
                        M.Lms.index[k]=which(base.lms[[p]]==M.Lms[k])
                    points(Lms, returny(A1, Lms), pch=19, col="red")
                    text(Lms, returny(A1, Lms), as.character(M.Lms.index), pos=3, col='red')
                    points(lms.original[[p]], returny(A1, lms.original[[p]]), pch=21, col="black") ##all the peaks
                    par(mfg=c(1, cn))
                    plot(A2, type="l", col="red", xlim=xlim, ylim=ylim, xlab=p)
                    points(M.Lms, returny(A2, M.Lms), pch=15, col="blue")
                    text(M.Lms, returny(A2, M.Lms), as.character(M.Lms.index), pos=3, col='blue')
                }
                cn=cn+1
            }
        }
    }
    if(debug){
        if(fname=='')
            readline()
        dev.off()
    }
    return (norm.data)

}


## returns the peaks (local maxima's) in the kernel density estimate of data.
landmarker <- function(data, channel.name, max.lms, span=3){
    d=data[,channel.name]
    A=density(na.omit(d))
    d=A$y
    pks=c()
    ## sliding a window of size span and returns locations where the middle point in the window is maximum.
    for( i in 1:(length(d)-span)){
        if(!is.na(d[i+span%/%2]) & (d[i+span%/%2]==max(d[c(i:(i+span))], na.rm=T)))
            if(!is.na(d[i]) & !is.na(d[i+span]) & d[i+span%/%2]!=d[i] & d[i+span%/%2]!=d[i+span])
                pks=append(pks, i+span%/%2)

    }
    return(A$x[pks])
}


extract.base.landmarks <- function(lms, channel.names, max.lms){
    lms.list=list()
    max.lms.num=0
    for( p in channel.names){
        if(max.lms[p]!=0){
            ## if the total number of landmarks for channel p is less than max.lms[p] return them as base landmarks.
            if(length(unlist(lms[p,])) <= max.lms[p]){
                lms.list[[p]]=sort(unlist(lms[p,]))
            }
            else{

                if (max.lms[p]==1){
                    lms.list[[p]]=median(unlist(lms[p,]))
                }
                else {
                    ## first identify samples that have exactly max.lms[p] landmarks on their channel p. These landmarks are labeled from 1 to max.lms
                    ## the landmarks samples with less than max.lms[p] landmarks are stored in short.lms vector.
                    short.lms=vector()
                    lms.class=list()
                    for(jj in 1:max.lms[p])
                        lms.class[[jj]]=vector()
                    for (ii in 1:length(lms[p,])){
                        lms.p.ii=lms[p,ii][[1]]
                        if(length(lms.p.ii)==max.lms[p]){
                            for (jj in 1:max.lms[p])
                                lms.class[[jj]]=append(lms.class[[jj]], lms.p.ii[jj])
                        }
                        else{
                            short.lms=append(short.lms, lms.p.ii)
                        }
                    }

                    if(length(lms.class[[1]])==0){
                        cat('No curve with ', max.lms[p], ' landmarks was found for channel ', p, '. Decrease max.lms[',p,'] and try again.\n')
                        stop()
                    }
                    mean.lms.class=0
                    for(jj in 1:max.lms[p]){
                        mean.lms.class[jj]=mean(lms.class[[jj]])
                    }

                    ## assigning short.lms landmarks to the class of labeled landmarks wich are closer to.
                    if(length(short.lms)!=0){
                        for(jj in 1:length(short.lms)){
                            kk=which(abs(mean.lms.class-short.lms[jj])==min(abs(mean.lms.class-short.lms[jj])))
                            kk=kk[1]
                            lms.class[[kk]]=append(lms.class[[kk]], short.lms[jj])
                        }
                    }
                    lms.class.len=0
                    lms.class.med=0
                    for(jj in 1:max.lms[p]){
                        lms.class.len[jj]=length(lms.class[[jj]])
                        lms.class.med[jj]=median(lms.class[[jj]], na.rm=T)
                    }
                    s=sd(lms.class.len, na.rm=T)
                    m=mean(lms.class.len, na.rm=T)
                    ll=which(m-lms.class.len > 2*s)
                    ## if a class of landmarks has too few landmarks in it just ignore this class.
                    if(length(ll)!=0){
                        cat('warning: fewer landmark classes found in channel ', p, '\n')
                        tmp.lms=list()
                        kk=1
                        for(jj in (1:max.lms[p])[-ll]){
                            tmp.lms[[kk]]=lms.class[[jj]]
                            kk=kk+1
                        }
                        lms.class=tmp.lms

                    }

                    ## returning the median of each class as the base landmark.
                    for(jj in 1:length(lms.class)){
                        lms.class.med[jj]=median(lms.class[[jj]], na.rm=T)
                    }

                    lms.list[[p]]=lms.class.med

                }
            }
        }
    }
    return(lms.list)
}

## returns the max.lms top score landmarks.
## the score of a landmarks is a funciton of its sharpness and density value
filter.lms <- function(lms, data, max.lms, peak.density.thr, peak.distance.thr){
    filtered=list()
    if(length(lms) == 0){
        filtered$lms=vector()
        filtered$score=vector()
        return(filtered)
    }
    filtered$score=score.lms(lms, data, max.lms, peak.density.thr, peak.distance.thr)
    lms.score=data.frame(score=filtered$score, ind=c(1:length(lms)))
    lms.score=lms.score[do.call(order, c(lms.score["score"], decreasing=T)), ]
    ind=which(lms.score$score>0)
    if(length(ind)==0){
        filtered$lms=vector()
        filtered$score=vector()
        return(filtered)
    }
    lms.score.ind=lms.score$ind[ind]
    if(length(lms.score.ind)<max.lms)
        filtered$lms=sort(lms[lms.score.ind], decreasing=F)
    else
        filtered$lms=sort(lms[lms.score.ind[c(1:max.lms)]], decreasing=F)
    return(filtered)
}


## assigns a score to each landmark. the score of a landmarks is a funciton of its sharpness and density value.
## the peaks with density value less than peak.density.thr*maximum.peak.density are discarded.
## of the peaks with distance less than peak.distance.thr*range.data only one is considered.

score.lms <- function(lms, data, max.lms, peak.density.thr, peak.distance.thr){
    bw=64
    score=vector()
    height.cutoff=peak.density.thr
    if(length(lms) == 0)
        return(score)
    A=density(na.omit(data))
    bw=min(64, length(A$x)/10)
    lms.max.height=max(returny(A, lms), na.rm=T)
    MIN.LMS.DIST=(max(A$x, na.rm=T)-min(A$x, na.rm=T))*peak.distance.thr
    last.lms=-1
    last.lms.i=-1
    last.lms.score=0
    for(i in 1:length(lms)){
        lms.ind=which(na.omit(abs(A$x-lms[i]))==min(na.omit(abs(A$x-lms[i])) ) )
        ind=(max(lms.ind-bw%/%2, 1)):(min(lms.ind+bw%/%2, length(A$x)))
        if(length(ind)==0)
            ind=1

        if(A$y[lms.ind] <   height.cutoff*lms.max.height){
            w=0
        }
        else{
            ## computing the sharpness
            w=A$y[lms.ind]-A$y[ind]
            w[which(w<0)]=3*w[which(w<0)]
        }
        ## computing final score
        score[i]=sum(w, na.rm=T)*A$y[lms.ind]
        if(score[i]<0)
            score[i]=0
        if(last.lms<0){
            last.lms=lms[i]
            last.lms.i=i

        }
        else{
            ##If two lms's are very close only choose one with the better score.
            if(lms[i]-last.lms < MIN.LMS.DIST){
                if(score[i]>score[last.lms.i]){
                    last.lms=lms[i]
                    score[last.lms.i]=0
                    last.lms.i=i
                }
                else{
                    score[i]=0
                }
            }
            else{
                last.lms=lms[i]
                last.lms.i=i
            }
        }
    }
    return(score)
}


match.all.lms <- function(lms, base.lms, channel.names, max.lms){
    n=length(lms$filter[channel.names[1],]) ##number of samples.
    matched.lms=matrix(vector("list"), length(channel.names), n)
    row.names(matched.lms)=channel.names
    lms.class=list()
    lms.class.median=list()
    for(p in channel.names){
        lms.class[[p]]=list()
        lms.class.median[[p]]=list()
    }
    for (p in channel.names){
        if(max.lms[[p]]==0)
            next
        for(i in 1:n){
            matched.lms[p, i][[1]]=best.match(lms$original[p, i][[1]], lms$score[p, i][[1]], base.lms[[p]], max.lms[[p]])
        }
    }
    return(matched.lms)
}

best.match <- function(lms, lms.score, base.lms, max.lms){
    comb=choose.nk(max(length(lms), length(base.lms)), min(length(lms), length(base.lms)))
    sc=list()
    k=1
    max.s=-1
    if(length(lms)<length(base.lms)){
        for(i in 1:(dim(comb)[1])){
            if(length(which(comb[i,]==0))!=0)
                c=comb[i,][-which(comb[i,]==0)]
            else
                c=comb[i,]
            d=combinations.itr(length(comb[i,]), length(c))
            for(j in 1:(dim(d)[1])){
                s=match.score(lms[d[j,]], base.lms[c], lms.score[d[j,]])
                if(max.s<s){
                    lms.index=d[j,]
                    base.lms.index=c
                    max.s=s
                }
                sc[[k]]=list(score=s, lms.index=d[j,], base.lms.index=c)
                k=k+1
            }
        }

    }
    else{
        for(i in 1:(dim(comb)[1])){
            if(length(which(comb[i,]==0))!=0)
                c=comb[i,][-which(comb[i,]==0)]
            else
                c=comb[i,]
            d=combinations.itr(length(comb[i,]), length(c))
            for(j in 1:(dim(d)[1])){
                s=match.score(lms[c], base.lms[d[j,]], lms.score[c])
                if(max.s<s){
                    lms.index=c
                    base.lms.index=d[j,]
                    max.s=s
                }

                k=k+1
            }
        }
    }
    res=1:(2*length(lms.index))
    res[seq(1,2*length(lms.index), by=2)]=lms[lms.index]
    res[seq(2,2*length(lms.index), by=2)]=base.lms[base.lms.index]
    return(res)
}

## defines the score of a matching between landmarks (lms) and base landmarks (base.lms)
match.score <- function(lms, base.lms, lms.score){
    d=abs(lms-base.lms)
    if(length(which(d==0))!=0)
        d[which(d==0)]=0.01
    return(sum(1/d*lms.score*lms.score))
}

choose.nk <- function(n, k){
    v=matrix(0, ncol=k, nrow=0)
    for(j in 1:k){
        c=combinations.itr(n, j)
        if(j<k)
            c=cbind(c, matrix(0, nrow=dim(c)[1], ncol=k-j))
        v=rbind(v, c)
    }
    return(v)
}

register.channel <- function(data, matched.lms){
    if(length(matched.lms)==0){
        cat('*')
        return (data)
    }
    s=m=shift=vector()
    lms=vector()
    for(i in seq(1,length(matched.lms), by=2)){
        shift=append(shift, matched.lms[i+1]-matched.lms[i])
        lms=append(lms, matched.lms[i])
        s=append(s, sd(na.omit(data)))
    }
    r.data=register.function(data, s, lms, shift)
    return(r.data)
}

returny <- function(A, X){
    y=vector()
    i=1;
    for( x in X){
        y[i]=A$y[which(abs(A$x-x)==min(abs(A$x-x)))[1]]
        i=i+1
    }
    return (y)
}


combinations.itr <- function(n, k){
    v=1
    res=NULL
    while(length(v)!=0){
        if(v[length(v)]>n){
            v=v[-length(v)]
            v[length(v)]=v[length(v)]+1
        }
        else if(length(v)==k){

            res=rbind(res, v)
            v[length(v)]=v[length(v)]+1
        }
        else{
            v=append(v, v[length(v)]+1)
        }
    }
    return(res)
}


## computes the confidence value based on the matched landmarks.
compute.confidence <- function(matched.lms, base.lms){
    confidence=rep(1, times=length(matched.lms[1,]))
    # for each channel c.
    for(c in 1:length(base.lms)){
        for(l in 1:length(base.lms[[c]])){
            lms.list=0
            for ( i in 1:length(matched.lms[1,])){
                ind=which(matched.lms[c,i][[1]]==base.lms[[c]][l])
                if(length(ind)==0){
                    lms.list[i]=NA
                }
                else{
                    if(ind[1] %% 2 != 0)
                        ind[1]=ind[1]+1
                    lms.list[i]=matched.lms[c,i][[1]][ind[1]-1]
                }
            }
            d=density(na.omit(lms.list))
            den.lms=0
            for(i in 1:length(lms.list)){
                if(is.na(lms.list[i]))
                    den.lms[i]=NA
                else
                    den.lms[i]=d$y[which(abs(d$x-lms.list[i])==min(abs(d$x-lms.list[i])))]
            }
            m.den.lms=mean(den.lms, na.rm=T)
            ind1=which(den.lms>m.den.lms/4 & den.lms<m.den.lms/3)
            ind2=which(den.lms<m.den.lms/4)
            if(length(ind1)!=0)
                confidence[ind1]=confidence[ind1]*0.8
            if(length(ind2)!=0)
                confidence[ind2]=confidence[ind2]*0.7
            if(length(which(is.na(lms.list)))!=0)
                confidence[which(is.na(lms.list))]=confidence[which(is.na(lms.list))]*0.6
        }
    }
    confidence
}



register.function <- function(data, s, m, shift){
    sum=0
    if(length(m)==1){
        return(data+shift)
    }
    if(length(m)==2){
        sh=(shift[1]-shift[2])
        data=data+gau(data, s[1], m[1])*(sh/2)
        data=data-gau(data, s[2], m[2])*(sh/2)
        return(data+shift[1]-sh/2)
    }
    max.shift=which(abs(shift)==max(abs(shift)))[1]
    if(shift[max.shift]>0){
        sh=(shift[max.shift]-(shift[max.shift]-min(shift[-max.shift]))/2)
    }
    else{
        sh=(shift[max.shift]-(shift[max.shift]-max(shift[-max.shift]))/2)
    }
    data=data+sh
    shift=shift-sh

    for(i in 1:length(m))
        data=data+gau(data, s[i], m[i])*shift[i]
    return (data)
}

## gaussian function used in shifting the data.
gau <- function(d, s, m){
    return(2.7182^(-((d-m)^2)/(2*s*s)))
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



## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}



