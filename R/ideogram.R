
.rwb <- c("#0000FF","#0B0BFF","#1515FF","#2020FF","#2B2BFF","#3535FF","#4040FF"
,"#4A4AFF","#5555FF","#6060FF","#6A6AFF","#7575FF","#8080FF","#8A8AFF"
,"#9595FF","#9F9FFF","#AAAAFF","#B5B5FF","#BFBFFF","#CACAFF","#D4D4FF"
,"#DFDFFF","#EAEAFF","#F4F4FF","#FFFFFF","#FFFFFF","#FFF4F4","#FFEAEA"
,"#FFDFDF","#FFD5D5","#FFCACA","#FFBFBF","#FFB5B5","#FFAAAA","#FF9F9F"
,"#FF9595","#FF8A8A","#FF8080","#FF7575","#FF6A6A","#FF6060","#FF5555"
,"#FF4A4A","#FF4040","#FF3535","#FF2B2B","#FF2020","#FF1515","#FF0B0B"
,"#FF0000")



setClass("cytoband",representation(stain="character",
                                   band="character",
                                   start="numeric",
                                   end="numeric",
                                   length="integer"))

buildChromLocation.2 <- function (dataPkg,major=NULL) 
{
    CHRLOC2chromLoc <- function(chrEnv) {
        chrLocs <- contents(chrEnv)
        chrLens <- sapply(chrLocs, length)
        multis <- split(chrLens, factor(chrLens))
        singleNames <- names(multis$"1")
        singleLocs <- chrLocs[singleNames]
        chromNames <- unlist(sapply(singleLocs, function(z) {
            if (is.na(z)) 
                z
            else names(z)
        }))
        chromNames <- factor(chromNames)
        a <- split(singleLocs, chromNames)
        chrLocList <- lapply(a, function(x) {
            g <- unlist(lapply(x, function(z) {
                names(z) <- NULL
                z
            }))
            g
        })
        if (length(multis) > 1) {
            for (i in 2:length(multis)) {
                curNames <- names(multis[[i]])
                curLocs <- chrLocs[curNames]
                for (j in 1:length(curLocs)) {
                  curGene <- curLocs[[j]]
                  curGeneChroms <- names(curGene)
                  names(curGene) <- rep(curNames[j], length(curGene))
                  for (k in 1:length(curGene)) chrLocList[[curGeneChroms[k]]] <- c(chrLocList[[curGeneChroms[k]]], 
                    curGene[k])
                }
            }
        }
        chrLocList
    }
    if (!require(dataPkg, character.only = TRUE)) 
        stop(paste("Package:", dataPkg, "is not available"))
    pEnv <- paste("package", dataPkg, sep = ":")
    chrLocList <- CHRLOC2chromLoc(get(paste(dataPkg, "CHRLOC", 
        sep = ""), pos = pEnv))

    species <- tolower(substr(get(paste(dataPkg, "ORGANISM", sep = "")), 1, 1))
    cytoEnv <- switch(species, h = get("Hs.cytoband", "package:ideogram"),
	r = get("Rn.cytoband", "package:ideogram"), m = get("Mm.cytoband", 
	"package:ideogram"), NULL)
    if (is.null(cytoEnv)) stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")    
    
    
    chrLocListNames <- names(chrLocList)
    badNames <- NULL
    for(i in 1:length(chrLocListNames)) if(inherits(try(get(chrLocListNames[i], pos = cytoEnv),silent=T), "try-error")) badNames <- c(badNames,i)
    chrLocListNames <- chrLocListNames[-badNames]
    
    
    if("arms" %in% major){

	chrLocList2 <- list()
	armList <- NULL
	for(i in chrLocListNames){
		cyto <- try(get(i, pos = cytoEnv),silent=T)
		if (inherits(cyto, "try-error")) next()
		breakPoint <- cyto@start[min(grep("q",cyto@band))]
		ps <- chrLocList[[i]][chrLocList[[i]] < breakPoint]
		qs <- chrLocList[[i]][chrLocList[[i]] > breakPoint]
		l <- length(chrLocList2)
		names <- names(chrLocList2)
		chrLocList2 <- c(chrLocList2,list(ps))
		chrLocList2 <- c(chrLocList2,list(qs))
		names(chrLocList2) <- c(names,paste(i,"p",sep=""),paste(i,"q",sep=""))
		armList <- c(armList,paste(i,"p",sep=""),paste(i,"q",sep=""))
	}
	chrLocList <- c(chrLocList,chrLocList2,armList=list(armList))
	
    }
    
    if("bands" %in% major){
	chrLocList2 <- list()
	bandList <- NULL
	for(i in chrLocListNames){
		cyto <- try(get(i, pos = cytoEnv),silent=T)
		if (inherits(cyto, "try-error")) next()
		bands <- try(gsub("\\..*", "", cyto@band),silent=T)
		if (inherits(bands, "try-error")) next()	
		ix <- !duplicated(bands)
		bands <- bands[ix]
		bands2 <- paste(i,bands,sep="")
		start <- cyto@start[ix]
		names <- names(chrLocList2)
		for(j in 1:length(bands)){
			if(j < length(bands)){
				toAdd <- chrLocList[[i]][(start[j] < chrLocList[[i]]) & (start[j+1] > chrLocList[[i]])]
			} else if(j==length(bands)){
				toAdd <- chrLocList[[i]][(start[j] < chrLocList[[i]])]
			}
			chrLocList2 <- c(chrLocList2,list(toAdd))
		}
		names(chrLocList2) <- c(names,bands2)
		bandList <- c(bandList,bands2)
	}
	chrLocList <- c(chrLocList,chrLocList2,bandList=list(bandList))
    }
    

    if("mb" %in% major){
	chrLocList2 <- list()
	mbNamesList <- NULL
	for(i in chrLocListNames){
		cyto <- try(get(i, pos = cytoEnv),silent=T)
		if (inherits(cyto, "try-error")) next()
		
		megaNames <- paste(i,as.character(as.integer(seq(0,floor(max(chrLocList[[i]])/1000000)*1000000,1000000))),sep="-")
		megaList <- split(chrLocList[[i]],seq(0,floor(max(chrLocList[[i]])/1000000)*1000000,1000000))
		names <- names(chrLocList2)
		chrLocList2 <- c(chrLocList2,megaList)
		names(chrLocList2) <- c(names,megaNames)
		mbNamesList <- c(mbNamesList,megaNames)
	}
	chrLocList <- c(chrLocList,chrLocList2,mbList=list(mbNamesList))
    }    

    newCC <- new("chromLocation", organism = get(paste(dataPkg, 
        "ORGANISM", sep = ""), pos = pEnv), dataSource = dataPkg, 
        chromLocs = chrLocList, chromInfo = get(paste(dataPkg, 
            "CHRLENGTHS", sep = ""), pos = pEnv), probesToChrom = get(paste(dataPkg, 
            "CHR", sep = ""), pos = pEnv), geneSymbols = get(paste(dataPkg, 
            "SYMBOL", sep = ""), pos = pEnv))
    return(newCC)
}




.combine <- function(x,simplify=TRUE) {
  return(list(x))
}

.parse.chr <- function(x) {
  ids <- try(eval(parse(text=x)),silent=TRUE)
  if(inherits(ids,"try-error")) {
    return(x)    
  } else {
    return(ids)
  }
}

.pick.element <- function(x,index) {
  return(x[index])
}

.naMean <- function(x) {
  mean(x,na.rm=TRUE)
}

.usedChromExprs <- function(exprs,genome,chr,aggrfun=NULL)  {
  if(!is.matrix(exprs))
    stop("requires a matrix")
  locs <- try(chromLocs(genome)[[chr]]) ## simple usedChromGenes
  if (inherits(locs,"try-error") || is.null(locs)) {
    cat("Warning: no expression values in region ",chr,"\n")
    return(NULL)
  }
  locs <- locs[!duplicated(names(locs))] ## get a single location for each gene, just pick the first on
  ids <- rownames(exprs)

  multi <- try(grep('c\\(.*\\)',ids))
  if(length(multi)) {  ## if we aggregated before, the locations may be complex
    sids <- ids
    ll.ids <- lapply(ids[multi],.parse.chr)
    ## sids[multi] <- lapply(ll.ids,.pick.element,1)
    sids[multi] <- unlist(lapply(ll.ids,.pick.element,1))
  } else {
    sids <- ids
  }
  keep <- sids %in% names(locs)
  exprs <- as.matrix(exprs[keep,])
  sids <- sids[keep]
  ## sometimes the ann env contains more data then the eset...make sure they agree
  locs <- locs[match(sids,names(locs))]
  locs <- abs(locs)
  ix <- order(locs)
  locs <- locs[ix]
  exprs <- as.matrix(exprs[ix,])
  ids <- rownames(exprs)
  rownames(exprs) <- locs
  if (!is.null(aggrfun) && sum(duplicated(locs)) > 0) {
    if(!is.function(aggrfun)) 
	stop("requires a summary function")
    begin.nrow <- nrow(exprs)
    f <- factor(locs)
    ids <- tapply(ids,f,.combine)
    exprs <- aggregate(exprs,by=list(f),aggrfun)
    locs <- levels(f)
    exprs <- as.matrix(exprs[,2:ncol(exprs)])
    rownames(exprs) <- locs
    sids <- unlist(lapply(ids,.pick.element,1))
    warning(begin.nrow - nrow(exprs), " genes map to the same location on chromosome ",chr," and were summarized.","\n")    
  } else {
    sids <- ids
  }
  return(list(exprs=as.matrix(exprs),locs=as.numeric(locs),geneIDs=ids,simpleIDs=sids))
}

ideogram <- function(data,genome,chr=NULL,organism=NULL,method=c("plot","matplot","image"),margin=c("ticks","ideogram"),grid.col=c("red","grey"),grid.lty=c(1,2),widths=c(1,2),relative=FALSE,dlim=NULL,main=NA,xlab=NA,ylab=NA,cex.axis=.7,na.color=par("bg"),...){
  method <- match.arg(method)
  margin <- match.arg(margin)
  
  ##Initial setup
  if (is.null(chr))
    stop("No chromosome chosen","\n")
  if(!is.character(chr))
    chr <- as.character(chr)
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  cytoEnv <- NULL
  cytoEnv <- switch(organism,
                    "h"=get("Hs.cytoband","package:ideogram"),
                    "r"=get("Rn.cytoband","package:ideogram"),
                    "m"=get("Mm.cytoband","package:ideogram"),
                    NULL)
                    
  if(is.null(cytoEnv))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")

  cyto <- try(get(paste(chr),pos=cytoEnv))
  if(inherits(cyto,"try-error"))
    stop("Chromosome ",chr," is not recognized")

  if(is.na(main)) main <- chr
  
  if(class(data) == "exprSet") data <- data@exprs
	  
  if(method == "plot" & !is.vector(data)) {
	warning(sQuote("data")," needs to be a vector for this plot method: resetting to image")
	method <- "image"
  }
  
  if(method =="image" & !is.matrix(data)) {
	warning(sQuote("data")," needs to be a matrix for this plot method: resetting to plot")
	method <- "plot"
  }
  
  chr.size <- as.integer(max(cyto@end,na.rm=TRUE))
  if(is.null(chr.size))
    stop("chromosome ",chr," does not contain size information")

   chrList <- switch(method,
                    "plot"= .usedChromExprs(as.matrix(data),genome,chr,NULL),
                    "matplot"=.usedChromExprs(as.matrix(data),genome,chr,NULL),
                    "image"=.usedChromExprs(as.matrix(data),genome,chr,.naMean))

  ids <- chrList$geneIDs
  locs <- chrList$locs
  z <- chrList$exprs
  
  if(method == "image" & is.null(dlim)) {
    dlim <- range(z,na.rm=TRUE,finite=TRUE)
    cat("Warning: zlim used is: ",dlim," Please see ",sQuote("dlim")," to adjust.\n")
  }
  if(!is.null(dlim)) {
    z[z < dlim[1]] <- dlim[1]
    z[z > dlim[2]] <- dlim[2]
  }

  y <- locs-chr.size
  ylim <- c(0,-chr.size)
  args <- list(...)

  col.axis <- par("col.axis")
  font.axis <- par("font.axis")
  if(!is.null(args$col.axis)) col.axis <- args$col.axis
  if(!is.null(args$font.axis)) font.axis <- args$font.axis    
  
  ## Set the size of the chromosome based on the largest
  if(relative){
    allChr <- ls(pos=cytoEnv)
    allCyto <- mget(allChr,cytoEnv)
    allSize <- sapply(allCyto,function(x) max(x@end))
    maxLength <- max(allSize)
    inLen <- par()$pin[2]
    a <- par()$mai[1]
    b <- par()$mai[2]
    c <- par()$mai[3]
    d <- par()$mai[4]
    cur <- max(get(chr,pos=cytoEnv)@end)
    mod <- inLen*(1-(cur/maxLength))
    tempVec <- c(a,b,c+(mod),d)
    tempPar <- par()$mai
    par(mai=tempVec)
  }
  
  if(margin=="ideogram"){
    ##Draw Ideogram
    oldMai <- par()$mai
    op <- par(no.readonly=TRUE)
    leftOver <- oldMai[2]/7
    newMai <- c(oldMai[1],leftOver,oldMai[3],leftOver + par()$pin[1] + par()$mai[4])
    par(mai=newMai)
    ann <- c("acen","gvar","stalk","gneg","gpos25","gpos50","gpos75","gpos100")
    bColor <- c("red","grey","blue",grey(4:0/4))  
    new.colors <- rev(cyto@stain)
    ord <- match((new.colors),ann)
    bColor <- bColor[ord]
    bands <- (cyto@end-cyto@start)
    barplot(matrix(rev(bands),ncol=1),border="black",col=bColor,axes=FALSE)
    par(mai=oldMai,new=TRUE)
  }

  ##Points plot
  if(method=="plot"){
    z <- as.vector(z)
    if(!is.null(args$col)) {
      if(length(args$cols) != 1) {
        if(length(args$col) != length(data)) {
          warning("length of ",sQuote("col")," should equal length of ",sQuote("data"),". May produce unexpected results\n")
        } else {
          ix <- match(ids,names(data))
          args$col <- args$col[ix]
        }
      }
    }
    if(!is.null(args$pch)) {
      if(length(args$pch) != 1) {
        if(length(args$pch) != length(data)) {
          warning("length of ",sQuote("pch")," should equal length of ",sQuote("data"),". May produce unexpected results\n")
        } else {
          ix <- match(ids,names(data))
          args$pch <- args$pch[ix]
        }
      }
    }
    def.args <- list(y=y,x=z,ylim=ylim,axes=FALSE,main=main,xlab=xlab,ylab=ylab)
    args <- c(def.args,args)

    do.call("plot",args)
    axis(1,cex.axis=cex.axis,font.axis=font.axis,col.axis=col.axis)
##    plot(y=y,x=z,ylim=ylim,axes=FALSE,main=main,xlab=xlab,ylab=ylab,col=col,pch=pch, ...)
  }
  if(method=="matplot") {
    matplot(y=y,x=z,ylim=ylim,axes=FALSE,main=main,xlab=xlab,ylab=ylab,...)
    axis(1,cex.axis=cex.axis,font.axis=font.axis,col.axis=col.axis)
  }
  if(method=="image") {
    graph.dlim <- seq(from=dlim[1],to=dlim[2],length=nrow(z))
    z <- cbind(z,graph.dlim)
    offset <- dim(z)[2]-1
    offset <- (1/offset)
    offset <- offset/2
	  
	qq <- duplicated(y)
	y <- y[!qq]
	z <- z[!qq,]
	y <- c(y,chr.size)
	  
    image(y=y,z=t(z),axes=FALSE,xlim=c(0-offset,1-offset),ylim=ylim,zlim=dlim,main=main,xlab=xlab,ylab=ylab,...)
    if(!is.na(na.color) & any(is.na(y))) {
	na.z <- ifelse(is.na(z),1,NA)
	    str(z)
	    str(na.z)
	    cat("asdf")
	try(image(y=y,z=t(na.z),axes=FALSE,xlim=c(0-offset,1-offset),ylim=ylim,col=na.color,add=T),silent=T)
    }
	  

  }
  if(margin=="ticks"){
    ##Draw tick marks
    lab <-  cyto@band
    lab2 <-  gsub("\\..*","",lab)
    dupes <- duplicated(lab2)
    lab2[dupes] <- ""
    tickLocs <- -(chr.size-cyto@start)
    tickLocs2 <- tickLocs
    tickLocs[dupes] <- NA
    tickLocs2[!dupes] <- NA
    ##gray ticks
    try(axis(2, at = tickLocs2,labels = FALSE,las=2,col="gray"),silent=TRUE)
    ##black ticks
    axis(2, at = c(tickLocs,0) ,labels = FALSE,col="black")
    ##fix up label locations
    labLocs <- tickLocs[!is.na(tickLocs)]
    for(i in 1:(length(labLocs)-1))
      labLocs[i] <- labLocs[i] + (labLocs[i+1]-labLocs[i])/2
    labLocs[length(labLocs)] <- labLocs[length(labLocs)] - (labLocs[i] / 2)
    ##draw in lables
    axis(2, at = c(labLocs) ,tick=FALSE, labels = lab2[!is.na(tickLocs)],las=2,line=NA,cex.axis=cex.axis,font.axis=font.axis,col.axis=col.axis)
    ##re-draw black line on axis
    axis(2,labels=FALSE,col="black",tick=FALSE,cex.axis=cex.axis,font.axis=font.axis,col.axis=col.axis)
    ##draw dashed lines @ main bands
    if(!is.na(grid.lty[2])) {
      #str(lab2)
	for(i in 1:length(cyto@start)) if(lab2[i] != "") abline(h=-(chr.size-cyto@start[i]),col=grid.col[2],lty=grid.lty[2])

    }
    ##draw in centromere
    if(!is.na(grid.lty[1])) {
      for(i in which(cyto@stain=="acen")) {
        abline(h=-(chr.size-cyto@start[i]),col=grid.col[1],lty=grid.lty[1])
      }
    }
  }
  if(relative) par(mai=tempPar)
  return(invisible(list(x=z,y=locs-chr.size,labels=ids)))
}

ideograb <- function(ideo,show.box=TRUE,brush=NULL,...){
  x <- ideo$x
  y <- ideo$y
  ids <- ideo$labels
  cat("Please click on two points to define a retangular region.","\n")
  first <- locator(1)
  second <- locator(1)
  
  if(first$x > second$x) {
    x2 <- first$x
    x1 <- second$x
  }else{
    x2 <- second$x
    x1 <- first$x
  }

  between <- (x > x1) & (x < x2)
  
  if(first$y > second$y) {
    y1 <- first$y
    y2 <- second$y
  }else{
    y1 <- second$y
    y2 <- first$y
  }
  between2 <- (y > y2) & (y < y1)
  
  if(!is.null(brush)) try(points(x=x[(between & between2)],y=y[(between & between2)],col=brush,...))
  if(show.box){
    x <- vector(length=5,mode="numeric")
    x[1] <- first$x
    x[2] <- second$x
    x[3] <- second$x
    x[4] <- first$x
    x[5] <- first$x
    y <- vector(length=5,mode="numeric")
    y[1] <- first$y
    y[2] <- first$y
    y[3] <- second$y
    y[4] <- second$y
    y[5] <- first$y
    lines(list(x=x,y=y))
  }
  a <- ids[(between & between2)]
  
  return(a[!is.na(a)])
}

mideogram <- function(data,genome,chr=NULL,organism=NULL,method=c("plot","matplot","image"),margin=c("ticks","ideogram"),grid.col=c("red","grey"),grid.lty=c(1,2),widths=c(1,2),relative=TRUE,dlim=NULL,main=NA,xlab=NA,ylab=NA,cex.axis=.7,...){
  op <- par(no.readonly = TRUE)
  layout(rbind(c(1:8),c(9:16),c(17:24)))
  par(mai= c(par()$mai[1]*.5, par()$mai[2]*.6, par()$mai[3]*.4, par()$mai[4]*.6))		
  
  if(is.null(dlim)) dlim <- range(data,na.rm=TRUE)
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  chroms <- NULL
  chroms <- switch(organism,
                   "h"=c(1:22,"X","Y"),
                   "r"=get("Rn.cytoband","package:ideogram"),
                   "m"=get("Mm.cytoband","package:ideogram"),
                   NULL)
  
  if(is.null(chroms))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")
  
  for(i in chroms) {
    try(  ideogram(data,genome,i,organism=organism,method=method,margin=margin,grid.col=grid.col,grid.lty=grid.lty,widths=widths,relative=relative,dlim=dlim,main=main,xlab=xlab,ylab=ylab,cex.axis=cex.axis,...))
  }
  
  if(is.null(dlim)) dlim <- range(data,na.rm=TRUE,finite=TRUE)
  try(for(i in chr) ideogram(data=data,genome=genome,chr=i,method=method,widths=widths,dlim=dlim,margin=margin,grid.col=grid.col,grid.lty=grid.lty,relative=relative,...))
  par(op)
}
