\name{idiogram}
\alias{idiogram}
\alias{midiogram}
\title{ Plotting of Genomic Data}
\description{
  Function for plotting genomic data along with corresponding
  cytogenetic banding information
}
\usage{
idiogram(data, genome, chr=NULL, organism=NULL,
method=c("plot","matplot","image"), margin=c("ticks","idiogram"),
grid.col=c("red","grey"), grid.lty=c(1,2), widths=c(1,2),
relative=FALSE, dlim=NULL, main=NA, xlab=NA,
ylab=NA,cex.axis=.7,na.color = par("bg"),cen.color="red",...)
}

\arguments{
  \item{data}{ a vector or matrix of numeric data to plot. The names/rownames
    attribute \emph{needs} to contain corresponding gene identifiers}
  \item{genome}{ a chromLocation object associated with the specified
  data set. See below for details. }
  \item{chr}{ which chromosome to plot }
  \item{organism}{ if NULL, determination of the host organism will be
  retrieved from the \code{organism} slot of the chromLocation
  object. Otherwise \code{"h"}, \code{"r"}, or \code{"m"} can be used to
  specify \bold{h}uman, \bold{r}at, or \bold{m}ouse chromosome
  information}
  \item{method}{ plotting method }
  \item{margin}{ type of banding information to display in the plot
    margin}
  \item{grid.col}{ a two element vector specifying the centromere and
    band grid colors.}
  \item{grid.lty}{ a two element vector specifying the centromere and
    band grid line type.}
  \item{widths}{ a two element vector specifying the relative width of
    the margin idiogram two the adjacent graph. This option is currently
    ignored.}
  \item{relative}{ If \code{relative} is TRUE, the vertical height of
  the plot is scaled relative to the size of largest chromosome.}
  \item{dlim}{ a two element vector specifying the minimum and maximum
    values for \code{x}.  Data in \code{x} that exceed the min and max
  limits will be set to the min/max limit values before plotting.}
  \item{main}{ an overall title for the plot. Defaults to the chromosome
    name.}
  \item{xlab}{ a title for the x axis.}
  \item{ylab}{ a title for the y axis.}
  \item{cex.axis}{ the magnification to be used for axis annotation
    relative to the current.}
  \item{na.color}{ color to be used for NA values, defaults to: par("bg")}
  \item{cen.color}{ color to be used for the centromere when margin="idiogram", defaults to: "red"}  
  \item{\dots}{ additional graphical parameters can be given as
    arguments.}
}
\details{
  This function displays cytogenetic banding information in the plot
  margin and calls a secondary plotting function to display associated
  data at the same relative position.  Cytogenetic data for human, mouse,
  and rat genomes are currently included.

  The data is arranged by associating gene identifiers to genomic
  location using a \code{chromLoc} annotation object built using the
  \code{buildChromLocatio}n function from the annotation package. As
  such, a vector of data is to be plotted, the \code{names} attribute of
  the vector \emph{needs} to contain the gene identifiers.  Likewise if
  a matrix of data is to be plotted, the \code{\link{rownames}}
  attribute of the matrix \emph{needs} to contain the gene identifiers.

  To date, \code{plot} can be called for vector data, while
  \code{\link{matplot}} and \code{\link{image}} can be called for matrix
  data.  Most additional plotting arguments can be passed down via
  \dots.  However, the idiogram function plots the axis
  independently. Currently, only the \code{cex.axis},\code{col.axis},
  and \code{font.axis} parameters are intercepted from \dots and
  redirected to the specialized \code{\link{axis}} call. Other
  parameters that effect the axis should be set via \code{\link{par}}.

  The function \code{midiogram} is a simple wrapper around
  \code{idiogram} to plot \bold{all the chromosomes from a particular organism}
  using sensible default values. The \code{"m"} refers to
  plotting \bold{m}ultiple idiograms.
}

\author{ Kyle Furge <kyle.furge@vai.org> and Karl Dykema
  <karl.dykema@vai.org>}

\seealso{\code{\link[annotate]{buildChromLocation}},
  \code{\link{Hs.cytoband}}, \code{\link{idiograb}}}

\examples{
library(idiogram)

## 
## NOTE:This requires an annotation package to work.
##      In this example packages "hu6800" and "golubEsets" are used.
##      They can be downloaded from http://www.bioconductor.org
##      "hu6800" is under MetaData, "golubEsets" is under Experimental
##      Data.


if(require(hu6800) && require(golubEsets)) {
   library(golubEsets)
   data(Golub_Train)

   hu.chr <- buildChromLocation("hu6800")
   ex <- assayData(Golub_Train)$exprs[,1]

   ## make sure the names() attribute is set correctly
   gN <- names(ex)
   gN[1:10]   

   idiogram(ex,hu.chr,chr="1")

   colors <- rep("black",times=length(ex))
   colors[ex > 10000] <- "red"
   pts <- rep(1,times=length(ex))
   pts[ex > 10000] <- 2
   idiogram(ex,hu.chr,chr="1",col=colors,pch=pts,font.axis=2,cex.axis=1)
   abline(v=0,col="darkgreen")

   ## An example of the dlim option. It is most useful for making
   ## consistant multi-panel plots
   colors <- rep("black",times=length(ex))
   colors[ex > 10000] <- "red"
   colors[ex < 0] <- "blue"

   idiogram(ex,hu.chr,chr="1",col=colors,xlim=c(-3000,21000))
   idiogram(ex,hu.chr,chr="1",col=colors,dlim=c(-100,7500),xlim=c(-3000,21000))
   idiogram(ex,hu.chr,chr="1",col=colors,dlim=c(-100,7500),xlim=c(-3000,10000))

   ## Using the identify function
   ip <- idiogram(ex,hu.chr,chr="1",col=colors,pch=19)
   #identify(ip$x,ip$y,labels=ip$labels)
   
} else print("This example requires the hu6800 and golubEsets data packages.")

## The example data is BAC array CGH data from J. Fridlyand's
## aCGH package

data(idiogramExample)

idiogram(colo.eset[,1],ucsf.chr,chr="1")

idiogram(colo.eset,ucsf.chr,chr="1",method="image")

idiogram(colo.eset,ucsf.chr,chr="1",method="image",col=topo.colors(50),grid.lty=c(1,NA))

idiogram(colo.eset,ucsf.chr,chr="1",method="mat",type="l")

## for a consistant multi-panel plot it can be helpful to force the data
## range within each panel to a defined range using 'dlim'
## This is similar to calling the 'midiogram' function

op <- par(no.readonly=TRUE)
par(mai=par("mai")*c(0.1,0.5,0.5,0.5))
layout(rbind(c(1:8),c(0,9:14,0),c(15:22)))

for(i in c(1:22)) {
  idiogram(colo.eset,ucsf.chr,chr=i,method="i",dlim=c(-1,1),margin="i",relative=TRUE)
}

par(op)
}

\keyword{hplot}

