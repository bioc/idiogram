\name{ideograb}
\alias{ideograb}
\title{ ideograb }
\description{
  \code{ideograb} reads the position of the graphics pointer
  when the (first) mouse button is pressed. \code{ideogram} also reads
  a \emph{second} position of the graphics pointer after another
  mouse button press. The two points selected are used to define a
  diagonal line from which a bounding box will be constructed. It then
  retrieves the gene identifiers of the points that lie within the
  bounding box.
}
\usage{
ideograb(ideo, show.box = TRUE, brush = NULL, ...)
}
\arguments{
  \item{ideo}{ point coordinates and cooresponding gene identifiers
    from an \code{ideogram} function call}
  \item{show.box}{ boolean. if TRUE, a box is drawn showing the selected
    region}
  \item{brush}{ a color to highlight the points within the selected region}
  \item{\dots}{ additional plotting paramters passed to \code{points} to
    modify the points within the selected region}
}

\details{
  Coordinates can be passed in a plotting structure (a list with
  \code{x}, \code{y}, and \code{labels} components). Typically this is
  generated from a call to \code{ideogram}.
  
  The points selected are used to define the top-left and bottom-right
  locations \emph{or} the bottom-left and top-right locations for
  bounding box. These locations can be selected in any order.
  A character vector of the labels of all the points that lie within the
  selected region is returned.
}
\value{
  A character vector of gene identifiers
}

\author{ Karl Dykema <karl.dykema@vai.org> }
\seealso{ \code{\link{ideogram}} }

\examples{
data(ideogramExample)
ip <- ideogram(colo.eset@exprs[,1],ucsf.chr,chr="1")
if(interactive()) ideograb(ip,brush="red")
}
\keyword{ iplot }
