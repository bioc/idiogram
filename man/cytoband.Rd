\name{cytoband-class}
\docType{class}
\alias{cytoband-class}

\title{Class "cytoband" }
\description{  Cytogenetic banding information }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cytoband", ...)}.
    
}
\section{Slots}{
  \describe{
    \item{\code{stain}:}{typical staining designation from classical cytogenetics}
    \item{\code{band}:}{character string representing the band name/number }
    \item{\code{start}:}{position in basepairs for the start of a given band }
    \item{\code{end}:}{position in basepairs for the end of a given band }
    \item{\code{length}:}{ length in basepairs for a given band }
  }
}

\author{ Karl Dykema <karl.dykema@vai.org> }
\seealso{
  See Also as \code{\link{Hs.cytoband}}
}
\examples{
## None
}
\keyword{classes}
