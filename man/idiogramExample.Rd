\name{idiogramExample}
\alias{colo.eset}
\alias{idiogramExample}
\alias{ucsf.chr}
\alias{vai.chr}
\docType{data}
\title{data included for idiogram package examples}
\description{
   colo.eset and ucsf.chr are included for use in example plots. 
   vai.chr is a chromLocation object included for examples in the package 'reb'.
   Please see the aCGH Bioconductor package for more information.
}
\usage{data(idiogramExample)}

\source{
  http://www.bioconductor.org/repository/devel/package/html/aCGH.html
}




\examples{
	library(idiogram)
	data(Hs.cytoband)
	data(idiogramExample)
	idiogram(colo.eset@exprs[,1],ucsf.chr,chr="1")

}

\keyword{datasets}
