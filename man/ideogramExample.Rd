\name{ideogramExample}
\alias{colo.eset}
\alias{ideogramExample}
\alias{ucsf.chr}
\alias{vai.chr}
\docType{data}
\title{data included for ideogram package examples}
\description{
   colo.eset and ucsf.chr are included for use in example plots. 
   vai.chr is a chromLocation object included for examples in the package 'reb'.
   
}
\usage{data(ideogramExample)}

\source{
  http://www.bioconductor.org/repository/devel/package/html/aCGH.html
}

\reference{aCGH Bioconductor package}


\examples{
	library(ideogram)
	data(Hs.cytoband)
	data(ideogramExample)
	ideogram(colo.eset@exprs[,1],ucsf.chr,chr="1")

}

\keyword{datasets}
