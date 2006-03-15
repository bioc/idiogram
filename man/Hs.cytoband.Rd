\name{Hs.cytoband}
\alias{Hs.cytoband}
\alias{Rn.cytoband}
\alias{Mm.cytoband}
\docType{data}
\title{Cytogenetic Banding information}
\description{
  Cytogenetic banding information for Homo sapiens (Hs) ,Mus musculus
  (Mm) and Rattus norvegicus (Rn)
}
\usage{data(Hs.cytoband)}
\format{
  Hs.cytoband is an environment containing 24 objects of class "cytoband" as defined by the idiogram package. Rn.cytoband and Mm.cytoband have also been included for the rat and mouse genomes.
}
\source{
  UCSC Genome Browser
  http://genome.ucsc.edu/
}
\examples{
data(Hs.cytoband)
cyto <- get("1",env=Hs.cytoband)
bands <- matrix(cyto@end-cyto@start,ncol=1)
barplot(bands,col="white")
}
\keyword{datasets}
