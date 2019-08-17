\name{x.dc}
\alias{x.dc}
\docType{data}
\title{
A simulated Dirichlet-multinomial based microbiome data existing covariates.
}
\description{

It is simulated by Dirichlet-multinomial model with 100 individuals and 6 taxa, the sequencing depth is 500. We divide the 100 individuals into 4 groups with 10, 20, 30, 40 in each group, and induce difference among taxa between four groups.
}
\usage{data(x.dc)}
\format{
  The dataset \code{x.dc} has 100 rows represent individuals, and 6 columns show OTUs.
}

\seealso{
\code{\link{sim_x_c}}
}

\examples{
data(x.dc)
}
\keyword{datasets}
