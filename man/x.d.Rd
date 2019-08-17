\name{x.d}
\alias{x.d}
\docType{data}
\title{
A simulated Dirichlet-multinomial based microbiome data.
}
\description{

It is simulated by Dirichlet-multinomial model with 100 individuals and 30 taxa, the sequencing depth is 1000 and dispersion is 0.2. We divide the 100 individuals into 4 groups with 10, 20, 30, 40 in each group, and induce abundance differences by increasing and decreasing the proportions of six taxa.
}
\usage{data(x.d)}
\format{
  The dataset \code{x.d} has 100 rows represent individuals, and 30 columns show OTUs.
}

\seealso{
\code{\link{sim_x}}
}

\examples{
data(x.d)
}
\keyword{datasets}
