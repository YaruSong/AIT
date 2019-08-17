\name{x.m}
\alias{x.m}
\docType{data}
\title{
A simulated multinomial based microbiome data.
}
\description{

It is simulated by multinomial model with 100 individuals and 30 taxa, the sequencing depth is 1000. We divide the 100 individuals into 4 groups with 10, 20, 30, 40 in each group, and induce abundance differences by increasing and decreasing the proportions of six taxa.
}
\usage{data(x.m)}
\format{
The dataset \code{x.m} has 100 rows represent individuals, and 30 columns show OTUs.
}
\seealso{\code{\link{sim_x}}}
\examples{
data(x.m)
}
\keyword{datasets}

