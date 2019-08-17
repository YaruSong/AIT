\name{individual_index}
\alias{individual_index}
\docType{data}
\title{
Index information for samples.
}
\description{

Yatsunenko et al. (2012) conducted a study to investigate whether there exist taxa abundance differences across age and geography, where fecal specimens were collected from 528 healthy infants, children, and adults living in Amazonas (Venezuela), rural Malawi, and USA metropolitan areas. This data set \code{individual_index} contains "Age", "Country", "BMI", "Gender" and "Breast.fed" information for 528 individuals.

}
\usage{data(individual_index)}
\format{
The data \code{individual_index} has 528 rows and 5 columns:
\describe{
    \item{\code{Age}}{age ranges from 0.03 to 83, with 39 NA.}
    \item{\code{Country}}{including three region: "Malawi", "USA", "Venezuela".}
    \item{\code{BMI}}{Body Mass Index, ranges from 11 to 66, with 94 NA.}
    \item{\code{Gender}}{195 "Male", 298 "Female" and 35 NA.}
    \item{\code{Breast.fed}}{whether the child is breast-feeding (100) or formula-feeding (27), and 401 NA.}
    }
}

\examples{
data(individual_index)
}
\keyword{datasets}
\source{
Yatsunenko, T., Rey, F. E., Manary, M. J., Trehan, I., Dominguez-Bello, M. G., Contreras, M., Magris, M., Hidalgo, G., Baldassano, R. N., Anokhin, A. P., et al. (2012). Human gut microbiome viewed across age and geography. Nature, 486(7402), 222–227.

}
