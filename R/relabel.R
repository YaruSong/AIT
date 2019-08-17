#' Reassigning values of a vector
#'
#' Reassigning values of a continuous or categorical vector. It is used to generate legal
#' value before applying independence test.
#'
#' @param y a continuous or categorical vector.
#' @return an integer vector with values range from 1 to k (k > 1).
#' @examples
#' n <- 10
#' y <- c(rep("G3", n),rep("G2", n), rep("G1", n/2))
#' y.r <- relabel(y)
#' @export
relabel <- function(y)
{
  if(length(which(is.na(y))) > 0){
    stop("There exists NA!")
  }
  uy <- unique(y[!is.na(y)])
  v <- vector(length = length(y[!is.na(y)]), mode = "integer")
  for(i in 1:length(uy)){
    v[which(y==uy[i])]=i
  }
  return(v)
}




