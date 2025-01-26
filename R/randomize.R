#' Randomize samples
#'
#' Randomize samples into two groups
#'
#' @param orderSamples a vector of sample IDs
#' @return A data frame of sample IDs and their associated designated group
#' @examples
#' randomize(c("1","2","3"))
#' randomize(data/orderSamples.rda)
randomize <- function(orderSamples){
  sample_samples <- sample(orderSamples,round(length(orderSamples)/2))
  test_samples <- prob::setdiff(orderSamples,sample_samples)
  randomized_samples <- dplyr::bind_rows(dplyr::tibble(sample=sample_samples,randomized="graph"),dplyr::tibble(sample=test_samples,randomized="test"))
  randomized_samples
}
