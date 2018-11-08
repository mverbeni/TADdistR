#' Put results in dataframe according to ascending p value
#' @param  gr GRanges object
#' @return data frame
#' @keywords results, annotation
#' @export
#'
resume <- function(gr){
  data.df <- as(gr,"data.frame")
  data.df <- data.df[order(data.df$p.value),]
  return(data.df)
}
