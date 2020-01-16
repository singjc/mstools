#// **********************************************************************************************
#//                         Mode.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title    Mode
#' @description Get the mode for a vector
#' 
#' @param x A vector 
#' @return The mode of the vector 
#' 
Mode <- function(x) {
  ux <- unique(x)
  ux[which(tabulate(match(x, ux)) == max(tabulate(match(x, ux))))]
}