#// **********************************************************************************************
#//                         pairwise_comparison.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Pair-Wise Comparison
#' @description This function can be used to perform a pairwise comparison
#' 
#' @param x a vector or list containing two comparisons 
#' @param out a character vector specifying method to perform
#' @return A value for comparison
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
pairwise_comparison <- function( x, out='delta' ){
  if ( length(x) >2 ){
    cat(red$bold$underline('Warning! There is more than 2 values, will only use first 2!\n'))
  }
  if (out=='delta'){
    return( abs( x[1] - x[2]) )
  } else if (out=='strcat') {
    return( paste(x[1],x[2],sep='_') )
  } else if (out=='isequal'){
    return( equals(x[1],x[2]))
  }  else if (out=='samePeak'){
    ifelse(x$rightWidth[1] >= x$leftWidth[2] & x$rightWidth[2] >= x$leftWidth[1], 'same', 'diff')
  } else {
    cat(red(underline(out), ' is not a valid option. One of: {"delta","strcat","isequal"}'))
  }
}