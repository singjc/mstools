#// **********************************************************************************************
#//                         pairwise_comparison_vectorized.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Pair-Wise Comparison vectorized
#' @description This function can be used to perform a pairwise comparison for vectors of boundaries
#' 
#' @param leftWidth.x left border of comparison 1
#' @param leftWidth.y left border of comparison 2
#' @param rightWidth.x right border of comparison 1
#' @param rightWidth.y right border of comparison 2
#' @param out a character vector specifying method to perform. Options: ('delta', 'strcat', 'isequal', 'samePeak', 'percentOverlap')
#' @param P_overlap_threshold a numeric value specifying the curret off to use for considering a peak overlap as bein the same or not
#' @return A logical value if boundaries are comparable
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
pairwise_comparison_vectorized <- function( leftWidth.x, leftWidth.y, rightWidth.x, rightWidth.y, out='samePeak', P_overlap_threshold=0.2 ){
  
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  if ( out=='samePeak' ){
    
    if ( F ){
      leftWidth.x <- 1970.4
      leftWidth.y <- 1973.24
      rightWidth.x <- 2043.6
      rightWidth.y <- 2005.95
    }
    
    MazamaCoreUtils::logger.trace( '[mstools::pairwise_comparison::samePeak] Calculating if two peaks are the same/overlapping' )
    
    boolean_overlap <-( (leftWidth.x  > leftWidth.y & rightWidth.x  < rightWidth.y) | (leftWidth.x  < leftWidth.y & rightWidth.x  > rightWidth.y) )
    
    A1 <- numeric(length(boolean_overlap))
    A2 <- numeric(length(boolean_overlap))
    
    A1[which(boolean_overlap)] <- abs((rightWidth.x[which(boolean_overlap)] - leftWidth.x[which(boolean_overlap)]) * 1000)
    A2[which(boolean_overlap)] <- abs((rightWidth.y[which(boolean_overlap)] - leftWidth.y[which(boolean_overlap)]) * 1000)
    
    ## Non-Perfect Overlaps
    A1[which(!boolean_overlap)] <- abs((rightWidth.x[which(!boolean_overlap)] - leftWidth.x[which(!boolean_overlap)]) * 1000)
    A2[which(!boolean_overlap)] <- abs((rightWidth.y[which(!boolean_overlap)] - leftWidth.y[which(!boolean_overlap)]) * 1000)
    
    ## NAs
    A1[which(is.na(boolean_overlap))] <- NA
    A2[which(is.na(boolean_overlap))] <- NA
    
    
    A_overlap = apply(cbind(0, (apply(cbind(rightWidth.x,rightWidth.y),1,min) - apply(cbind(leftWidth.x, leftWidth.y),1,max))), 1, max) * max(0,  max(1000))
    
    P_overlap = abs( (A_overlap) / (A1+A2-A_overlap) )
    
    
    MazamaCoreUtils::logger.trace( sprintf('[mstools::pairwise_comparison::samePeak] Done..\n') )
    
    return( ifelse(P_overlap > P_overlap_threshold, 1, 0) )
  }
}