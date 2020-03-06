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
#' @param out a character vector specifying method to perform. Options: ('delta', 'strcat', 'isequal', 'samePeak', 'percentOverlap')
#' @param P_overlap_threshold a numeric value specifying the curret off to use for considering a peak overlap as bein the same or not
#' @return A value for comparison
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
pairwise_comparison <- function( x, out='delta', P_overlap_threshold=0.2 ){
  
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  DEBUG=FALSE
  if ( length(x) >2 ){
    MazamaCoreUtils::logger.warn(red$bold$underline('Warning! There is more than 2 values, will only use first 2!\n'))
  }
  if (out=='delta'){
    return( abs( x[1] - x[2]) )
  } else if (out=='strcat') {
    return( paste(x[1],x[2],sep='_') )
  } else if (out=='isequal'){
    return( equals(x[1],x[2]))
  }  else if (out=='samePeak'){
    
    MazamaCoreUtils::logger.trace( '[mstools::pairwise_comparison::samePeak] Calculating if two peaks are the same/overlapping' )
    
    
    MazamaCoreUtils::logger.trace( sprintf("[mstools::pairwise_comparison::samePeak] Input: %s", as.character(str(x))) )

    ## Check if any of the entries are NA, if so return single
    if ( any(is.na(x$rightWidth)) | any(is.na(x$leftWidth)) | length(x$leftWidth)==1 | length(x$rightWidth)==1  ){
      # cat("Returning NA\n")
      return(NA)
      }
     
    ## Check for cases where either boundary is completely withing the other boundary
    if ( x$leftWidth[1]  > x$leftWidth[2] & x$rightWidth[1]  < x$rightWidth[2] | x$leftWidth[1]  < x$leftWidth[2] & x$rightWidth[1]  > x$rightWidth[2] ){
      P_overlap = 1
    } else {
      
      ## Check for cases where either boundary is completely withing the other boundary
      if ( x$leftWidth[1]  > x$leftWidth[2] & x$rightWidth[1]  < x$rightWidth[2] | x$leftWidth[1]  < x$leftWidth[2] & x$rightWidth[1]  > x$rightWidth[2] ){
        P_overlap = 1
      } else {
        
        A1 <- abs((x$rightWidth[1] - x$leftWidth[1]) * 1000)
        A2 <- abs((x$rightWidth[2] - x$leftWidth[2]) * 1000)
        # A_overlap = abs( (min(c(x$rightWidth[1], x$rightWidth[2])) - max(c(x$leftWidth[1], x$leftWidth[2]))) * (1000-0) )
        
        A_overlap = max(0, (min(c(x$rightWidth[1],x$rightWidth[2])) - max(c(x$leftWidth[1], x$leftWidth[2])))) * max(0,  max(1000))
        
        P_overlap = abs( (A_overlap) / (A1+A2-A_overlap) )
        
      }
      
    }
    
    MazamaCoreUtils::logger.trace( sprintf('[mstools::pairwise_comparison::samePeak] Border overlap: %s with overlap threshold being: %s\n', P_overlap, P_overlap_threshold) )
    
    return( ifelse(P_overlap > P_overlap_threshold, 1, 0) )
    
  } else if (out=='percentOverlap'){
    if ( DEBUG ){
      x<- list()
      x$rightWidth[1]<-7008.0#728.4#90
      x$leftWidth[1]<-6775.8#5104.2#667.8#60
      x$rightWidth[2]<-6908.011#5254.063#792.6#95
      x$leftWidth[2]<-6765.386#742.8#65
    }
    # ( x$rightWidth[1] >= x$leftWidth[2] & x$rightWidth[2] >= x$leftWidth[1] )
    # 
    # ( x$rightWidth[2] >= x$leftWidth[1] & x$rightWidth[1] >= x$leftWidth[2] )
    # cat( sprintf( 'Region1-right(%s) >= Region2-left(%s) & Region2-right(%s) >= Region1-left(%s)\n' , x$rightWidth[1], x$leftWidth[2], x$rightWidth[2], x$leftWidth[1] )     )
    
    ## Check if any of the entries are NA, if so return single
    if ( any(is.na(x$rightWidth)) | any(is.na(x$leftWidth)) | length(x$leftWidth)==1 | length(x$rightWidth)==1 ){
      return(-1)
    }
    
    ## Check for cases where either boundary is completely withing the other boundary
    if ( x$leftWidth[1]  > x$leftWidth[2] & x$rightWidth[1]  < x$rightWidth[2] | x$leftWidth[1]  < x$leftWidth[2] & x$rightWidth[1]  > x$rightWidth[2] ){
      P_overlap = 1
    } else {
      
      A1 <- abs((x$rightWidth[1] - x$leftWidth[1]) * 1000)
      A2 <- abs((x$rightWidth[2] - x$leftWidth[2]) * 1000)
      # A_overlap = abs( (min(c(x$rightWidth[1], x$rightWidth[2])) - max(c(x$leftWidth[1], x$leftWidth[2]))) * (1000-0) )
      
      A_overlap = max(0, (min(c(x$rightWidth[1],x$rightWidth[2])) - max(c(x$leftWidth[1], x$leftWidth[2])))) * max(0,  max(1000))
      
      P_overlap = abs( (A_overlap) / (A1+A2-A_overlap) )
      
    }
    return( P_overlap )
  } else {
    MazamaCoreUtils::logger.error(red(underline(out), ' is not a valid option. One of: {"delta","strcat","isequal, samePeak, percentOverlap"}'))
  }
}
