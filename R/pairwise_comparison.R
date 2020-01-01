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
#' @import ggplot2
pairwise_comparison <- function( x, out='delta', P_overlap_threshold=0.2 ){
  DEBUG=FALSE
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
    
    # Decaprecated
    # ifelse(x$rightWidth[1] >= x$leftWidth[2] & x$rightWidth[2] >= x$leftWidth[1], 'same', 'diff')
    
    ## Check if any of the entries are NA, if so return single
    if ( any(is.na(x$rightWidth)) | any(is.na(x$leftWidth)) | length(x$leftWidth)==1 | length(x$rightWidth)==1  ){
      return('N/A')
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
        
        if ( DEBUG ){
          ggplot2::ggplot( data.frame(x=c(x$leftWidth, x$rightWidth), region=as.factor(c(1,2,1,2))) ) +
            ggplot2::geom_vline( ggplot2::aes(xintercept=x, col=region) )
        }
      }
      
      if ( DEBUG ){
        ggplot2::ggplot( data.frame(x=c(x$leftWidth, x$rightWidth), region=as.factor(c(1,2,1,2))) ) +
          ggplot2::geom_vline( ggplot2::aes(xintercept=x, col=region) )
      }
    }
    
    return( ifelse(P_overlap > P_overlap_threshold, 'same', 'diff') )
    
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
      
      if ( DEBUG ){
        ggplot2::ggplot( data.frame(x=c(x$leftWidth, x$rightWidth), region=as.factor(c(1,2,1,2))) ) +
          ggplot2::geom_vline( ggplot2::aes(xintercept=x, col=region) )
      }
    }
    return( P_overlap )
  } else {
    cat(red(underline(out), ' is not a valid option. One of: {"delta","strcat","isequal, samePeak, percentOverlap"}'))
  }
}
