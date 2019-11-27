#// **********************************************************************************************
#//                         getOSWReportsVisualization.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title mstools osw reports visualization environment for various stats plotting functions
#' @description This function encapsulates various plotting functions that can be used for plotting
#' various reports from the OSW file.
#' 
#' @author Justin Sing \url{https://github.com/singjc}
osw_reports_visualization <- new.env()

## Dummy Documentation
#' @title Plot distribution of d-score
#' @usage osw_reports_visualization$d_score_hist( osw, bins=20, position= 'dodge2' )
#' @param osw A data.table/data.frame of the OSW results
#' @param bins A numeric vector for number of bins to bin data into
#' @param position A character vector indicating how groupings should be plot. (Default: 'dodge2') See ggplot2:geom_histogram
#' @inheritParams ggplot2::geom_histogram
#' @name d_score_hist
NULL

#' @rdname d_score_hist
osw_reports_visualization$d_score_hist <- function( osw, bins = 20, position = 'dodge2' ) { 
  
  ## TODO: Add Checks for columns not found in osw
  ggplot2::ggplot( osw, ggplot2::aes(x=d_score, group=factor(decoy), fill=factor(decoy) )) + 
    ggplot2::geom_histogram( ggplot2::aes(y = (..count..) ), position = position, bins = bins,  ) +
    ggplot2::scale_fill_manual(values=c("darkgreen", "red")) +
    ggplot2::ggtitle( sprintf("Distribution of d-score (bins=%s)", bins) ) +
    ggplot2::labs( x = 'd-score', y = 'Number of Peptides') +
    ggplot2::theme( plot.title = ggplot2::element_text(hjust=0.5, size = 20),
                    axis.title = ggplot2::element_text(size = 18),
                    axis.text = ggplot2::element_text(size=16)) 
}


