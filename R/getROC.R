#// **********************************************************************************************
#//                                       getROC.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title get Receiver Operator Characterostocs plot and data
#' @description This function can be used to generate the pseudo-recall plot
#' 
#' @param data A data.table containing p-value/q-value score and true class labels
#' @param scores A character vector of the data column name scores to compute the pseudo-recall for
#' @param p.title A character to label plot
#' @param fdr_threshold A numeric vector to add fdr threshold lines to the plot. Default: fdr_threshold = c(0.01, 0.05)
#' @return A list containing the ROC data and  a ggplot
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
#' @importFrom ROCR prediction performance
#' @importFrom RColorBrewer brewer.pal
getROC <- function( data, scores, p.title = NULL, fdr_threshold = c(0.01, 0.05) ){
  require(ROCR)
  require(RColorBrewer)
  if ( F ){
  data <- bm_struct$bm
  scores <- c( "osw_fdr", "ipf_fdr")
  scores <- c( "osw_fdr", "osw_aligned_fdr")
  p.title <- NULL
  fdr_threshold = c(0.01, 0.05)
  current_score <- scores[1]
  }
  ## Prepare ROC data
  roc_data <- lapply(scores, function(current_score){
    current_metric <- paste(gsub("^(\\w+)_.*", '\\1', current_score), '_TP', sep='')
    pred <- prediction(1-dplyr::select(data, current_score), dplyr::select(data, current_metric ) )
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    x<-perf@x.values[[1]]
    y<-perf@y.values[[1]]
    y[length(y)]<-y[length(y)-1]
    results <- data.table(fpr=x, tpr=y)
    results$method <- current_score
    cat( sprintf( "Evaluating %s at %s resulting AUC: %s\n", current_score, current_metric, performance(pred, measure = "auc")@y.values[[1]] ) )
    return(results)
  })
  roc_data <- rbindlist(roc_data)
  
  ## Plot ROC data
  g <- ggplot( ) +
    geom_line( data=roc_data, mapping=aes(y=tpr, x=fpr, col=method ), show.legend = T ) +
    geom_point( data=roc_data, mapping=aes(y=tpr, x=fpr, col=method ), size=1.3, show.legend = T )

  ## Add FDR Thresholds
  fdr_vline_values <- rep(2, length(fdr_threshold))
  names(fdr_vline_values) <- fdr_threshold
  fdr_colors <- suppressWarnings( brewer.pal(n = length(fdr_threshold), name = "Set1") )
  fdr_colors <- fdr_colors[1:length(fdr_threshold)]
  g <- g + 
    geom_vline( data=data.frame(fdr=fdr_threshold, linetype=as.character(fdr_threshold)), mapping=aes(xintercept = fdr, linetype=linetype), color=fdr_colors, show.legend = F)  +
    scale_linetype_manual( name = "FDR Threshold", labels = as.character(fdr_threshold), values = fdr_vline_values) 
  
  
  g <- g + ggtitle( p.title ) +
    labs( y = 'True positive rate', x = 'False positive rate') +
    theme( plot.title = element_text(hjust=0.5, size = 20),
           plot.subtitle = element_text(hjust=0.5),
           axis.title = element_text(size = 18),
           axis.text.y = element_text(size = 14),
           axis.text.x = element_text(angle = 90, size = 12),
           legend.text = element_text(size=14),
           legend.title = element_text(size = 16),
           legend.position = "right") +
    guides( col = guide_legend( override.aes = list(linetype = 0 ), nrow=5 ),
            linetype = guide_legend( override.aes = list(linetype = 0, color = fdr_colors) )
    )
  
  return( list(data=roc_data, plot=g))
}