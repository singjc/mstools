#// **********************************************************************************************
#//                                       getFDREstimatePerformance.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title get qunatification per run plot
#' @description This function can be used to generate a qunatification box plot per run
#' 
#' @param data A data.table containing intensity information
#' @param p.title A character to label plot
#' @param true_label A character vector indicating the column header name containing true class labels
#' @param global_fdr_label A character vector indicating the column header name containing the global fdr estimates
#' @param local_fdr_label A chracter vector indicating the column header name containing the local fdr estimates
#' @return A list containg a data.table of the data used to make performance and plot, and a ggplot object of the performance plot
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
#' @importFrom dplyr %>% select filter 
#' @importFrom data.table data.table
#' @importFrom ROCR prediction performance
#' @importFrom zoo rollapply
#' @importFrom ggplot2 ggplot geom_line geom_abline ggtitle theme element_text
#' @importFrom RColorBrewer brewer.pal
getFDREstimatePerformance <- function( data, p.title = NULL, true_label="TP", global_fdr_label='ipf_fdr', local_fdr_label='ipf_lfdr' ){
  
  # require(ROCR)
  # require(zoo)
  
  if ( F ){
    dim(ipf_ptm_fdr)
    data <- merge(ipf_df, benchmark, by=c("Run_Grouping", "FullPeptideName"), all.x=T, all.y=F)
    # data <- data[!is.na(data$ipf_lfdr),]
    data$GT[is.na(data$GT)] <- 0
    true_label <- "GT"
    global_fdr_label <- 'ipf_fdr'
    local_fdr_label <- 'ipf_lfdr'
    
    data <- ipf_df_fdr
    true_label <- "TP"
  }
  
  pred <- ROCR::prediction(1-dplyr::select(data, global_fdr_label), dplyr::select(data, true_label ) )
  perf <- ROCR::performance(pred, measure = "pcfall")
  datat <- data[order(dplyr::select(data, local_fdr_label)), ]
  
  global.data <- data.table()
  global.data$`Estimated FDR/fdr` <- 1-perf@x.values[[1]]
  global.data$`FDR/fdr` <- perf@y.values[[1]]
  global.data$measure <- "global FDR"
  
  local.data <- data.table()
  local.data$`Estimated FDR/fdr` <- zoo::rollapply(dplyr::select(datat,local_fdr_label), width = 500, mean, by = 1)
  local.data$`FDR/fdr` <- 1-zoo::rollapply(dplyr::select(datat,true_label), width = 500, mean, by = 1)
  local.data$measure <- "local FDR"
  
  p.data <- merge(global.data, local.data, by = c("Estimated FDR/fdr", "FDR/fdr", "measure"), all=TRUE)
  
  p.out <- ggplot() + 
    geom_line(data=p.data, aes(x=`Estimated FDR/fdr`, y=`FDR/fdr`, col=measure)) +
    geom_abline(intercept=0, slope=1, col="grey", alpha=0.9, linetype='dashed') +
    ggtitle( p.title ) +
    theme( plot.title = element_text(hjust=0.5, size = 20),
           plot.subtitle = element_text(hjust=0.5),
           axis.title = element_text(size = 18),
           axis.text.y = element_text(size = 14),
           axis.text.x = element_text(angle = 90, size = 12),
           legend.text = element_text(size=14),
           legend.title = element_text(size = 16),
           legend.position = "right"  ) 
    
    # print( p.out )
  
  
  # pred <- prediction(1-dplyr::select(data, global_fdr_label), dplyr::select(data, true_label ))
  # perf <- performance(pred, measure = "pcfall")
  # datat <- data[order(dplyr::select(data, local_fdr_label)), ]
  # 
  # plot(1-perf@x.values[[1]],perf@y.values[[1]], type="l",xlab="Estimated FDR/fdr",ylab="FDR/fdr",xlim=c(0,0.2),ylim=c(0,0.2), lwd=2,col=rainbow(2)[1])
  # lines(rollapply(dplyr::select(datat,local_fdr_label), width = 500, mean, by = 1) ,1-rollapply(dplyr::select(datat,true_label), width = 500, mean, by = 1), type="l", lwd=2,col=rainbow(2)[2])
  # 
  # legend("bottomright", c("global FDR", "local fdr"), pch = 0, col=rainbow(2))
  # abline(a=0, b=1, col = "gray40", lty=2)
  
  return(list( p.data=p.data, p.plot=p.out ))
  
}