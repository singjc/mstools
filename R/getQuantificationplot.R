#// **********************************************************************************************
#//                                       getQuantificationplot.R
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
#' @param fdr_threshold A numeric vector to add fdr threshold lines to the plot. Default: fdr_threshold =  0.05
#' @param scores A character vector of the data column name scores to filter data on                        
#' @param intensities A character vector inficating the column to use for intentisty information
#' @param true_peptides A character vector indicating the column header name containing true class labels
#' @param keep.only.TP A logical value indicating to only show intensity information for true positive peptides
#' @param p.title A character to label plot
#' @return A  a ggplot object
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
#' @importFrom dplyr %>% select filter group_by group_by ungroup
#' @importFrom ROCR prediction performance
#' @importFrom ggplot2 ggplot stat_boxplot geom_boxplot stat_summary geom_label ggtitle theme element_text 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang sym
getQuantificationplot <- function( data, fdr_threshold=0.05, scores="ipf_fdr", intensities="ipf_Intensity", true_peptides="ipf_TP", normalzied_relative_intensity=FALSE, keep.only.TP=TRUE, p.title=NULL, dataset="" ) {
  
  if ( F ){
    fdr_threshold <- 0.05
    data <- bm
    scores <- "ipf_fdr"
    intensities <- "ipf_Intensity"
    true_peptides <- "ipf_TP"
    normalzied_relative_intensity=T
    data <- bm_struct$bm
    
    scores=c("ipf_fdr", "ipf_aligned_fdr")
    
    
  }
  
  data %>%
    dplyr::select( Run, FullPeptideName, !!rlang::sym(true_peptides), !!rlang::sym(scores), !!rlang::sym(intensities)  ) %>%
    dplyr::filter( !!rlang::sym(scores) < fdr_threshold ) -> data_subset
  
  if ( keep.only.TP ){
    data_subset %>%
      dplyr::filter( !!rlang::sym(true_peptides)==1 ) -> data_subset
  }
  
  data_subset$`log2(Intensity)` <- log2(data_subset[[intensities]])
  
  if ( normalzied_relative_intensity ){
    data_subset <- ddply( data_subset,
                          .(FullPeptideName),
                          function( X ){
                            if ( 'Run_013' %in% X$Run ){
                              # print(X)
                              X$relative_Intensity <- X[[intensities]] /X[which(X$Run=='Run_013'),intensities]
                              # print(X)
                              return( X )# End Return
                            }
                          }
    )
    data_subset$`log2(Intensity / Intensity [1:0])` <- log2(data_subset$relative_Intensity)
    data_subset$y <- data_subset$`log2(Intensity / Intensity [1:0])`
    y_label <- "log2(Intensity / Intensity [1:0])"
  } else {
    data_subset$y <- data_subset$`log2(Intensity)`
    y_label <- "log2(Intensity)"
  }
  
  
  
  if (  dataset=="Synth_Dilution"  ){
    data_subset$Run_ID <- plyr::mapvalues(x = data_subset$Run, from = c('Run_013', 'Run_012', 'Run_011', 'Run_010', 'Run_009', 'Run_008', 'Run_007', 'Run_006', 'Run_005', 'Run_004', 'Run_003', 'Run_002', 'Run_001'), to = c("1:0", "1:1", "1:3", "1:4", "1:7", "1:9", "1:15", "1:19", "1:31", "1:39", "1:63", "1:79", "1:127") )
    data_subset$Run_ID <- factor( data_subset$Run_ID, levels=c("1:0", "1:1", "1:3", "1:4", "1:7", "1:9", "1:15", "1:19", "1:31", "1:39", "1:63", "1:79", "1:127"))
  } else {
    data_subset$Run_ID <- gsub(".*_(Run_[:alnum:]*)", "\\1", data_subset$Run)
    data_subset$Run_ID <- factor( data_subset$Run_ID, levels=unique(data_subset$Run_ID))
  }
  data_subset %>%  dplyr::add_count( Run_ID ) -> data_subset
  
  p.out <- ggplot( data = data_subset, aes(x=Run_ID, y = y) ) + 
    stat_boxplot(geom ='errorbar', width = 0.3) +
    geom_boxplot( width = 0.3, outlier.alpha = 0.35 ) +
    stat_summary(fun.y=median, geom="smooth", aes(group=1), lwd=1, color="red") +
    geom_label( data = dplyr::select( data_subset , Run, Run_ID, n ) %>% unique(), aes(x=Run_ID, y = max(data_subset$y)+0.5, label = n ), size = 6 ) +
    ggtitle( p.title ) +
    labs( y=y_label, x="" ) +
    theme( plot.title = element_text(hjust=0.5, size = 20),
           plot.subtitle = element_text(hjust=0.5),
           axis.title = element_text(size = 18),
           axis.text.y = element_text(size = 14),
           axis.text.x = element_text(angle = 90, size = 12),
           legend.text = element_text(size=14),
           legend.title = element_text(size = 16),
           legend.position = "right") 
  
  return( p.out )
  
}