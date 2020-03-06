#// **********************************************************************************************
#//                                       getIdentificationsBarplot.R
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
#' @param osw_df_bm A data.table containing osw results with bitmap information (TPs, FPs, TNs, FNs)
#' @param ipf_df_bm A data.table containing ipf results with bitmap information
#' @param fdr_threshold A numeric vector to add fdr threshold lines to the plot. Default: fdr_threshold =  0.05
#' @param p.title A character to label plot
#' @param dataset A character vector specifying special label mappings for specific dataset. Options: "", "Synth_Dilution"
#' @return A ggplot object
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom crayon blue red underline magenta bold 
#' @importFrom dplyr %>% select filter group_by summarise group_by ungroup
#' @importFrom plyr mapvalues
#' @importFrom ROCR prediction performance
#' @importFrom ggplot2 ggplot geom_bar ggtitle theme element_text scale_y_continuous scale_fill_manual
#' @importFrom RColorBrewer brewer.pal
getIdentificationsBarplot <- function( data, data_to_show=NULL, fdr_threshold = 0.05, p.title = NULL, dataset=""  ) {
  if (F){
    osw_df_bm <- bm_struct$osw_df_bm
    ipf_df_bm <- bm_struct$ipf_df_bm
    fdr_threshold = 0.05
    p.title = NULL
    dataset=NULL
    data <- bm_struct$bm_list
  }
  
  if ( is.null( data_to_show ) ) data_to_show <- names(data)
  
  results.out <- lapply(data_to_show, function( data_id ){
    current_data <- data[[data_id]]
    ### Filter at fdr threshold
    fdr_col_name <- grep(".*_fdr", colnames(current_data), value = TRUE)
    current_data  %>%
      dplyr::filter( !!rlang::sym(fdr_col_name) < fdr_threshold ) %>% as.data.table() -> current_data
    ### Remove method labels from class label prefix
    colnames(current_data)[which(grepl(".*_TP", colnames(current_data)))] <- "TP"
    colnames(current_data)[which(grepl(".*_FP", colnames(current_data)))] <- "FP"
    colnames(current_data)[which(grepl(".*_TN", colnames(current_data)))] <- "TN"
    colnames(current_data)[which(grepl(".*_FN", colnames(current_data)))] <- "FN"
    ### Set Run ID and Method
    current_data$method <- gsub("(\\w*)_fdr", "\\1", fdr_col_name)
    current_data$Run_ID <- gsub(".*_(Run_.*)$", "\\1", current_data$Run_Grouping)
    dplyr::select( current_data, Run_ID, method, TP, FP, TN, FN ) -> confusion_dt
    return( confusion_dt )
  })
  confusion_dt <- data.table::rbindlist( results.out )
  
  # ## OSW
  # ### Filter at fdr threshold
  # fdr_col_name1 <- grep(".*_fdr", colnames(osw_df_bm), value = TRUE)
  # osw_df_bm %>%
  #   dplyr::filter( !!rlang::sym(fdr_col_name1) < fdr_threshold ) %>% as.data.table() -> osw_df_bm
  # 
  # osw_df_bm$method <- gsub("(\\w*)_fdr", "\\1", fdr_col_name1)
  # osw_df_bm$Run_ID <- gsub(".*_(Run_.*)$", "\\1", osw_df_bm$Run_Grouping)
  # 
  # ## IPF
  # ### Filter at fdr threshold
  # fdr_col_name2 <- grep(".*_fdr", colnames(ipf_df_bm), value = TRUE)
  # ipf_df_bm %>%
  #   dplyr::filter( !!rlang::sym(fdr_col_name2) < fdr_threshold ) %>% as.data.table() -> ipf_df_bm
  # 
  # ipf_df_bm$method <- gsub("(\\w*)_fdr", "\\1", fdr_col_name2)
  # ipf_df_bm$Run_ID <- gsub(".*_(Run_.*)$", "\\1", ipf_df_bm$Run_Grouping)
  # 
  # confusion_dt <- merge( dplyr::select( osw_df_bm, Run_ID, method, TP, FP, TN, FN ), 
  #                        dplyr::select( ipf_df_bm, Run_ID, method, TP, FP, TN, FN ), 
  #                        by=c("Run_ID", "method", "TP", "FP", "TN", "FN"), all=T )
  
  confusion_dt %>%
    dplyr::group_by( method, Run_ID ) %>%
    dplyr::summarise( TP = sum(TP, na.rm = T),
                      FP = sum(FP, na.rm = T),
                      TN = sum(TN, na.rm = T),
                      FN = sum(FN, na.rm = T),
    ) %>%
    dplyr::ungroup() -> confusion_dt
  
  confusion_dt <- melt(confusion_dt)
  confusion_dt %>% dplyr::filter( grepl("*TP|*FP", variable)) -> confusion_dt
  confusion_dt$variable  <- paste(confusion_dt$method, confusion_dt$variable , sep='_')
  confusion_dt$variable <- as.character(confusion_dt$variable)
  confusion_dt$variable  <- factor( confusion_dt$variable, levels = rev(sort(unique(confusion_dt$variable))) )
  if (  dataset=="Synth_Dilution"  ){
    confusion_dt$Run_ID <- plyr::mapvalues(x = confusion_dt$Run, from = c('Run_013', 'Run_012', 'Run_011', 'Run_010', 'Run_009', 'Run_008', 'Run_007', 'Run_006', 'Run_005', 'Run_004', 'Run_003', 'Run_002', 'Run_001'), to = c("1:0", "1:1", "1:3", "1:4", "1:7", "1:9", "1:15", "1:19", "1:31", "1:39", "1:63", "1:79", "1:127") )
    confusion_dt$Run_ID <- factor( confusion_dt$Run_ID, levels=c("1:0", "1:1", "1:3", "1:4", "1:7", "1:9", "1:15", "1:19", "1:31", "1:39", "1:63", "1:79", "1:127"))
  } 
  
  TP_col_pallete <- colorRampPalette(c("green", "darkgreen"))
  
  FP_col_pallete <- colorRampPalette(c("red", "darkred"))
  
  fill_cols <- c(TP_col_pallete(length(unique(confusion_dt$method))), FP_col_pallete(length(unique(confusion_dt$method))))
  
  # fill_cols <- c("green1", "green4", "tomato", "red")
  
  names(fill_cols) <- unique(confusion_dt$variable)
  
  p.out <- ggplot( ) + 
    geom_bar( data=dplyr::select(confusion_dt, c(Run_ID, variable, value)) %>% dplyr::filter(grepl("*TP", variable)), 
              aes(y=value, x=Run_ID, fill=variable), stat="identity", position="dodge2", col="black"
    ) +
    geom_bar( data=dplyr::select(confusion_dt, c(Run_ID, variable, value)) %>% dplyr::filter(grepl("*FP", variable)), 
              aes(y=value*-1, x=Run_ID, fill=variable), stat="identity", position="dodge2", col="black"
    ) +
    ggtitle( p.title ) +
    scale_y_continuous(breaks = seq(-(ceiling(max(confusion_dt$value)/50)*50), (ceiling(max(confusion_dt$value)/50)*50), 50), labels = abs(seq(-(ceiling(max(confusion_dt$value)/50)*50), (ceiling(max(confusion_dt$value)/50)*50), 50))) +
    scale_fill_manual("legend", values = fill_cols ) +
    labs( x="", y= paste("Detected Peptides (", fdr_threshold, "FDR )", sep=" " ) ) +
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