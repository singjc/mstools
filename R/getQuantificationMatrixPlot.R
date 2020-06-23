getQuantificationMatrixPlot <- function( data_list, filter_threshold, filter_scores, filter_GT_regex=NULL, use_seed=2020, sample_n_peptides=200, in_atleast_x_runs=1  ){
  
  if ( F ){
    data_list <- bm_struct$org_data
    filter_threshold <- 0.05
    filter_scores <- list(osw="osw_fdr", ipf="ipf_fdr", ipf_aligned="ipf_aligned_fdr" )
    filter_GT_regex <- "TP==1"
    use_seed=2020
    sample_n_peptides=100 
    in_atleast_x_runs=1
  }
  
  best_data_list <- lapply(names(data_list), function(data_i){
    
    current_data <- data_list[[data_i]] 
    filter_score <- filter_scores[[data_i]]
    
    current_data_best <- filterbestpeak(current_data, which_score = filter_score, fdr_threshold = filter_threshold)
    
    if ( !is.null(filter_GT_regex) ){
      current_data_best %>%
        dplyr::filter( eval(parse(text=filter_GT_regex)) ) -> current_data_best
    }
    ## Get column scores
    score_columns <- grep(".*lfdr|fdr", colnames(current_data_best), value = T)
    ## Make Sequence Column
    current_data_best$Sequence <- gsub("\\.?\\(([^()]+)\\)", "", current_data_best$FullPeptideName) 
    ## Actual RT and boundaries, not including GT
    RT_boundary_columns <- grep("^RT$|^leftWidth$|^rightWidth$|^RT.y$|^leftWidth.y$|^rightWidth.y$", colnames(current_data_best), value = T)
    ## Add Run_ID Column
    current_data_best$Run_ID <- gsub(".*_(Run_[:alnum:]*)", "\\1", current_data_best$Run_Grouping)
    
    current_data_best %>%
      dplyr::select( Run_Grouping, FullPeptideName, Sequence, RT_boundary_columns, Intensity, score_columns,  contains("method"), Run_ID ) -> current_data_best
    
    ## Rename scoring columns and Rt boundary columns for easier row binding
    setnames(current_data_best, old=score_columns, new=c("lfdr", "fdr"))
    setnames(current_data_best, old=RT_boundary_columns, new=c("RT", "leftWidth", "rightWidth"))
    
    ## Add current data separator
    current_data_best$method <- toupper(data_i)
    return(current_data_best)
    
  })
  data <- data.table::rbindlist(best_data_list)
  
  set.seed(use_seed)
  
  data %>%
    dplyr::group_by( FullPeptideName ) %>%
    dplyr::add_count() %>%
    dplyr::filter(n>in_atleast_x_runs) -> tmp
  
  
  rand_peptides <- sample( unique(tmp$FullPeptideName), size = ifelse(sample_n_peptides>length(unique(tmp$FullPeptideName)), length(unique(tmp$FullPeptideName)), sample_n_peptides) )
  # rand_peptides <- unique(data$FullPeptideName)
  
  data %>%
    dplyr::filter( FullPeptideName %in% rand_peptides ) %>%
    dplyr::group_by( method, Run_ID, FullPeptideName ) %>%
    dplyr::add_count() %>%
    dplyr::ungroup() -> tmp
  
  tmp %>%
    dplyr::filter( ifelse(n>1, ifelse(fdr==min(fdr), T, F),T) ) %>%
    dplyr::select( FullPeptideName, method, Run_ID, Intensity ) %>%
    tidyr::pivot_wider( values_from = "Intensity", names_from = "Run_ID") %>%
    melt() -> data_subset
  
  data_subset %>%
    tidyr::pivot_wider( values_from = "value", names_from = "variable") -> quant_matrix
  
  q1 <- merge(data.table(FullPeptideName=rand_peptides), subset(quant_matrix, method=="OSW"), by = "FullPeptideName", all.x=T)
  q1$method <- "OSW"
  rownames(q1) <- q1$FullPeptideName
  
  q2 <- merge(data.table(FullPeptideName=rand_peptides), subset(quant_matrix, method=="IPF"), by = "FullPeptideName", all.x=T)
  q2$method <- "IPF"
  rownames(q2) <- q2$FullPeptideName
  
  q3 <- merge(data.table(FullPeptideName=rand_peptides), subset(quant_matrix, method=="IPF_ALIGNED"), by = "FullPeptideName", all.x=T)
  q3$method <- "IPF_ALIGNED"
  rownames(q3) <- q3$FullPeptideName
  
  
  mat_data <- (log2(as.matrix(q3[,-c(1:2)])))
  mat_data[is.na(mat_data)] <- 0
  
  data_scaled <- scale( mat_data )
  ord <- hclust( dist(data_scaled, method = "euclidean"), method = "ward.D" )$order
  
  q1$FullPeptideName <- factor( q1$FullPeptideName, levels = rownames(q1)[ord],  labels = rownames(q1)[ord] )
  q2$FullPeptideName <- factor( q2$FullPeptideName, levels = rownames(q2)[ord],  labels = rownames(q2)[ord] )
  q3$FullPeptideName <- factor( q3$FullPeptideName, levels = rownames(q3)[ord],  labels = rownames(q3)[ord] )
  
  quant_matrix <- rbindlist(list(q1,q2,q3))
  quant_matrix %>%
    melt() -> data_subset
  
  data_subset$variable <- factor( data_subset$variable, levels = stringr::str_sort( unique(data_subset$variable), numeric = T ), labels = stringr::str_sort( unique(data_subset$variable), numeric = T ))
  data_subset$method <- factor( data_subset$method, levels = c("OSW", "IPF", "IPF_ALIGNED"), labels = c("OSW", "IPF", "IPF_ALIGNED") )
  
  p <- ggplot(data_subset, aes( x=variable, y=FullPeptideName, fill=log2(value) ) ) + geom_tile() + facet_grid( ~ method) +
    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", na.value = "grey92", midpoint = median(log2(data_subset$value), na.rm = T) ) + 
    labs(x="Run Condition/Replicate", y="Peptide precursors") +
    theme(axis.text.x = element_text(angle=90),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank() ) 
  p
  
}