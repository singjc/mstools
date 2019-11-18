getdeltaRT_ <- function( osw_df ){
  "This function is used to get the difference in retention time between isomers."
  delta_RT_list <- list()
  unique_pep <- unique(osw_df$Sequence)
  for (pep in unique_pep){
    osw_df %>%
      filter( Sequence==pep ) -> pep_subset
    delta_RT_list[['deltaRT']] <- c(delta_RT_list[['deltaRT']], abs(combn(pep_subset$RT, m=2, diff)))
    delta_RT_list[['peptide']] <- c(delta_RT_list[['peptide']], unique(pep_subset$Sequence))
  }
  df <- data.table(matrix(unlist(delta_RT_list), nrow=length(delta_RT_list$deltaRT), byrow=F),stringsAsFactors=FALSE)
  colnames(df) <- c('deltaRT', 'Sequence')
  return( df )
}