nIdentificationsComparison <- function( sqMass_files, man_annotated_df, type_pool_id, man_annotated_col, in_osw){
  
  list_comparisons <- lapply(sqMass_files, function( in_sqMass ){
    
    run_name <- gsub('_osw_chrom[.]sqMass$', '', basename(in_sqMass))
    pool_id <- gsub('_-_SW.mzML.gz', '', gsub('lgillet_L\\d+_\\d+-Manchester_dirty_phospho_-_Pool_', '', run_name))
    
    
    man_annotated_df %>%
      dplyr::filter( get(type_pool_id, envir=as.environment(man_annotated_df))==pool_id ) -> df_sub
    
    n_spiked_in_peps <- length(unique(df_sub$`Sequence (stripped)`))
    n_spiked_in_mod_peps <- length( gsub('\\[CAM\\]','\\(UniMod:4\\)', gsub('\\[Oxi\\]','\\(UniMod:35\\)', gsub('\\[Pho\\]','\\(UniMod:21\\)', (df_sub$Peptide)))) )
    
    # Manual Annotation
    df_sub %>%
      dplyr::filter( get(man_annotated_col, envir=as.environment(df_sub))=='yes' ) -> df_sub_truth
    
    n_spiked_in_truth_peps <- length(unique(df_sub_truth$`Sequence (stripped)`))
    n_spiked_in_truth_mod_peps <- length( gsub('\\[CAM\\]','\\(UniMod:4\\)', gsub('\\[Oxi\\]','\\(UniMod:35\\)', gsub('\\[Pho\\]','\\(UniMod:21\\)', (df_sub_truth$Peptide)))) )
    
    
    # Extract OpenSwath REsults for Specfific run
    osw_df <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id='', mod_residue_position='', peak_group_rank_filter=T, pep_list='' )
    # Original OSW Peptide Names
    osw_pep_names <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', osw_df$FullPeptideName)))))
    # Keep only Rows that correspond to the correct Assay
    osw_df %>% dplyr::filter( osw_pep_names == osw_df$ipf_FullPeptideName ) -> osw_df
    
    # Pre-Processing of Results.
    osw_df %>%
      dplyr::filter( m_score < 0.05 ) %>% # Keep peptides that have an m_score less than 0.05
      dplyr::group_by(FullPeptideName) %>% # Group by FullPeptideName to further trim data
      dplyr::filter( m_score==min(m_score) ) %>% # Keep result of multiple form entries that passed intial m_score filtering. Keep results with the lowest m_score
      ungroup() %>%
      dplyr::add_count( FullPeptideName )  %>% # Add column with counts for entries for each peptidoform. There might be cases where a peptidoform has to peak group ranks that have the same m_score (i.e. 0)
      dplyr::filter( ifelse( n==2, ifelse(peak_group_rank==1, T, F), T) ) -> osw_df_fil2 # if the former ends up being true, only keep the results with peakgroup rank 1 annotation
    
    
    n_osw_peps <- length(unique(osw_df_fil2$Sequence)) # 47 Unique for Pool_M1
    n_osw_mod_peps <- length(osw_df_fil2$FullPeptideName) # 84 Unique for Pool_M1
    
    num_comparison_df <- data.frame(  run = pool_id,
                                      n_spiked_in = n_spiked_in_peps,
                                      n_spiked_in_mod = n_spiked_in_mod_peps,
                                      n_truth = n_spiked_in_truth_peps,
                                      n_truth_mod = n_spiked_in_truth_mod_peps,
                                      n_osw = n_osw_peps,
                                      n_osw_mod = n_osw_mod_peps
    )
  })
  
  
  Results_Comparison <- do.call( rbind, list_comparisons )
  
  return( Results_Comparison )
}




