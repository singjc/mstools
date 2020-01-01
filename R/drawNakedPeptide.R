drawNakedPeptide_ <- function(df_lib, 
                              mod, 
                              pep,
                              in_sqMass,  
                              plotPrecursor=T,
                              plotIntersectingDetecting=T, 
                              plotIdentifying=T,
                              plotUniqueDetecting=T,
                              plotIdentifying.Unique=NULL, 
                              plotIdentifying.Shared=NULL, 
                              plotIdentifying.Against=NULL,
                              intersecting_mz=NULL, 
                              uni_mod_list=NULL, 
                              max_RT, 
                              min_RT, 
                              max_Int, 
                              in_osw, 
                              smooth_chromatogram=NULL, 
                              doFacetZoom=NULL, 
                              top_trans_mod_list=NULL, 
                              show_all_pkgrprnk=F, 
                              show_n_transitions=6,
                              FacetFcnCall=NULL, 
                              verbosity=-1, 
                              show_legend = TRUE ){
  
  cat( green('   --- Peptidoform: ', mod), '\n', sep='' )
  
  # Load OSW Merged df
  osw_df <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id=pep, mod_residue_position='', peak_group_rank_filter=T, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=F )
  if ( dim(osw_df)[1]==0 ){ cat(red(pep, ' was not found as a peak_rank_group=1 in osw file!!!, skipping...\n'),sep=''); return(list()) }
  
  osw_df %>%
    dplyr::filter( FullPeptideName %in% mod) %>%
    select( Charge ) %>%
    as.matrix() %>%
    Mode() -> Isoform_Target_Charge
  
  # Display other peak group rank features
  if ( show_all_pkgrprnk==T ){
    osw_df_all <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id=pep, mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=F )
    
    osw_df_all %>%
      dplyr::filter( FullPeptideName %in% mod) %>%
      dplyr::filter( Charge %in% Isoform_Target_Charge )-> osw_df_all_filtered
    
    RT_Table <- table( osw_df_all_filtered$RT )
    
    RT_pkgrps <- as.numeric(names(RT_Table)[RT_Table==2])
    if ( length(RT_pkgrps)==0 ){
      cat(red('WARNING: There were no common RT pkgrps found, will plot all pkgrps...\n'))
      RT_pkgrps <- as.numeric(names(RT_Table))
    }
    rm(osw_df_all, osw_df_all_filtered, RT_Table)
    
  } else {
    RT_pkgrps <- NULL
  }
  

  
  # uni_mod_list <- uni_mod
  # mod=uni_mod
  # Isoform_Target_Charge=3
  plot_list <- list()
  # max_Int <- 0
  ###########################
  ##     PLOT PRECURSOR    ##
  ###########################
  if ( plotPrecursor==T ){
    
    g <- ggplot()
    g <- getXIC_( g, df_lib, mod, in_sqMass, transition_type='precursor', intersecting_mz=NULL, uni_mod_list, max_RT, min_RT, max_Int, in_osw=NULL, smooth_chromatogram=smooth_chromatogram, doFacetZoom=F, top_trans_mod_list=NULL, Isoform_Target_Charge=Isoform_Target_Charge, verbosity=verbosity )
    max_Int <- g$max_Int
    g <- g$graphic_obj
  } else {
    g <- ggplot()
  }
  
  #################################
  ##   DETECTING TRANSITIONS     ##
  #################################
  
  if ( plotIntersectingDetecting==T |  plotUniqueDetecting==T ){
    g <- getXIC_( g, df_lib, mod, in_sqMass, transition_type='detecting', intersecting_mz=NULL, uni_mod_list, max_RT, min_RT, max_Int, in_osw=NULL, smooth_chromatogram=smooth_chromatogram, doFacetZoom=F, top_trans_mod_list=NULL, show_n_transitions=show_n_transitions, verbosity=verbosity )
    max_Int <- g$max_Int
    g <- g$graphic_obj
  }
  
  ###################################
  ##    IDENTIFYING TRANSITIONS   ###
  ###################################
  if (plotIdentifying==T){
    g <- getXIC_( g, df_lib, mod, in_sqMass, transition_type='identifying', intersecting_mz=NULL, uni_mod_list, max_RT, min_RT, max_Int, in_osw=NULL, smooth_chromatogram=smooth_chromatogram, doFacetZoom=F, top_trans_mod_list=NULL, plotIdentifying.Unique=plotIdentifying.Unique, plotIdentifying.Shared=plotIdentifying.Shared, plotIdentifying.Against=plotIdentifying.Against, show_n_transitions=show_n_transitions, show_legend=show_legend, verbosity=verbosity )
    max_Int <- g$max_Int
    g <- g$graphic_obj
  } else {
    cat(red('-- Identifying Transitions were not found for: ', underline(mod)), '\n', sep='')
  }
  
  ###################################
  ##     ADD OSW RESULTS INFO     ###
  ###################################
  g <- getXIC_( g, df_lib, mod, in_sqMass, transition_type='none', intersecting_mz=NULL, uni_mod_list, max_RT, min_RT, max_Int, in_osw, doFacetZoom=doFacetZoom, top_trans_mod_list=NULL, Isoform_Target_Charge=Isoform_Target_Charge, RT_pkgrps=RT_pkgrps, FacetFcnCall=FacetFcnCall, verbosity=verbosity, show_legend = show_legend  )
  max_Int <- g$max_Int
  g <- g$graphic_obj
  
  plot_list[[mod]] <- g
  
  graphics.off()
  final_g <- (arrangeGrob(grobs=plot_list, nrow=length(mod)))
  final_g_list <- list()
  final_g_list[[1]] <- final_g
  # grid.draw(final_g)
  
  return( final_g_list )
}
