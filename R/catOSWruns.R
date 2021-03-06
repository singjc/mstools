#// **********************************************************************************************
#//                         catOSWruns.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Concatenate different run OSW tables into one table
#' @description This function can be used to extract OSW data for several runs from an OSW file, and concatenate
#' them into one table with a run id
#' 
#' @param sqMass_files An array of character paths to the chromatograms that have corresponding names in the OSW file.
#' @param in_osw A character vector of the absolute path and filename of the osw file. (Must be .osw format)
#' @param which_m_score A character vector indicating which m_score to use. (Options: m_score or ms2-m_score. Default: 'm_score')
#' @param m_score_filter A numeric vector indicating the threshold m_score cut-off value for filtering. (Defualt: 0.05)
#' @param  report_top_single_result A logical value indificating to keep only the results with the lowest m_score. (Default: TRUE)
#' @param mod_exclusion A vector of unimod record id modifications to exclude, and replace by nothing. Example: c("35")
#' @param mod_grouping A dataframe contained uninod ids and random characters to group isomers by. Example: data.frame( record_id=c(4, 259, 267), rand_char_id="B", stringsAsFactors = F )
#' @return A data.table containing Run ID information
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom crayon blue bold underline 
#' @importFrom dplyr %>% select filter group_by ungroup add_count 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
catOSWruns_ <- function( sqMass_files, in_osw, which_m_score='m_score', m_score_filter=0.05, report_top_single_result=T, run_sub_expression=NULL, mod_exclusion=NULL, mod_grouping=NULL, ... ){
  
  DEBUG=F
  if ( DEBUG ) {
    in_sqMass="/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/phospho_enriched_U2OS/Results/Results_New_Lib/yanliu_I170114_016_PhosNoco4_SW/yanliu_I170114_016_PhosNoco4_SW_osw_chrom.sqMass"
    in_osw = "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/phospho_enriched_U2OS/Results/Results_New_Lib/pyprophet/merged_runs.osw"
    
    ## Synth_Phso
    in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/lower_product_mz_threshold/pyprophet/group_id/merged_runs_group_id_MS1MS2_intergration_ipf.osw"
    in_sqMass <-  "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/lower_product_mz_threshold//Dilution_1_0/chludwig_K150309_013_SW_0_osw_chrom.sqMass" 
    
    ## Christian sgolay 24
    in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/pyprophet_parametric/merged_runs_group_id_MS1MS2_intergration.osw"
    in_sqMass <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/sqMass/lgillet_L160915_002-Manchester_dirty_phospho_-_Pool_U1_-_SW.mzML.gz_osw_chrom.sqMass"
    in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_rt_extraction_window_750_extra_rt_extraction_window_300/pyprophet/LDA/merged_runs_group_id_MS1MS2_intergration_contexts_experiment_wide.osw"
    which_m_score="m_score"
    m_score_filter=0.01
    report_top_single_result=T
    decoy_filter = T; ipf_score=T
    
    ## Georges Synth_Phos Results, run by me
    in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/From_George/02032020/ipfm_nl_JS/psgs.osw"
    in_sqMass <- "/project/def-hroest/data/synth_phospho_pep/mzML/chludwig_K150309_001b_SW_1_64.mzXML.gz"
    m_score_filter=1
    report_top_single_result=T
    decoy_filter = T; ipf_score=F
    which_m_score="ms2_m_score"
  }
  
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  list_comparisons <- lapply(sqMass_files, function( in_sqMass ){
    
    run_name <- gsub('[_osw_]*chrom[.]sqMass|[_osw.]*chrom[.]sqMass|\\..*', '', basename(in_sqMass))
    
    MazamaCoreUtils::logger.info( paste('Reading in results for: ', crayon::blue$bold$underline(run_name), '\n', sep='' ))
    # Extract OpenSwath REsults for Specfific run
    osw_df <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id='', mod_residue_position='', peak_group_rank_filter=F, report_top_single_result=report_top_single_result, pep_list='', mod_exclusion=mod_exclusion, mod_grouping=mod_grouping, ... )
    
    # osw_df <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id='', mod_residue_position='', peak_group_rank_filter=F, report_top_single_result=report_top_single_result, pep_list='', decoy_filter = decoy_filter, ipf_score = ipf_score, mod_exclusion=mod_exclusion, mod_grouping=mod_grouping )
    
    # # Original OSW Peptide Names
    # osw_pep_names <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', osw_df$FullPeptideName)))))
    # 
    # #********************#
    # #***     DEBUG    ***#
    # #********************#
    # if ( DEBUG ){
    #   dim( osw_df )
    #   
    #   osw_df %>%
    #     dplyr::select( id_ipf_peptide, id_peptide, id_precursor, RT, leftWidth, rightWidth, Sequence, FullPeptideName, ipf_FullPeptideName, Charge, mz, ipf_pep, peak_group_rank, ms2_m_score, m_score ) %>%
    #     dplyr::filter( Sequence=="ADEICIAGSPLTPR")
    #   
    # }
    # 
    # if ( which_m_score=='m_score' ){
    #   # Keep only Rows that correspond to the correct Assay
    #   osw_df %>% dplyr::filter( osw_pep_names == osw_df$ipf_FullPeptideName ) -> osw_df
    # }
    
    #********************#
    #***     DEBUG    ***#
    #********************#
    if ( DEBUG ){
      dim( osw_df )
      
      osw_df %>%
        dplyr::select( id_ipf_peptide, id_peptide, id_precursor, RT, leftWidth, rightWidth, Sequence, FullPeptideName, ipf_FullPeptideName, Charge, mz, ipf_pep, peak_group_rank, ms2_m_score, m_score ) %>%
        dplyr::filter( Sequence=="ADEICIAGSPLTPR")
    }
    
    # Pre-Processing of Results.
    osw_df %>%
      dplyr::filter( !!as.name(which_m_score) < m_score_filter ) -> osw_df_fil2
    
    osw_df_fil2$Sequence_Group_id <- paste( osw_df_fil2$FullPeptideName, osw_df_fil2$Charge, osw_df_fil2$mz, sep='_')
    
    #********************#
    #***     DEBUG    ***#
    #********************#
    if ( DEBUG ){
      dim( osw_df_fil2 )
      
      osw_df_fil2 %>%
        dplyr::select( id_ipf_peptide, id_peptide, id_precursor, RT, leftWidth, rightWidth, Sequence, FullPeptideName, ipf_FullPeptideName, Charge, mz, ipf_pep, peak_group_rank, ms2_m_score, m_score ) %>%
        dplyr::filter( Sequence=="EDTEEISCR")
    }
    
    if ( report_top_single_result==T ){
      osw_df_fil2 %>% # Keep peptides that have an m_score less than 0.05
        # dplyr::select( Sequence_Group_id, id_ipf_peptide, id_peptide, id_precursor, RT, leftWidth, rightWidth, Sequence, FullPeptideName, ipf_FullPeptideName, Charge, mz, ipf_pep, peak_group_rank, ms2_m_score, m_score ) %>%
        # dplyr::filter( FullPeptideName=="ALTPERNT(Phospho)VPLKNNDSR(Label:13C(6)15N(4))") %>% # FOR DEBUGGING
        dplyr::group_by( Sequence_Group_id ) %>% # Group by FullPeptideName to further trim data
        dplyr::filter( !!as.name(which_m_score)==min(!!as.name(which_m_score)) ) %>% # Keep result of multiple form entries that passed intial m_score filtering. Keep results with the lowest m_score
        dplyr::ungroup() %>%
        dplyr::add_count( Sequence_Group_id ) -> osw_df_fil2  # Add column with counts for entries for each peptidoform. There might be cases where a peptidoform has two peak group ranks that have the same m_score (i.e. 0)
      osw_df_fil2 %>% 
        dplyr::group_by( Sequence_Group_id ) %>% # Group by FullPeptideName for a second pass of trimming the data
        dplyr::filter( ifelse( n>1, ifelse(peak_group_rank==min(peak_group_rank), T, F), T) ) %>%
        dplyr::ungroup() -> osw_df_fil2 # if the former ends up being true, only keep the results with peakgroup rank 1 annotation
      ## Remove n count column
      osw_df_fil2$n <- NULL
      ## Check to see if there are still more than one entry per run-peptide. 
      ## Depending on pyprophet scoring, peak_group_rank could be all run if
      ## using feature_id grouping.
      osw_df_fil2 %>%
        dplyr::group_by( Sequence_Group_id ) %>%
        dplyr::add_count() %>%
        dplyr::ungroup() -> osw_df_fil2
      
      osw_df_fil2 %>% 
        dplyr::group_by( Sequence_Group_id ) %>% # Group by FullPeptideName for a second pass of trimming the data
        dplyr::filter( ifelse( n>1, ifelse(d_score==max(d_score), T, F), T) ) %>%
        dplyr::ungroup() -> osw_df_fil2 # if the former ends up being true, only keep the results with highest d_score
      
      #********************#
      #***     DEBUG    ***#
      #********************#
      if ( DEBUG ){
        dim( osw_df_fil2 )
        
        osw_df_fil2 %>%
          dplyr::select( id_ipf_peptide, id_peptide, id_precursor, RT, leftWidth, rightWidth, Sequence, FullPeptideName, Charge, mz, ipf_pep, peak_group_rank, ms2_m_score, m_score ) %>%
          dplyr::filter( Sequence=="EDTEEISCR")
        
      }
      
    }
    if ( is.null(run_sub_expression) ){
      osw_df_fil2$RunID <- basename(osw_df_fil2$filename)
    } else {
      # 'lgillet_L\\d+_\\d+-Manchester_dirty_phospho_-_Pool_|_-_SW.mzML.gz'
      osw_df_fil2$RunID <- gsub(run_sub_expression, '', basename(osw_df_fil2$filename))
    }
    return( osw_df_fil2 )
  })
  
  
  Results_Comparison <- do.call( rbind, list_comparisons )
  
  return( Results_Comparison )
}




