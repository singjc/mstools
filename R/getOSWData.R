#// **********************************************************************************************
#//                         getOSWData.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing


#' @export
#' @title Extract data from OpenSwathWorkflow .osw results file
#' @description This function can be used to extract information from the OSW results file obtained from
#' running OpenSwathWorkflow
#' 
#' @param oswfile A character vector of the absolute path and filename of the osw results file. (Must be .osw format)
#' @param run_name A character vector for extraction of a specific run, this should be the same as the file name in the .OSW RUN table. (i.e. run_name='yanliu_I170114_016_PhosNoco4_SW.mzXML.gz', Default: '')
#' @param precursor_id A numeric value for a specific precursor id to extract information for. (Default: '')
#' @param peptide_id A string vector indicating a specific peptide to extract information for. (Default: '')
#' @param mod_peptide_id An array of two string vectors indicating a specific modified peptide sequence with both UniMod annotation and actual modification name to extract information for. I.e. c(ANS(Phos)SNSLK, ANS(UniMod:21)SNSLK) (Default: '') 
#' @param mod_residue_position A numeric value indicating the position of a modification. (Default: '')
#' @param peak_group_rank_filter A Logical value for filtering OSW results by Peak-Group Rank of 1. (Default: FALSE)
#' @param pep_list An arrary of string modified peptide sequences to extract information for. (Default: '')
#' @param mscore_filter A numeric value to filter results by MS2 q-value. (Default: '')
#' @param ipf_filter A numeric value to filter results by IPF PEP. (Default: '')
#' @param ipf_score A logical value to extract data using IPF Score results. (Default: FALSE)
#' @param ms2_score A logical value to extract data using MS2 Score results. (Default: TRUE)
#' @param  decoy_filter A logical value to filter decoys out of final results. (Default: TRUE)
#' @return A data.table containing OpenSwath Results information. 
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @import DBI
#' @import RSQLite
#' @import dbplyr
#' @import dplyr
getOSWData_ <- function ( oswfile,
                          run_name='',
                          precursor_id='',
                          peptide_id='',
                          mod_peptide_id='',
                          mod_residue_position='',
                          peak_group_rank_filter=FALSE, 
                          pep_list='',
                          mscore_filter='',
                          ipf_filter='',
                          ipf_score=FALSE,
                          ms2_score=TRUE,
                          decoy_filter=TRUE){
  
  DEBUG=F
  if ( DEBUG ) {
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/lower_product_mz_threshold/pyprophet/group_id/merged_runs_group_id_MS1MS2_intergration_ipf.osw"
    # run_name <- '/project/def-hroest/data/synth_phospho_pep/mzML/chludwig_K150309_013_SW_0.mzXML'
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/pyprophet_parametric/merged_runs_group_id_MS1MS2_intergration.osw"
    # run_name <- "/home/singjust/projects/def-hroest/data/christian_doerig_dataset/U_pools_mzML/sw/lgillet_L160915_024-Manchester_dirty_phospho_-_Pool_U12_-_SW.mzML.gz"
    #  run_name = "lgillet_L160915_028-Manchester_dirty_phospho_-_Pool_U14_-_SW.mzML.gz"
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/M_pools_re-run/pyprophet_parametric_no_precursor/merged_runs_group_id_MS1MS2_intergration.osw"
    # run_name <- "lgillet_L160918_002-Manchester_dirty_phospho_-_Pool_M1_-_SW.mzML.gz"
    # Load Requried Libraries
    # library(dplyr)
    # library(dbplyr)
  }
  
  # Connect to database
  osw_db <- DBI::dbConnect( RSQLite::SQLite(), oswfile )
  
  # Get RUN ID from database
  if ( run_name != '' ){
    run_id_df = getRunID_(oswfile, run_name)
    run_id_query = sprintf("INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID AND FEATURE.RUN_ID=(%s)", run_id_df$ID)
  } else {
    run_id_query = 'INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID'
  }
  if (precursor_id != ''){
    precursor_query = sprintf("INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID AND PRECURSOR.ID=(%s)", precursor_id)
  } else {
    precursor_query = 'INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID'
  }
  if( peptide_id != ''){
    peptide_query = sprintf("INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID AND PEPTIDE.UNMODIFIED_SEQUENCE=('%s')", peptide_id)
  } else {
    peptide_query = 'INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID'
  }
  if (mod_peptide_id!=''){
    if( ipf_score == TRUE){
      mod_peptide_query = sprintf("WHERE PEPTIDE_IPF.MODIFIED_SEQUENCE=('%s') OR PEPTIDE_IPF.MODIFIED_SEQUENCE=('%s')", mod_peptide_id[1], mod_peptide_id[2]) # PEPTIDE.MODIFIED_SEQUENCE use to correspond with PEPTIDE_IPF.MODIFIED_SEQUENCE. @ Justin
    } else {
      mod_peptide_query = sprintf("WHERE PEPTIDE.MODIFIED_SEQUENCE=('%s') OR PEPTIDE.MODIFIED_SEQUENCE=('%s')", mod_peptide_id[1], mod_peptide_id[2])
    }
  } else {
    mod_peptide_query = ''
  }
  if (mod_residue_position!=''){
    if ((mscore_filter != '') && (ipf_filter != '')){
      mod_residue_position_query = sprintf("AND INSTR(FullPeptideName,'(Phospho)')=(%s) OR INSTR(FullPeptideName,'(UniMod:21)')=(%s)", mod_residue_position, mod_residue_position)
    } else {
      mod_residue_position_query = sprintf("WHERE INSTR(FullPeptideName,'(Phospho)')=(%s) OR INSTR(FullPeptideName,'(UniMod:21)')=(%s)", mod_residue_position, mod_residue_position)
    }
  } else {
    mod_residue_position_query = ''
  }
  if (ipf_score == TRUE && ms2_score == TRUE){
    select_as_stmt = sprintf("SELECT FEATURE.ID AS feature_id,
       PEPTIDE_IPF.ID AS id_ipf_peptide,
       PEPTIDE.ID AS id_peptide,
       PRECURSOR.ID AS id_precursor,
       --PEPTIDE_IPF.MODIFIED_SEQUENCE || '_' || PRECURSOR.ID AS transition_group_id,
       PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.ID AS transition_group_id,
       PRECURSOR.DECOY AS decoy,
       RUN.ID AS run_id,
       RUN.FILENAME AS filename,
       FEATURE.EXP_RT AS RT,
       FEATURE.EXP_RT - FEATURE.DELTA_RT AS assay_rt,
       FEATURE.DELTA_RT AS delta_rt,
       FEATURE.NORM_RT AS iRT,
       PRECURSOR.LIBRARY_RT AS assay_iRT,
       FEATURE.NORM_RT - PRECURSOR.LIBRARY_RT AS delta_iRT,
       FEATURE.ID AS id,
       --PEPTIDE_IPF.UNMODIFIED_SEQUENCE AS Sequence,
       --PEPTIDE_IPF.MODIFIED_SEQUENCE AS FullPeptideName,
       PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
       PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
       PEPTIDE_IPF.MODIFIED_SEQUENCE AS ipf_FullPeptideName,
       PRECURSOR.CHARGE AS Charge,
       PRECURSOR.PRECURSOR_MZ AS mz,
       FEATURE_MS2.AREA_INTENSITY AS Intensity,
       FEATURE_MS1.AREA_INTENSITY AS aggr_prec_Peak_Area,
       FEATURE_MS1.APEX_INTENSITY AS aggr_prec_Peak_Apex,
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       --SCORE_MS1.PEP AS ms1_pep,
       SCORE_MS2.PEP AS ms2_pep,
       SCORE_IPF.PRECURSOR_PEAKGROUP_PEP AS precursor_pep,
       SCORE_IPF.PEP AS ipf_pep,
       SCORE_MS2.RANK AS peak_group_rank,
       SCORE_MS2.SCORE AS d_score,
       SCORE_MS2.QVALUE AS ms2_m_score,
       SCORE_IPF.QVALUE AS m_score")
    include_ipf_score = 'LEFT JOIN SCORE_IPF ON SCORE_MS2.FEATURE_ID = SCORE_IPF.FEATURE_ID'
    if (peak_group_rank_filter == TRUE){
      pk_grp_rnk_fil_query = 'INNER JOIN PEPTIDE AS PEPTIDE_IPF ON SCORE_IPF.PEPTIDE_ID = PEPTIDE_IPF.ID'
    } else {
      pk_grp_rnk_fil_query = 'LEFT JOIN PEPTIDE AS PEPTIDE_IPF ON SCORE_IPF.PEPTIDE_ID = PEPTIDE_IPF.ID'
    }
  } else if ((ipf_score == FALSE && ms2_score == TRUE)){
    select_as_stmt = sprintf("SELECT FEATURE.ID AS feature_id,
       PEPTIDE.ID AS id_peptide,
       PRECURSOR.ID AS id_precursor,
       PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.ID AS transition_group_id,
       PRECURSOR.DECOY AS decoy,
       RUN.ID AS run_id,
       RUN.FILENAME AS filename,
       FEATURE.EXP_RT AS RT,
       FEATURE.EXP_RT - FEATURE.DELTA_RT AS assay_rt,
       FEATURE.DELTA_RT AS delta_rt,
       FEATURE.NORM_RT AS iRT,
       PRECURSOR.LIBRARY_RT AS assay_iRT,
       FEATURE.NORM_RT - PRECURSOR.LIBRARY_RT AS delta_iRT,
       FEATURE.ID AS id,
       PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
       PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
       PEPTIDE.MODIFIED_SEQUENCE AS ipf_FullPeptideName,
       PRECURSOR.CHARGE AS Charge,
       PRECURSOR.PRECURSOR_MZ AS mz,
       FEATURE_MS2.AREA_INTENSITY AS Intensity,
       FEATURE_MS1.AREA_INTENSITY AS aggr_prec_Peak_Area,
       FEATURE_MS1.APEX_INTENSITY AS aggr_prec_Peak_Apex,
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       --SCORE_MS1.PEP AS ms1_pep,
       SCORE_MS2.PEP AS ms2_pep,
       SCORE_MS2.RANK AS peak_group_rank,
       SCORE_MS2.SCORE AS d_score,
       SCORE_MS2.QVALUE AS ms2_m_score")
    # print(select_as_stmt)
    if (peak_group_rank_filter == TRUE){
      pk_grp_rnk_fil_query = 'WHERE SCORE_MS2.RANK=1'
    } else {
      pk_grp_rnk_fil_query = ''
    }
    include_ipf_score = ''
    
  } else {
    cat( 'There was an error with the ms2_score or ipf_score argument!!!\n' )
    cat( sprintf("ms2_score: %s\n", ms2_score) )
    cat( sprintf("ipf_score: %s\n", ipf_score) )
  }
  
  if (pep_list != ''){
    filter_multiple_peps = sprintf("WHERE PEPTIDE.MODIFIED_SEQUENCE IN (%s)", paste("'" + paste(pep_list, collapse=",") + "'") )
  } else {
    filter_multiple_peps = ''
  }
  
  if ((mscore_filter != '') && (ipf_filter != '')){
    qval_filter_query = sprintf("WHERE SCORE_MS2.QVALUE < %s AND SCORE_IPF.PEP < %s", mscore_filter, ipf_filter)
  } else {
    qval_filter_query = ''
  }
  
  ## filter Out Decoys
  if ( decoy_filter==TRUE ){
    decoy_filter_query = "AND PRECURSOR.DECOY=0"
  } else {
    decoy_filter_query = ""
  }
  
  ## Construct Query Statement
  # print(filter_multiple_peps)
  stmt = sprintf(
    "%s
FROM PRECURSOR
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID %s
%s
%s
%s
LEFT JOIN FEATURE_MS1 ON FEATURE_MS1.FEATURE_ID = FEATURE.ID
LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
--LEFT JOIN SCORE_MS1 ON SCORE_MS1.FEATURE_ID = FEATURE.ID
LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
%s
%s
%s
%s
%s
%s
ORDER BY transition_group_id,
         peak_group_rank
", select_as_stmt, decoy_filter_query, peptide_query, precursor_query, run_id_query, include_ipf_score, pk_grp_rnk_fil_query, filter_multiple_peps, qval_filter_query, mod_peptide_query, mod_residue_position_query)
  
  # cat(stmt)
  
  # Query Databasse
  df_osw <- dplyr::collect( dplyr::tbl(osw_db, dbplyr::sql(stmt)) )
  
  cat("Dimensions of OSW Results file: ", dim(df_osw), "\n")
  
  # Disconnect from database
  DBI::dbDisconnect(osw_db)
  
  return( df_osw )
}