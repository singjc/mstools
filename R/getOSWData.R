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
#' @param peptide_id A numeric value for a specific peptide id to extract information for. (Default: '')
#' @param peptide_unmodified A string vector indicating a specific unmodified peptide sequence to extract information for. (Default: '')
#' @param peptide_modified A string vector indicating a specific modified peptide sequence to extract information for. (Default: '')
#' @param mod_peptide_id An array of two string vectors indicating a specific modified peptide sequence with both UniMod annotation and actual modification name to extract information for. I.e. c(ANS(Phos)SNSLK, ANS(UniMod:21)SNSLK) (Default: '') 
#' @param mod_residue_position A numeric value indicating the position of a modification. (Default: '')
#' @param peak_group_rank_filter A Logical value for filtering OSW results by Peak-Group Rank of 1. (Default: FALSE)
#' @param pep_list An arrary of string modified peptide sequences to extract information for. (Default: '')
#' @param mscore_filter A named numeric value to filter results by named q-value. (Default: '') Options: SCORE_MS2, SCORE_IPF, SCORE_PEPTIDE. Example: c(SCORE_MS2=0.01)
#' @param ipf_filter A numeric value to filter results by IPF PEP. (Default: '')
#' @param ipf_score A logical value to extract data using IPF Score results. (Default: FALSE)
#' @param ms2_score A logical value to extract data using MS2 Score results. (Default: TRUE)
#' @param  decoy_filter A logical value to filter decoys out of final results. (Default: TRUE)
#' @return A data.table containing OpenSwath Results information. 
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
getOSWData_ <- function ( oswfile,
                          run_name='',
                          precursor_id='',
                          peptide_id='',
                          peptide_unmodified='',
                          peptide_modified='',
                          mod_peptide_id=c('',''),
                          mod_residue_position='',
                          peak_group_rank_filter=FALSE, 
                          pep_list='',
                          mscore_filter=c(SCORE_=1),
                          ipf_filter='',
                          ipf_score=FALSE,
                          ms2_score=TRUE,
                          decoy_filter=TRUE,
                          inference_level='peptide_query'){
  
  if ( F ) {
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/lower_product_mz_threshold/pyprophet/group_id/merged_runs_group_id_MS1MS2_intergration_ipf.osw"
    # run_name <- '/project/def-hroest/data/synth_phospho_pep/mzML/chludwig_K150309_013_SW_0.mzXML'
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/pyprophet_parametric/merged_runs_group_id_MS1MS2_intergration.osw"
    # run_name <- "/home/singjust/projects/def-hroest/data/christian_doerig_dataset/U_pools_mzML/sw/lgillet_L160915_024-Manchester_dirty_phospho_-_Pool_U12_-_SW.mzML.gz"
    #  run_name = "lgillet_L160915_028-Manchester_dirty_phospho_-_Pool_U14_-_SW.mzML.gz"
    # oswfile <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/M_pools_re-run/pyprophet_parametric_no_precursor/merged_runs_group_id_MS1MS2_intergration.osw"
    # run_name <- "lgillet_L160918_002-Manchester_dirty_phospho_-_Pool_M1_-_SW.mzML.gz"
    # run_name <- "lgillet_L160915_002-Manchester_dirty_phospho_-_Pool_U1_-_SW.mzML.gz"
    # Load Requried Libraries
    # library(dplyr)
    # library(dbplyr)
  }
  
  # Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  # Connect to database
  MazamaCoreUtils::logger.trace(sprintf("[mstools::getOSWData_] Connecting To Database: %s\n", oswfile))
  osw_db <- DBI::dbConnect( RSQLite::SQLite(), oswfile )
  
  ## Filter for a specific run
  if ( run_name != '' ){
    run_id_df = mstools::getRunID_(oswfile, run_name)
    run_id_query = sprintf("AND RUN.ID=(%s)", run_id_df$ID)
  } else {
    run_id_query = ''
  }
  
  ## Check if FEATURE_MS1 table exits
  if(  DBI::dbExistsTable( osw_db, "FEATURE_MS1" ) ){
    join_feature_ms1 <- sprintf("LEFT JOIN FEATURE_MS1 ON FEATURE_MS1.FEATURE_ID = FEATURE.ID")
    select_feature_ms1 <- sprintf("FEATURE_MS1.AREA_INTENSITY AS aggr_prec_Peak_Area,
       FEATURE_MS1.APEX_INTENSITY AS aggr_prec_Peak_Apex,")
  } else {
    join_feature_ms1 <- ""
    select_feature_ms1 <- ""
  }
  
  ## Check is SCORE_MS1 table exists
  if(  DBI::dbExistsTable( osw_db, "SCORE_MS1" ) ){
    join_score_ms1 <- sprintf("LEFT JOIN SCORE_MS1 ON SCORE_MS1.FEATURE_ID = FEATURE.ID")
    select_score_ms1 <- sprintf("SCORE_MS1.PEP AS ms1_pep,")
  } else {
    join_score_ms1 <- ""
    select_score_ms1 <- ""
  }
  
  ## Use IPF scores
  if ( ipf_score ){
    ### Check if SCORE_IPF table exits
    if(  DBI::dbExistsTable( osw_db, "SCORE_IPF" ) ){
      join_score_ipf <- sprintf("INNER JOIN SCORE_IPF ON SCORE_IPF.FEATURE_ID = FEATURE.ID")
      select_score_ipf <- sprintf(",
	   SCORE_IPF.PRECURSOR_PEAKGROUP_PEP AS precursor_pep,
       SCORE_IPF.PEP AS ipf_pep,
       SCORE_IPF.QVALUE AS m_score")
      join_peptide <- sprintf("INNER JOIN PEPTIDE ON ( PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID OR PEPTIDE.ID = SCORE_IPF.PEPTIDE_ID )")
      miscellaneous_score_ipf_control <- "AND SCORE_IPF.QVALUE IS NOT NULL AND PEPTIDE.MODIFIED_SEQUENCE NOT LIKE '%UniMod%'"
    } else {
      MazamaCoreUtils::logger.error(sprintf("[mstools::getOSWData_] There was no SCORE_IPF table found in: %s\nFalling back to SCORE_MS2 scores only", oswfile))
      join_score_ipf <- ""
      select_score_ipf <- ""
      join_peptide <- sprintf("INNER JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID")
      miscellaneous_score_ipf_control <- ""
    }
  } else {
    join_score_ipf <- ""
    select_score_ipf <- ""
    join_peptide <- sprintf("INNER JOIN PEPTIDE ON PEPTIDE.ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID")
    miscellaneous_score_ipf_control <- ""
  }
  
  ## Filter for a specific precursor id
  if ( precursor_id != '' ){
    precursor_query = sprintf("AND PRECURSOR.ID = (%s)", precursor_id)
  } else {
    precursor_query = ''
  }
  
  ## Filter for a specific peptide id
  if( peptide_id != ''){
    peptide_query = sprintf("AND PEPTIDE.ID = (%s)", peptide_id)
  } else {
    peptide_query = ''
  }
  
  ## Filter for a specific unmodified peptide sequence
  if( peptide_unmodified != ''){
    peptide_unmodified_query = sprintf("AND PEPTIDE.UNMODIFIED_SEQUENCE = ('%s')", peptide_unmodified)
  } else {
    peptide_unmodified_query = ''
  }
  
  ## Filter for a specific modified peptide sequence
  if( peptide_modified != ''){
    peptide_modified_query = sprintf("AND PEPTIDE.MODIFIED_SEQUENCE = ('%s')", peptide_modified)
  } else {
    peptide_modified_query = ''
  }
  
  ## Filter Out Decoys
  if ( decoy_filter==TRUE ){
    decoy_filter_query <- "AND PRECURSOR.DECOY=0"
  } else {
    decoy_filter_query <- ""
  }
  
  ## Filter Peak Group Rank
  if ( peak_group_rank_filter ){
    peak_group_rank_filter_query <- "AND SCORE_MS2.RANK=1"
  } else {
    peak_group_rank_filter_query <- ""
  }
  
  ## Filter for a specifc Q_VALUE/M_SCORE
  if ( names(mscore_filter)=="SCORE_MS2" ){
    m_score_filter_query <- sprintf( "AND SCORE_MS2.QVALUE<%s", mscore_filter )
  } else if ( names(mscore_filter)=="SCORE_IPF" ){
    m_score_filter_query <- sprintf( "AND SCORE_IPF.QVALUE<%s", mscore_filter )
  } else if ( names(mscore_filter)=="SCORE_PEPTIDE" ){
    m_score_filter_query <- sprintf( "AND SCORE_PEPTIDE.QVALUE<%s", mscore_filter )
  } else {
    m_score_filter_query <- ""
  }
  
  if ( mod_peptide_id[1]!='' & mod_peptide_id[2]!='' ){
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
  if ( inference_level=="peptide_inference" ){
    inference_level_query <- "INNER JOIN SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID=PEPTIDE.ID"
    select_level_query <- sprintf(", SCORE_PEPTIDE.CONTEXT as context,
                                  SCORE_PEPTIDE.PEP as peptide_pep,
                                  SCORE_PEPTIDE.SCORE as peptide_d_score,
                                  SCORE_PEPTIDE.QVALUE as peptide_m_score")
  } else if ( inference_level=="protein_inference" ){
    inference_level_query <- "INNER JOIN SCORE_PROTEIN ON SCORE_PROTEIN.PEPTIDE_ID=PEPTIDE.ID"
  } else {
    inference_level_query <- ""
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
       %s
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       %s
       SCORE_MS2.PEP AS ms2_pep,
       SCORE_IPF.PRECURSOR_PEAKGROUP_PEP AS precursor_pep,
       SCORE_IPF.PEP AS ipf_pep,
       SCORE_MS2.RANK AS peak_group_rank,
       SCORE_MS2.SCORE AS d_score,
       SCORE_MS2.QVALUE AS ms2_m_score,
       SCORE_IPF.QVALUE AS m_score", select_feature_ms1, select_score_ms1 )
    include_ipf_score = 'LEFT JOIN SCORE_IPF ON SCORE_MS2.FEATURE_ID = SCORE_IPF.FEATURE_ID'
    
    if (peak_group_rank_filter == TRUE){
      pk_grp_rnk_fil_query = 'INNER JOIN PEPTIDE AS PEPTIDE_IPF ON SCORE_IPF.PEPTIDE_ID = PEPTIDE_IPF.ID '
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
       %s
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       %s
       SCORE_MS2.PEP AS ms2_pep,
       SCORE_MS2.RANK AS peak_group_rank,
       SCORE_MS2.SCORE AS d_score,
       SCORE_MS2.QVALUE AS ms2_m_score", select_feature_ms1, select_score_ms1 )
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
  
  
  
  ## Construct Query Statement
  #   stmt = sprintf(
  #     "%s 
  # FROM PRECURSOR
  # INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID %s
  # %s
  # %s
  # %s
  # %s -- Join FEATURE_MS1 table
  # LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
  # %s -- Join SCORE_MS1
  # LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
  # %s
  # %s
  # %s
  # %s
  # %s
  # %s
  # ORDER BY transition_group_id,
  #          peak_group_rank
  # ", select_as_stmt, decoy_filter_query, peptide_query, precursor_query, run_id_query, join_feature_ms1, join_score_ms1, include_ipf_score, pk_grp_rnk_fil_query, filter_multiple_peps, qval_filter_query, mod_peptide_query, mod_residue_position_query)
  #   
  stmt = sprintf(
    "
     SELECT 
FEATURE.ID AS feature_id,
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
       PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
       PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
       PRECURSOR.CHARGE AS Charge,
       PRECURSOR.PRECURSOR_MZ AS mz,
       FEATURE_MS2.AREA_INTENSITY AS Intensity,
       %s -- FEATURE_MS1 select statements #select_feature_ms1
       FEATURE.LEFT_WIDTH AS leftWidth,
       FEATURE.RIGHT_WIDTH AS rightWidth,
       %s -- SCORE_MS1 select statements #select_score_ms1
       SCORE_MS2.PEP AS ms2_pep,
       SCORE_MS2.RANK AS peak_group_rank,
       SCORE_MS2.SCORE AS d_score,
       SCORE_MS2.QVALUE AS ms2_m_score
       %s -- SCORE_IPF select statements #select_score_ipf
FROM FEATURE
INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
%s -- Join FEATURE_MS1 Table if available #join_feature_ms1
INNER JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
%s -- Join SCORE_MS1 Table if available #join_score_ms1
INNER JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
%s -- Join SCORE_IPF Table if available #join_score_ipf
INNER JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID %s -- Filter for Decoyrs #decoy_filter_query
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID = FEATURE.PRECURSOR_ID
%s -- Join PEPTIDE Table #join_peptide
WHERE FEATURE.ID IS NOT NULL -- Default WHERE Being clause
%s -- Miscellaneous control statements for if SCORE_IPF is used #miscellaneous_score_ipf_control
%s -- Filter for a specific PRECURSOR.ID #precursor_query
%s -- Filter for a specific PEPTIDE.ID #peptide_query
%s -- Filter for a specific unmodified peptide sequence #peptide_unmodified_query
%s -- Filter for a specific modified peptide sequence #peptide_modified_query
%s -- Filter for specific RUN.ID #run_id_query
%s -- Filter for PeakGroupRank=1 #peak_group_rank_filter_query
%s -- Filter for a specific level of Q-Value. (SCORE_MS2, SCORE_IPF, SCORE_PEPTIDE)
ORDER BY transition_group_id,
         peak_group_rank
", select_feature_ms1, select_score_ms1, select_score_ipf, join_feature_ms1, join_score_ms1, join_score_ipf, decoy_filter_query, join_peptide, 
    miscellaneous_score_ipf_control, precursor_query, peptide_query, peptide_unmodified_query, peptide_modified_query, run_id_query, peak_group_rank_filter_query, m_score_filter_query )
  
  # Query Databasse
  MazamaCoreUtils::logger.trace(sprintf("[mstools::getOSWData_] QueryingDatabase: %s\n", stmt))
  df_osw <- dplyr::collect( dplyr::tbl(osw_db, dbplyr::sql(stmt)) )
  
  ## Get Precursor Table for Precursor_Peptide_Mapping
  stmt2 <- sprintf("SELECT 
PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID AS id_peptide,
PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AS id_precursor,
PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
PEPTIDE.DECOY AS decoy,
PRECURSOR.TRAML_ID AS traml_id,
PRECURSOR.GROUP_LABEL AS group_label,
PRECURSOR.PRECURSOR_MZ AS precursor_mz,
PRECURSOR.CHARGE AS precursor_charge,
PRECURSOR.LIBRARY_INTENSITY AS library_int,
PRECURSOR.LIBRARY_RT AS library_rt
FROM PEPTIDE
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
INNER JOIN PRECURSOR ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID" )
  df_precursor_peptide_mapping_table <- dplyr::collect( dplyr::tbl(osw_db, dbplyr::sql(stmt2)) )
  
  ## Get Peptide Table
  df_peptide_table <- dplyr::collect( dplyr::tbl(osw_db, dbplyr::sql("SELECT * FROM PEPTIDE")) )
  
  ## TODO: Need to make this more robust in the case that TRAML_ID is no longer valid
  # df_osw %>%
  #   dplyr::filter( paste(df_osw$FullPeptideName, df_osw$id_precursor, sep="_") %in% unlist(lapply(seq(1:dim(df_precursor_peptide_mapping_table)[1]), function(x){ gsub("*\\d+$", df_precursor_peptide_mapping_table$ID[x], df_precursor_peptide_mapping_table$TRAML_ID[x]) } )) ) -> df_osw
  # 
  df_osw$FullPeptideName.unimod <- unlist(lapply(df_osw$FullPeptideName, mstools::codenameTounimod))
  
  df_osw %>%
    dplyr::filter( paste(df_osw$FullPeptideName.unimod, df_osw$id_precursor, sep='_') %in% paste(df_precursor_peptide_mapping_table$FullPeptideName, df_precursor_peptide_mapping_table$id_precursor, sep='_') ) -> df_osw
  ## Remove tmp column
  df_osw$FullPeptideName.unimod <- NULL
  
  ## Context Level Inference
  
  if ( any(inference_level %in% "peptide_inference") ){
    inference_level_query <- sprintf( "SELECT
                                        SCORE_PEPTIDE.CONTEXT AS context,
                                        SCORE_PEPTIDE.RUN_ID AS run_id,
                                        SCORE_PEPTIDE.PEPTIDE_ID AS id_peptide,
                                        PEPTIDE.UNMODIFIED_SEQUENCE AS Sequence,
                                        PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName,
                                        PEPTIDE.DECOY as decoy,
                                        SCORE_PEPTIDE.PEP AS peptide_pep,
                                        SCORE_PEPTIDE.SCORE AS peptide_d_score,
                                        SCORE_PEPTIDE.QVALUE AS peptide_m_score
                                        FROM PEPTIDE
                                        INNER JOIN SCORE_PEPTIDE ON PEPTIDE.ID=SCORE_PEPTIDE.PEPTIDE_ID" )
    
    df_peptide_inference <- dplyr::collect( dplyr::tbl(osw_db, dbplyr::sql(inference_level_query)) )
    
    if ( unique(df_peptide_inference$context)=='global' ){
      MazamaCoreUtils::logger.trace(sprintf("[mstools::getOSWData_] Results contain GLOBAL level contexts.. Returning single data.frame for GLOBAL level contexts"))
      ## Check if remove decoys is set to true
      if ( decoy_filter==TRUE ){
        df_peptide_inference %>%
          dplyr::filter( decoy==0 ) -> df_peptide_inference
      }
      df_osw <- df_peptide_inference
    } else if ( unique(df_peptide_inference$context)=='experiment-wide' ){
      ## Check if we need to filter for a specific run id
      if ( run_name != '' ){
        df_peptide_inference %>%
          dplyr::filter( run_id==run_id_df$ID ) -> df_peptide_inference
      }
      ## Check if remove decoys is set to true
      if ( decoy_filter==TRUE ){
        df_peptide_inference %>%
          dplyr::filter( decoy==0 ) -> df_peptide_inference
      }
      
      df_peptide_inference$FullPeptideName <- unlist(lapply( df_peptide_inference$FullPeptideName, mstools::unimodTocodename )) # TODO: Need to make unimodTocodename work without doing a list loop
      
      df_peptide_inference <- merge(df_peptide_inference, dplyr::select(df_peptide_table, c("ID", "MODIFIED_SEQUENCE")), by.x="FullPeptideName", by.y="MODIFIED_SEQUENCE", all.x = T)
      df_osw <- merge(df_osw, dplyr::select(df_peptide_inference, c("FullPeptideName", "ID", "context", "peptide_pep", "peptide_d_score", "peptide_m_score")), by.x=c("FullPeptideName", "id_peptide"), by.y = c("FullPeptideName", "ID"), all.x=T)
      
    } else if ( unique(df_peptide_inference$context)=='run-specific' ){
      MazamaCoreUtils::logger.error( sprintf("[mstools::getOSWData_] run-specific context level is not yet supported... Returning Default MS2_Score/IPF_Score results.." ) )
    } else {
      MazamaCoreUtils::logger.error( sprintf("[mstools::getOSWData_] Unknown '%s' context level... Returning Default MS2_Score/IPF_Score results..",  unique(df_peptide_inference$context) ) )
    }
    
  }
  
  if ( any(inference_level %in% "protein_inference") ){
    MazamaCoreUtils::logger.trace("[mstools::getOSWData_] Protein level inference is not yet supported. Submit Feature Request on GitHub to singjc/mstools")
  } 
  
  MazamaCoreUtils::logger.trace(paste("[mstools::getOSWData_] Dimensions of OSW Results file: ", dim(df_osw), "\n"))
  
  # Disconnect from database
  MazamaCoreUtils::logger.trace(sprintf("[mstools::getOSWData_] Disconnecting From Database: %s\n", oswfile))
  DBI::dbDisconnect(osw_db)
  
  return( df_osw )
}