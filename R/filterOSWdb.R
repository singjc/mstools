#// **********************************************************************************************
#//                         filterOSWdb.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Filter and OSW db file 
#' @description This function can be used to filter an osw db file given a list of unmodified sequences
#' 
#' @param osw_file A character vector of the absolute path and filename of the osw file. (Must be .osw format)
#' @param unmodified_sequence_filter A character vector for extraction of specific peptide(s). I.e. c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
#' @return A data.table containing spectral library information
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom tools file_ext
#' @import MazamaCoreUtils
#' @import DBI
#' @import RSQLite
#' @import dplyr
#' @import dbplyr
filterOSWdb <- function( osw_file, unmodified_sequence_filter) {
  ## TODO add controls tatements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    osw_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/osw/merged.osw"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
  }
  
  # ## Setup Logging
  # mstools:::log_setup()
  # ## Get and Evaluate File Extension Type to ensure an osw file was supplied
  # fileType <- tools::file_ext(osw_file)
  # if( tolower(fileType)!='osw' ){
  #   MazamaCoreUtils::logger.error( "The supplied file was not a valid OSW database file!\n You provided a file of type: %s", fileType)
  # }
  
  ##************************************************
  ##    Establiash Connection to DB
  ##************************************************
  
  db <- DBI::dbConnect( RSQLite::SQLite(), osw_file )
  
  ##************************************************
  ##    Peptide Filter Selection
  ##************************************************
  ## query statement to get a table of only desired unmodified sequences
  peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.UNMODIFIED_SEQUENCE in ('%s')", paste(unmodified_sequence_filter, collapse="','") )
  ## Send query to database
  peptide_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_filter_stmt )) )
  ## Delete Query
  peptide_delete_stmt <- sprintf( "DELETE FROM PEPTIDE WHERE PEPTIDE.ID NOT IN (%s)", paste(peptide_table$ID, collapse = ","))
  ## Execute delete query
  DBI::dbExecute( db, peptide_delete_stmt )
  
  ##***********************************************
  ##    Precursor Filter Selecion
  ##***********************************************
  ## precursor to peptide mapping table query
  precursor_peptide_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PRECURSOR_PEPTIDE_MAPPING on PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
  ## send query to database to get filtered table
  precursor_peptide_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( precursor_peptide_mapping_table_stmt ) ) )
  ## Delete Query
  precursor_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
  ## Execute delete query
  DBI::dbExecute( db, precursor_peptide_mapping_table_delete_stmt )
  ## Delete Query for Precursor table
  precursor_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR WHERE PRECURSOR.ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
  ## Execute delete query
  DBI::dbExecute( db, precursor_table_delete_stmt )
  
  ##***********************************************
  ##    Transition Filter Selection
  ##***********************************************
  ## precursor to peptide mapping table query
  transition_precursor_mapping_table_stmt <- "SELECT * FROM PRECURSOR INNER JOIN TRANSITION_PRECURSOR_MAPPING on TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID=PRECURSOR.ID"
  ## send query to database to get filtered table
  transition_precursor_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( transition_precursor_mapping_table_stmt ) ) )
  ## Delete Query
  transition_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PEPTIDE_MAPPING WHERE TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
  ## Execute delete query
  DBI::dbExecute( db, transition_peptide_mapping_table_delete_stmt )
  ## Delete Query for transition to precursor mapping table
  transition_precursor_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
  ## Execute delete query
  DBI::dbExecute( db, transition_precursor_table_delete_stmt )
  ## Delete query for transition table
  transition_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION WHERE TRANSITION.ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = ","))
  ## Execute delete query
  DBI::dbExecute( db, transition_table_delete_stmt )
  
  ##***********************************************
  ##    Protein Filter 
  ##***********************************************
  ## protein to peptide mapping table query
  peptide_protein_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PEPTIDE_PROTEIN_MAPPING on PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
  ## send query to database to get filtered table
  peptide_protein_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_protein_mapping_table_stmt ) ) )
  ## Delete query for peptide_protein_protein_mapping
  peptide_protein_mapping_table_delete_stmt <- sprintf( "DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
  ## Execute delte query
  DBI::dbExecute( db, peptide_protein_mapping_table_delete_stmt )
  ## Delete query for peptide_protein_protein_mapping
  protein_table_delete_stmt <- sprintf( "DELETE FROM PROTEIN WHERE PROTEIN.ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
  ## Execute delte query
  DBI::dbExecute( db, protein_table_delete_stmt )
  
  ##***********************************************
  ##    Feature Filter 
  ##***********************************************
  ## feature table query
  feature_table_stmt <- "SELECT * FROM FEATURE INNER JOIN PRECURSOR ON PRECURSOR.ID=FEATURE.PRECURSOR_ID"
  ## send query to database to get filtered table
  feature_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( feature_table_stmt ) ) )
  ## delete query for feature table
  feature_table_delete_stmt <- sprintf( "DELETE FROM FEATURE WHERE FEATURE.ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
  ## Exectute delete query
  DBI::dbExecute( db, feature_table_delete_stmt )
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     FEATURE_MS1
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "FEATURE_MS1" ) ){
    ## delete query for feature ms1 table
    feature_ms1_table_delete_stmt <- sprintf( "DELETE FROM FEATURE_MS1 WHERE FEATURE_MS1.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, feature_ms1_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no FEATURE_MS1 table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     FEATURE_MS2
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "FEATURE_MS2" ) ){
    ## delete query for feature ms2 table
    feature_ms2_table_delete_stmt <- sprintf( "DELETE FROM FEATURE_MS2 WHERE FEATURE_MS2.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, feature_ms2_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no FEATURE_MS2 table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     FEATURE_TRANSITION
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "FEATURE_TRANSITION" ) ){
    ## delete query for feature transition table
    feature_transition_table_delete_stmt <- sprintf( "DELETE FROM FEATURE_TRANSITION WHERE FEATURE_TRANSITION.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, feature_transition_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no FEATURE_TRANSITION table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     SCORE_IPF
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "SCORE_IPF" ) ){
    ## delete query for score ipf table
    score_ipf_table_delete_stmt <- sprintf( "DELETE FROM SCORE_IPF WHERE SCORE_IPF.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, score_ipf_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no SCORE_IPF table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     SCORE_MS1
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "SCORE_MS1" ) ){
    ## delete query for score ms1 table
    score_ms1_table_delete_stmt <- sprintf( "DELETE FROM SCORE_MS1 WHERE SCORE_MS1.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, score_ms1_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no SCORE_MS1 table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     SCORE_MS2
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "SCORE_MS2" ) ){
    ## delete query for score ms2 table
    score_ms2_table_delete_stmt <- sprintf( "DELETE FROM SCORE_MS2 WHERE SCORE_MS2.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, score_ms2_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no SCORE_MS2 table." )
  }
  
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###     SCORE_TRANSITION
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if ( DBI::dbExistsTable( db, "SCORE_TRANSITION" ) ){
    ## delete query for score transition table
    score_transition_table_delete_stmt <- sprintf( "DELETE FROM SCORE_TRANSITION WHERE SCORE_TRANSITION.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
    ## Exectute delete query
    DBI::dbExecute( db, score_transition_table_delete_stmt )
  } else {
    MazamaCoreUtils::logger.warn( "There was no SCORE_TRANSITION table." )
  }
  
  ##***********************************************
  ##    Clear unused space in db
  ##***********************************************
  DBI::dbExecute(db, "VACUUM")
  
  ##***********************************************
  ##    Disconnect fom DB 
  ##***********************************************
  DBI::dbDisconnect( db )
}