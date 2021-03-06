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
#' @param modified_sequence_filter A character vector for extraction of specific peptide(s) with modifications
#' @param random_seed (numeric) Set the seed for sampling
#' @param ratio_keep (float) Set the fraction of sequences to keep
#' @return A data.table containing spectral library information
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbExistsTable dbSendQuery dbBind
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
#' @importFrom tools file_ext
filterOSWdb <- function( osw_file, unmodified_sequence_filter=NULL, modified_sequence_filter=NULL, random_seed=NULL, ratio_keep=NULL) {
  ## TODO add controls tatements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    osw_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/osw/merged.osw"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
    osw_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/George_lib_repeat2/pyprophet/subsample_test2/test_13_2.osw"
    osw_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR_Test_Data/Subsample_Test/test_13.osw"
    random_seed <- 1
    ratio_keep <- 0.5
  }
  
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  tryCatch(
    expr = {

      ## Get and Evaluate File Extension Type to ensure an osw file was supplied
      fileType <- tools::file_ext(osw_file)
      if( tolower(fileType)!='osw' ){
        MazamaCoreUtils::logger.error( "[mstools::filterOSWdb] The supplied file was not a valid OSW database file!\n You provided a file of type: %s", fileType)
      }
      
      ##************************************************
      ##    Establiash Connection to DB
      ##************************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Connecting to Database: %s", osw_file)
      db <- DBI::dbConnect( RSQLite::SQLite(), osw_file )
      
      
      
      ##************************************************
      ##    Peptide Filter Selection
      ##************************************************
      ## query statement to get a table of only desired unmodified sequences
      if ( !is.null(unmodified_sequence_filter) & is.null(modified_sequence_filter) ) {
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Using Unmodified Sequences Input")
        peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.UNMODIFIED_SEQUENCE in ('%s')", paste(unmodified_sequence_filter, collapse="','") )
      } else if ( is.null(unmodified_sequence_filter) & !is.null(modified_sequence_filter) ) {
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Using Modified Sequences Input")
        peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.MODIFIED_SEQUENCE in ('%s')", paste(modified_sequence_filter, collapse="','") )
      } else if (  is.null(unmodified_sequence_filter) & is.null(modified_sequence_filter) & !is.null(random_seed) & !is.null(ratio_keep) ){
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Using random seed (%s) selection for fraction %s", random_seed, ratio_keep)
        tmp_peptide_table_query <- "SELECT 
PEPTIDE.ID AS PEPTIDE_ID,
PEPTIDE.UNMODIFIED_SEQUENCE AS UNMODIFIED_SEQUENCE,
PEPTIDE.MODIFIED_SEQUENCE AS MODIFIED_SEQUENCE,
PEPTIDE.DECOY AS PEPTIDE_DECOY,
PRECURSOR.ID AS PRECURSOR_ID,
PRECURSOR.PRECURSOR_MZ AS PRECURSOR_MZ,
PRECURSOR.CHARGE AS PRECURSOR_CHARGE,
PRECURSOR.LIBRARY_RT as PRECURSOR_LIBRARY_RT,
PRECURSOR.DECOY AS PRECURSOR_DECOY
FROM PEPTIDE
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
INNER JOIN PRECURSOR ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID"
        tmp_peptide_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( tmp_peptide_table_query )) )
        
        tmp_peptide_table %>%
          dplyr::select( PEPTIDE_ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, PRECURSOR_DECOY ) %>%
          unique() -> write_peptide_table
        colnames(write_peptide_table)[c(1,4)] <- c("ID", "DECOY")
        
        update <- DBI::dbSendQuery(db, 'update PEPTIDE set "ID"=?, "UNMODIFIED_SEQUENCE"=?, "MODIFIED_SEQUENCE"=?, "DECOY"=?  WHERE _rowid_=?')
        
        DBI::dbWriteTable(db, "PEPTIDE", write_peptide_table, overwrite=TRUE)  # send the updated data
        
        ## Set a Random Seed
        set.seed( random_seed )
        ## Subsample Targets
        tmp_peptide_table %>%
          dplyr::filter( PRECURSOR_DECOY==0 ) %>%
          dplyr::sample_frac( size = ratio_keep ) -> tmp_peptide_table_subsampled
        
        ## Get Decoy Counter Parts
        tmp_peptide_table %>%
          dplyr::filter( PRECURSOR_DECOY==1 ) %>%
          dplyr::filter( PRECURSOR_MZ %in% tmp_peptide_table_subsampled$PRECURSOR_MZ & PRECURSOR_CHARGE %in% tmp_peptide_table_subsampled$PRECURSOR_CHARGE & PRECURSOR_LIBRARY_RT %in% tmp_peptide_table_subsampled$PRECURSOR_LIBRARY_RT ) -> tmp_peptide_table_subsampled_decoys
        
        tmp_peptide_table_subsampled <- data.table::rbindlist(list(tmp_peptide_table_subsampled, tmp_peptide_table_subsampled_decoys))
        
        peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.MODIFIED_SEQUENCE in ('%s')", paste(tmp_peptide_table_subsampled$MODIFIED_SEQUENCE, collapse="','") )
        ## Disconnect from database
        DBI::dbDisconnect( db )
        ## Reconnect to database to reflect changes
        db <- DBI::dbConnect( RSQLite::SQLite(), osw_file )
      }
      ## Send query to database
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", peptide_filter_stmt)
      peptide_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_filter_stmt )) )
      ## Delete Query
      peptide_delete_stmt <- sprintf( "DELETE FROM PEPTIDE WHERE PEPTIDE.ID NOT IN (%s)", paste(peptide_table$ID, collapse = ","))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", peptide_delete_stmt)
      DBI::dbExecute( db, peptide_delete_stmt )
      
      ##***********************************************
      ##    Precursor Filter Selecion
      ##***********************************************
      ## precursor to peptide mapping table query
      precursor_peptide_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PRECURSOR_PEPTIDE_MAPPING on PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", precursor_peptide_mapping_table_stmt)
      precursor_peptide_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( precursor_peptide_mapping_table_stmt ) ) )
      ## Delete Query
      precursor_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", precursor_peptide_mapping_table_delete_stmt)
      DBI::dbExecute( db, precursor_peptide_mapping_table_delete_stmt )
      ## Delete Query for Precursor table
      precursor_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR WHERE PRECURSOR.ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", precursor_table_delete_stmt)
      DBI::dbExecute( db, precursor_table_delete_stmt )
      
      ##***********************************************
      ##    Transition Filter Selection
      ##***********************************************
      ## precursor to peptide mapping table query
      transition_precursor_mapping_table_stmt <- "SELECT * FROM PRECURSOR INNER JOIN TRANSITION_PRECURSOR_MAPPING on TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID=PRECURSOR.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", transition_precursor_mapping_table_stmt)
      transition_precursor_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( transition_precursor_mapping_table_stmt ) ) )
      ## Delete Query
      transition_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PEPTIDE_MAPPING WHERE TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", transition_peptide_mapping_table_delete_stmt)
      DBI::dbExecute( db, transition_peptide_mapping_table_delete_stmt )
      ## Delete Query for transition to precursor mapping table
      transition_precursor_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", transition_precursor_table_delete_stmt)
      DBI::dbExecute( db, transition_precursor_table_delete_stmt )
      ## Delete query for transition table
      transition_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION WHERE TRANSITION.ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = ","))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", transition_table_delete_stmt)
      DBI::dbExecute( db, transition_table_delete_stmt )
      
      ##***********************************************
      ##    Protein Filter 
      ##***********************************************
      ## protein to peptide mapping table query
      peptide_protein_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PEPTIDE_PROTEIN_MAPPING on PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", peptide_protein_mapping_table_stmt)
      peptide_protein_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_protein_mapping_table_stmt ) ) )
      ## Delete query for peptide_protein_protein_mapping
      peptide_protein_mapping_table_delete_stmt <- sprintf( "DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
      ## Execute delte query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", peptide_protein_mapping_table_delete_stmt)
      DBI::dbExecute( db, peptide_protein_mapping_table_delete_stmt )
      ## Delete query for peptide_protein_protein_mapping
      protein_table_delete_stmt <- sprintf( "DELETE FROM PROTEIN WHERE PROTEIN.ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
      ## Execute delte query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", protein_table_delete_stmt)
      DBI::dbExecute( db, protein_table_delete_stmt )
      
      ##***********************************************
      ##    Feature Filter 
      ##***********************************************
      ## feature table query
      feature_table_stmt <- "SELECT * FROM FEATURE INNER JOIN PRECURSOR ON PRECURSOR.ID=FEATURE.PRECURSOR_ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", feature_table_stmt)
      feature_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( feature_table_stmt ) ) )
      ## delete query for feature table
      feature_table_delete_stmt <- sprintf( "DELETE FROM FEATURE WHERE FEATURE.ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
      ## Exectute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", feature_table_delete_stmt)
      DBI::dbExecute( db, feature_table_delete_stmt )
      
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ###     FEATURE_MS1
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if ( DBI::dbExistsTable( db, "FEATURE_MS1" ) ){
        ## delete query for feature ms1 table
        feature_ms1_table_delete_stmt <- sprintf( "DELETE FROM FEATURE_MS1 WHERE FEATURE_MS1.FEATURE_ID NOT IN (%s)", paste(feature_table$ID, collapse = ","))
        ## Exectute delete query
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", feature_ms1_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", feature_ms2_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", feature_transition_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", score_ipf_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", score_ms1_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", score_ms2_table_delete_stmt)
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
        MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Querying Database: %s", score_transition_table_delete_stmt)
        DBI::dbExecute( db, score_transition_table_delete_stmt )
      } else {
        MazamaCoreUtils::logger.warn( "There was no SCORE_TRANSITION table." )
      }
      
      ##***********************************************
      ##    Clear unused space in db
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Vacuuming Database")
      DBI::dbExecute(db, "VACUUM")
      
      ##***********************************************
      ##    Disconnect fom DB 
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterOSWdb] Disconnecting From Database: %s", osw_file)
      DBI::dbDisconnect( db )
    },
    error = function(e){
      MazamaCoreUtils::logger.error(sprintf("[mstools::filterOSWdb] There was the following error that occured during function call: %s\n", e$message))
    }
  ) # End tryCatch
} # End function