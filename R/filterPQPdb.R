#// **********************************************************************************************
#//                         filterPQPdb.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Filter and PQP db file 
#' @description This function can be used to filter an PQP db file given a list of unmodified sequences
#' 
#' @param pqp_file A character vector of the absolute path and filename of the sqMass file. (Must be .osw format)
#' @param unmodified_sequence_filter A character vector for extraction of specific peptide(s). I.e. c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
#' @param modified_sequence_filter A character vector for extraction of specific peptide(s) with modifications
#' @return A data.table containing spectral library information
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom DBI dbConnect dbDisconnect dbExecute
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
#' @importFrom tools file_ext
filterPQPdb <- function( pqp_file, unmodified_sequence_filte, modified_sequence_filter=NULL) {
  ## TODO add control statements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    pqp_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/pqp/test.pqp"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
  }
  
  tryCatch(
    expr = {
      
      ## Check if logging has been initialized
      if( !MazamaCoreUtils::logger.isInitialized() ){
        log_setup()
      }
      
      ## Get and Evaluate File Extension Type to ensure an osw file was supplied
      fileType <- tools::file_ext(pqp_file)
      if( tolower(fileType)!='pqp' ){
        MazamaCoreUtils::logger.error( "[mstools::filterPQPdb] The supplied file was not a valid OSW database file!\n You provided a file of type: %s", fileType)
      }
      
      ##************************************************
      ##    Establiash Connection to DB
      ##************************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Connecting to Database: %s", pqp_file)
      db <- DBI::dbConnect( RSQLite::SQLite(), pqp_file )
      
      ##************************************************
      ##    Filter Peptide Table
      ##************************************************
      ## query statement to get a table of only desired unmodified sequences
      if ( is.null(modified_sequence_filter) ) {
      peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.UNMODIFIED_SEQUENCE in ('%s')", paste(unmodified_sequence_filter, collapse="','") )
      } else {
        peptide_filter_stmt <- sprintf( "SELECT * FROM PEPTIDE WHERE PEPTIDE.MODIFIED_SEQUENCE in ('%s')", paste(modified_sequence_filter, collapse="','") )
      }
      ## Send query to database
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", peptide_filter_stmt)
      peptide_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_filter_stmt )) )
      ## Delete Query
      peptide_delete_stmt <- sprintf( "DELETE FROM PEPTIDE WHERE PEPTIDE.ID NOT IN (%s)", paste(peptide_table$ID, collapse = ","))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", peptide_delete_stmt)
      DBI::dbExecute( db, peptide_delete_stmt )
      
      ##***********************************************
      ##    Precursor Filter Selecion
      ##***********************************************
      ## precursor to peptide mapping table query
      precursor_peptide_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PRECURSOR_PEPTIDE_MAPPING on PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", precursor_peptide_mapping_table_stmt)
      precursor_peptide_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( precursor_peptide_mapping_table_stmt ) ) )
      ## Delete Query
      precursor_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", precursor_peptide_mapping_table_delete_stmt)
      DBI::dbExecute( db, precursor_peptide_mapping_table_delete_stmt )
      ## Delete Query for Precursor table
      precursor_table_delete_stmt <- sprintf( "DELETE FROM PRECURSOR WHERE PRECURSOR.ID NOT IN (%s)", paste(precursor_peptide_mapping_table$PRECURSOR_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", precursor_table_delete_stmt)
      DBI::dbExecute( db, precursor_table_delete_stmt ) 
      
      ##***********************************************
      ##    Transition Filter Selection
      ##***********************************************
      ## precursor to peptide mapping table query
      transition_precursor_mapping_table_stmt <- "SELECT * FROM PRECURSOR INNER JOIN TRANSITION_PRECURSOR_MAPPING on TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID=PRECURSOR.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", transition_precursor_mapping_table_stmt)
      transition_precursor_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( transition_precursor_mapping_table_stmt ) ) )
      ## Delete Query
      transition_peptide_mapping_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PEPTIDE_MAPPING WHERE TRANSITION_PEPTIDE_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", transition_peptide_mapping_table_delete_stmt)
      DBI::dbExecute( db, transition_peptide_mapping_table_delete_stmt )
      ## Delete Query for transition to precursor mapping table
      transition_precursor_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = "," ))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", transition_precursor_table_delete_stmt)
      DBI::dbExecute( db, transition_precursor_table_delete_stmt )
      ## Delete query for transition table
      transition_table_delete_stmt <- sprintf( "DELETE FROM TRANSITION WHERE TRANSITION.ID NOT IN (%s)", paste(transition_precursor_mapping_table$TRANSITION_ID, collapse = ","))
      ## Execute delete query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", transition_table_delete_stmt)
      DBI::dbExecute( db, transition_table_delete_stmt )
      
      ##***********************************************
      ##    Protein Filter 
      ##***********************************************
      ## protein to peptide mapping table query
      peptide_protein_mapping_table_stmt <- "SELECT * FROM PEPTIDE INNER JOIN PEPTIDE_PROTEIN_MAPPING on PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID=PEPTIDE.ID"
      ## send query to database to get filtered table
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", peptide_protein_mapping_table_stmt)
      peptide_protein_mapping_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( peptide_protein_mapping_table_stmt ) ) )
      ## Delete query for peptide_protein_protein_mapping
      peptide_protein_mapping_table_delete_stmt <- sprintf( "DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
      ## Execute delte query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", peptide_protein_mapping_table_delete_stmt)
      DBI::dbExecute( db, peptide_protein_mapping_table_delete_stmt )
      ## Delete query for peptide_protein_protein_mapping
      protein_table_delete_stmt <- sprintf( "DELETE FROM PROTEIN WHERE PROTEIN.ID NOT IN (%s)", paste(peptide_protein_mapping_table$PROTEIN_ID, collapse = ",") )
      ## Execute delte query
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Querying Database: %s", protein_table_delete_stmt)
      DBI::dbExecute( db, protein_table_delete_stmt )
      
      
      
      ##***********************************************
      ##    Clear unused space in db
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Vacuuming Database")
      DBI::dbExecute(db, "VACUUM")
      
      ##***********************************************
      ##    Disconnect fom DB 
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterPQPdb] Diconnecting From Database: %s", pqp_file)
      DBI::dbDisconnect( db )
    },
    error = function(e){
      MazamaCoreUtils::logger.error("[mstools::filterPQPdb] filterPQPdb.R: There was the following error that occured during function call...\n", e)
    }
  ) # End tryCatch
} # End Function
