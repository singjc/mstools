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
filterPQPdb <- function( pqp_file, unmodified_sequence_filter) {
  ## TODO add controls tatements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    pqp_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/pqp/test.pqp"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
  }
  
  # ## Setup Logging
  # mstools:::log_setup()
  # ## Get and Evaluate File Extension Type to ensure an osw file was supplied
  # fileType <- tools::file_ext(pqp_file)
  # if( tolower(fileType)!='pqp' ){
  #   MazamaCoreUtils::logger.error( "The supplied file was not a valid OSW database file!\n You provided a file of type: %s", fileType)
  # }
  
  ##************************************************
  ##    Establiash Connection to DB
  ##************************************************
  
  db <- DBI::dbConnect( RSQLite::SQLite(), pqp_file )
  
  ##************************************************
  ##    Filter Peptide Table
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
  ##    Clear unused space in db
  ##***********************************************
  DBI::dbExecute(db, "VACUUM")
  
  ##***********************************************
  ##    Disconnect fom DB 
  ##***********************************************
  DBI::dbDisconnect( db )
}
