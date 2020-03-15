#// **********************************************************************************************
#//                         getPepLibData.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Extract data from a PQP library file
#' @description This function can be used to extract information from the PQP library file. This can also be used to extract the transition informaiton from and OSW results file.
#' 
#' @param libfile A character vector of the absolute path and filename of the library file. (Must be .pqp format)
#' @param peptide_id A character vector for extraction of a specific peptide. I.e. 'ANSSPTTNIDHLK'
#' @param mod_peptide_id An array of two string vectors indicating a specific modified peptide sequence with both UniMod annotation and actual modification name to extract information for. I.e. c(ANS(Phos)SNSLK, ANS(UniMod:21)SNSLK) (Default: '')
#' @return A data.table containing spectral library information
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
getPepLibData_ <- function( libfile, peptide_id = '', mod_peptide_id=c('','') ){
  
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  tryCatch( 
    expr = {
      # Connect to database
      MazamaCoreUtils::logger.trace(sprintf("[mstools::getPepLibData_] Connecting To Database: %s\n", libfile))
      lib_db <- DBI::dbConnect( RSQLite::SQLite(), libfile )
      # Add Query statement to extract a specific peptide
      if ( peptide_id !='' ){
        peptide_query = sprintf( "INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID AND PEPTIDE.UNMODIFIED_SEQUENCE=('%s')", peptide_id )	
      } else if ( mod_peptide_id[1]!='' & mod_peptide_id[2]!='' ){
        peptide_query = sprintf("INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID AND ( PEPTIDE.MODIFIED_SEQUENCE=('%s') OR PEPTIDE.MODIFIED_SEQUENCE=('%s') )", mod_peptide_id[1], mod_peptide_id[2])
      } else {
        peptide_query = 'INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID'
      }
      # Create Query
      lib_query =  sprintf( "
SELECT TRANSITION.ID AS TRANSITION_ID,
    TRANSITION.TRAML_ID,
    TRANSITION.PRODUCT_MZ,
    TRANSITION.CHARGE,
    TRANSITION.TYPE,
    TRANSITION.ORDINAL,
    TRANSITION.DETECTING,
    TRANSITION.IDENTIFYING,
    TRANSITION.QUANTIFYING,
    TRANSITION.LIBRARY_INTENSITY,
    TRANSITION.DECOY,
    TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID,
    PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID,
    PEPTIDE.MODIFIED_SEQUENCE,
    PEPTIDE.UNMODIFIED_SEQUENCE,
    PRECURSOR.LIBRARY_RT AS PRECURSOR_LIBRARY_RT,
    PRECURSOR.PRECURSOR_MZ,
    PRECURSOR.CHARGE AS PRECURSOR_CHARGE 
FROM TRANSITION
INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION.ID=TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AND TRANSITION.DECOY=0
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID=PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
%s
INNER JOIN PRECURSOR ON PRECURSOR.ID=PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID", peptide_query)
      # Query Databasse
      MazamaCoreUtils::logger.trace(sprintf("[mstools::getPepLibData_] Querying Database: %s\n", lib_query))
      df_lib <- dplyr::collect( dplyr::tbl(lib_db, dbplyr::sql(lib_query)) )
      # Disconnect from database
      MazamaCoreUtils::logger.trace(sprintf("[mstools::getPepLibData_] Disconnecting From Database: %s\n", libfile))
      DBI::dbDisconnect(lib_db)
      return( df_lib )
    },
    error = function(e){
      MazamaCoreUtils::logger.error(sprintf("[mstools::getPepLibData_] There was the following error that occured during function call: %s\n", e$message))
    } )
  
}