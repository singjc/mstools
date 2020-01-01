#// **********************************************************************************************
#//                         filterSQMASSdb.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Filter and sqMass db file 
#' @description This function can be used to filter an sqMass db file given a list of unmodified sequences
#' 
#' @param sqmass_file A character vector of the absolute path and filename of the sqMass file. (Must be .osw format)
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
filterSQMASSdb <- function( sqmass_file, unmodified_sequence_filter) {
  ## TODO add controls tatements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    sqmass_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/sqmass/test.sqMass"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
  }
  
  ## Setup Logging
  mstools:::log_setup()
  ## Get and Evaluate File Extension Type to ensure an osw file was supplied
  fileType <- tools::file_ext(osw_file)
  if( tolower(fileType)!='sqmass' ){
    MazamaCoreUtils::logger.error( "The supplied file was not a valid OSW database file!\n You provided a file of type: %s", fileType)
  }
  
  ##************************************************
  ##    Establiash Connection to DB
  ##************************************************
  
  db <- DBI::dbConnect( RSQLite::SQLite(), sqmass_file )
  
  ##************************************************
  ##    Filter Precursor Table
  ##************************************************
  ## query statement to get a table of non desired peptide sequences to delete
  precursor_filter_stmt <- sprintf( "SELECT * FROM PRECURSOR
WHERE PRECURSOR.PEPTIDE_SEQUENCE NOT IN ('%s')", paste(unmodified_sequence_filter, collapse="','") )
  
  ## Send query to database
  precursor_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( precursor_filter_stmt )) )
  
  ## Split chromatogram ids into partitions of 100. Large queries don't work well in sqlite
  chromatogram_ids_list <- split(precursor_table$CHROMATOGRAM_ID, ceiling(seq_along(precursor_table$CHROMATOGRAM_ID)/10000))
  
  for ( chromatogram_ids_sub_list in chromatogram_ids_list ){
    
    ## Delete Query
    precursor_delete_stmt <- sprintf( "DELETE FROM PRECURSOR WHERE PRECURSOR.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
    ## Execute delete query
    DBI::dbExecute( db, precursor_delete_stmt )
    
    ##************************************************
    ##    Filter Chromatogram Table
    ##************************************************
    ## chromatogram table delte query
    chromatogram_table_delete_stmt <- sprintf( "DELETE FROM CHROMATOGRAM WHERE CHROMATOGRAM.ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
    ## Execute delete query
    DBI::dbExecute( db, chromatogram_table_delete_stmt )
    
    ##************************************************
    ##    Filter Data Table
    ##************************************************
    ## chromatogram table delte query
    data_table_delete_stmt <- sprintf( "DELETE FROM DATA WHERE DATA.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
    ## Execute delete query
    DBI::dbExecute( db, data_table_delete_stmt )
    
    ##************************************************
    ##    Filter Product Table
    ##************************************************
    ## chromatogram table delte query
    product_table_delete_stmt <- sprintf( "DELETE FROM PRODUCT WHERE PRODUCT.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
    ## Execute delete query
    DBI::dbExecute( db, product_table_delete_stmt )
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
