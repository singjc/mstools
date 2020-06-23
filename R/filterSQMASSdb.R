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
#' @importFrom DBI dbConnect dbDisconnect dbExecute
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn logger.trace
#' @importFrom tools file_ext
filterSQMASSdb <- function( sqmass_file, unmodified_sequence_filter, recreate_indexes=FALSE) {
  ## TODO add controls tatements for check tables being present
  DEBUG=FALSE
  if ( DEBUG ){
    sqmass_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/Synthetic_Dilution_Phosphoproteomics/sqmass/test.sqMass"
    unmodified_sequence_filter <- c('ANSSPTTNIDHLK', 'ESTAEPDSLSR', 'NLSPTKQNGKATHPR', 'KDSNTNIVLLK', 'NKESPTKAIVR')
    sqmass_file <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/test_data/Synthetic_Dilution_Phosphoproteomics/sqmass/test.sqMass"
    
    
    library(mstools)
    ##*********************************************
    ##      Python Setup
    ##*********************************************
    
    if ( !reticulate::py_available()  ){
      find_python()
      install_python_dependencies()
      MazamaCoreUtils::logger.trace( "[mstools::getChromatogramDataPoints_] ** Loading Python Modules **")
      .onload()
    }
    
    db <- DBI::dbConnect( RSQLite::SQLite(), sqmass_file )
  
      
    meta <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( "select * from RUN_EXTRA" )) )
    
    tmp <- pyzlib$decompress(  pybuiltins$bytes( reticulate::r_to_py(meta$DATA[[1]]) ) )
    
    # tmp2 <- reticulate::py_to_r( tmp )
    
    h <- xml2::read_html( tmp )
    
    tmp_outfile <- tempfile( pattern = "meta_data", tmpdir = getwd(), fileext = ".xml")
    xml2::write_xml(h, tmp_outfile, options = "format")
    
  }
  
  tryCatch(
    expr = {
      
      ## Check if logging has been initialized
      if( !MazamaCoreUtils::logger.isInitialized() ){
        log_setup()
      }
      
      ## Get and Evaluate File Extension Type to ensure an osw file was supplied
      fileType <- tools::file_ext(sqmass_file)
      if( tolower(fileType)!='sqmass' ){
        MazamaCoreUtils::logger.error( "[mstools::filterSQMASSdb] The supplied file was not a valid SQMASS database file!\n You provided a file of type: %s", fileType)
      }
      
      ##************************************************
      ##    Establiash Connection to DB
      ##************************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Connecting to Database: %s", sqmass_file)
      db <- DBI::dbConnect( RSQLite::SQLite(), sqmass_file )
      
      ##************************************************
      ##    Filter Precursor Table
      ##************************************************
      ## query statement to get a table of non desired peptide sequences to delete
      precursor_filter_stmt <- sprintf( "SELECT * FROM PRECURSOR
WHERE PRECURSOR.PEPTIDE_SEQUENCE NOT IN ('%s')", paste(unmodified_sequence_filter, collapse="','") )
      
      ## Send query to database
      MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Querying Database: %s", precursor_filter_stmt)
      precursor_table <- dplyr::collect( dplyr::tbl( db, dbplyr::sql( precursor_filter_stmt )) )
      
      ## Split chromatogram ids into partitions of 100. Large queries don't work well in sqlite
      chromatogram_ids_list <- split(precursor_table$CHROMATOGRAM_ID, ceiling(seq_along(precursor_table$CHROMATOGRAM_ID)/10000))
      
      for ( chromatogram_ids_sub_list in chromatogram_ids_list ){
        
        ## Delete Query
        precursor_delete_stmt <- sprintf( "DELETE FROM PRECURSOR WHERE PRECURSOR.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
        ## Execute delete query
        MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Querying Database: %s", precursor_delete_stmt)
        DBI::dbExecute( db, precursor_delete_stmt )
        
        ##************************************************
        ##    Filter Chromatogram Table
        ##************************************************
        ## chromatogram table delte query
        chromatogram_table_delete_stmt <- sprintf( "DELETE FROM CHROMATOGRAM WHERE CHROMATOGRAM.ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
        ## Execute delete query
        MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Querying Database: %s", chromatogram_table_delete_stmt)
        DBI::dbExecute( db, chromatogram_table_delete_stmt )
        
        ##************************************************
        ##    Filter Data Table
        ##************************************************
        ## chromatogram table delte query
        data_table_delete_stmt <- sprintf( "DELETE FROM DATA WHERE DATA.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
        ## Execute delete query
        MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Querying Database: %s", data_table_delete_stmt)
        DBI::dbExecute( db, data_table_delete_stmt )
        
        ##************************************************
        ##    Filter Product Table
        ##************************************************
        ## chromatogram table delte query
        product_table_delete_stmt <- sprintf( "DELETE FROM PRODUCT WHERE PRODUCT.CHROMATOGRAM_ID IN (%s)", paste(chromatogram_ids_sub_list, collapse = ","))
        ## Execute delete query
        MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Querying Database: %s", product_table_delete_stmt)
        DBI::dbExecute( db, product_table_delete_stmt )
      }
      
      if ( recreate_indexes ){
        ##***********************************************
        ##    Drop and Recreate Indexes
        ##***********************************************
        
        drop_idxs <- "DROP INDEX 'main'.'chrom_run_idx';
      DROP INDEX IF EXISTS 'main'.'data_chr_idx';
      DROP INDEX IF EXISTS 'main'.'data_sp_idx';
      DROP INDEX IF EXISTS 'main'.'precursor_chr_idx';
      DROP INDEX IF EXISTS 'main'.'precursor_sp_idx';
      DROP INDEX IF EXISTS 'main'.'product_chr_idx';
      DROP INDEX IF EXISTS 'main'.'product_sp_idx';
      DROP INDEX IF EXISTS 'main'.'run_extra_idx';
      DROP INDEX IF EXISTS 'main'.'spec_mslevel_idx';
      DROP INDEX IF EXISTS 'main'.'spec_rt_idx';
      DROP INDEX IF EXISTS 'main'.'spec_run_idx';"
        ## Execute drop index query
        MazamaCoreUtils::logger.trace( '[mstools::filterSQMASSdb] Querying Database: %s', drop_idxs)
        DBI::dbExecute( db, drop_idxs )
        
        create_idxs <- "CREATE INDEX IF NOT EXISTS 'chrom_run_idx' ON 'CHROMATOGRAM' ('RUN_ID');
      CREATE INDEX IF NOT EXISTS 'data_chr_idx' ON 'DATA' ('CHROMATOGRAM_ID');
      CREATE INDEX IF NOT EXISTS 'data_sp_idx' ON 'DATA' ('SPECTRUM_ID');
      CREATE INDEX IF NOT EXISTS 'precursor_chr_idx' ON 'DATA' ('CHROMATOGRAM_ID');
      CREATE INDEX IF NOT EXISTS 'precursor_sp_idx' ON 'DATA' ('SPECTRUM_ID');
      CREATE INDEX IF NOT EXISTS 'product_chr_idx' ON 'DATA' ('CHROMATOGRAM_ID');
      CREATE INDEX IF NOT EXISTS 'product_sp_idx' ON 'DATA' ('SPECTRUM_ID');
      CREATE INDEX IF NOT EXISTS 'run_extra_idx' ON 'RUN_EXTRA' ('RUN_ID');
      CREATE INDEX IF NOT EXISTS 'spec_mslevel_idx' ON 'SPECTRUM' ('MSLEVEL');
      CREATE INDEX IF NOT EXISTS 'spec_rt_idx' ON 'SPECTRUM' ('RETENTION_TIME');
      CREATE INDEX IF NOT EXISTS 'spec_run_idx' ON 'SPECTRUM' ('RUN_ID');"
        ## Execute create index query
        MazamaCoreUtils::logger.trace( '[mstools::filterSQMASSdb] Querying Database: %s', create_idxs)
        DBI::dbExecute( db, create_idxs )
      }
      
      ##***********************************************
      ##    Clear unused space in db
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Vacuuming Database")
      DBI::dbExecute(db, "VACUUM")
      
      ##***********************************************
      ##    Disconnect fom DB 
      ##***********************************************
      MazamaCoreUtils::logger.trace( "[mstools::filterSQMASSdb] Disconnecting From Database: %s", sqmass_file)
      DBI::dbDisconnect( db )
    },
    error = function(e){
      MazamaCoreUtils::logger.error("[mstools::filterSQMASSdb] There was the following error that occured during function call...\n", e$message)
    }
  ) # End tryCatch
}# End function

