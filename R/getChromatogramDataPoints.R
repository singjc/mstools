#// **********************************************************************************************
#//                         getChromatogramDataPoints.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' 
#' @export
#' @title Extract Chromatographic data (Retention Time and Intensity) from mzML or sqMass chromatograms.
#' @description This function can be used to extract chromatogram retention and intensity values for a given
#' list of fragment ids.
#' 
#' @param filename A character vector of the absolute path and filename of the chromatogram file. (Must be .mzML or sqMass format)
#' @param frag_ids A list of a vector containing fragment ids
#' @param id_type A character indicating if the passed frag_ids is of type 'transition_ids' or 'chromatogramIndex'
#' @param name_time A character vector to name retention column
#' @param name_intensity A character vector to name Intensity column. You can pass a character vectory of a function that will be evaluated to a string. i.e. "paste0('X',  data_row$FRAGMENT_ID )"
#' @param mzPntrs A list object containing cached mzR objects.
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @return A list of fragment ids containing 2 arrays for Retention time and Intensity
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#'
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite 
#' @importFrom dplyr collect tbl %>% filter
#' @importFrom dbplyr sql 
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
#' @importFrom tools file_ext
#' @importFrom crayon blue bold underline red 
#' @import reticulate
#' @importFrom mzR openMSfile chromatogramHeader chromatograms 
getChromatogramDataPoints_ <- function( filename, frag_ids, id_type='transition_ids', name_time = 'RT', name_intensity = "Int", mzPntrs=NULL, SgolayFiltOrd=NULL, SgolayFiltLen=NULL  ){
  ## Check if logging has been initialized
  if( !MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  ## Get File Extension Type
  fileType <- (tools::file_ext(filename))
  ## Extract Chromatogram Data
  if ( tolower(fileType)=='sqmass' ){
    # Read in an sqMass Chromatogram ------------------------------------------
    
    if ( F ){
      filename <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/M_pools_re-run//lgillet_L160920_012-Manchester_dirty_phospho_-_Pool_M6_-_SW/lgillet_L160920_012-Manchester_dirty_phospho_-_Pool_M6_-_SW.mzML.gz_osw_chrom.sqMass"
      frag_ids <- list()
      frag_ids[[1]] <- c("7788", "7789", "7790", "7791", "7792", "7793")
      data <- mzPntrs$run11$mz
      id_type = 'chromatogramIndex'
      frag_ids= chromIndices
      name_time = 'RT'
      name_intensity = "Int"
    }
    
    MazamaCoreUtils::logger.trace( sprintf( '[mstools::getChromatogramDataPoints_]\tReading in chromatogram of sqmass type.\n' ) )
    
    ##*********************************************
    ##      Python Setup
    ##*********************************************
    
    if ( !reticulate::py_available()  ){
      if ( Sys.getenv("RETICULATE_PYTHON")=="" & reticulate::py_available()==FALSE ){
        find_python()
      }
      install_python_dependencies()
      MazamaCoreUtils::logger.trace( "[mstools::getChromatogramDataPoints_] ** Loading Python Modules **")
      .onload()
    }
    
    ##********************************************
    ##    Do Actual Int and RT Extraction 
    ##********************************************
    if ( is.null(mzPntrs) ) {
      # Connect to database
      MazamaCoreUtils::logger.trace( sprintf("[mstools::getChromatogramDataPoints_] Connecting to Database: %s", filename) )
      sqmass_db <- DBI::dbConnect( RSQLite::SQLite(), filename )
    }
    
    ##*********************************************************************
    ##    Get Chromatogram ID to Transition Fragment IDs mapping
    ##*********************************************************************
    if ( is.null(mzPntrs) & id_type=='transition_ids' ){
      
      # Query statement
      sql_query <- "SELECT * FROM CHROMATOGRAM WHERE NATIVE_ID in ("
      for (precursor in frag_ids){
        for (current_id in precursor){
          sql_query <- paste( sql_query, "'", current_id, "', ", sep='')
        }
      }
      sql_query <- substr(sql_query,1,nchar(sql_query)-2)
      sql_query <- paste(sql_query, ')', sep='')
      
      # Query Databasse
      MazamaCoreUtils::logger.trace( sprintf("[mstools::getChromatogramDataPoints_] Querying Database: %s", sql_query) )
      chrom_index_df <- dplyr::collect( dplyr::tbl(sqmass_db, dbplyr::sql(sql_query)) )
      colnames(chrom_index_df) <- c('CHROMATOGRAM_ID', 'RUN_ID', 'FRAGMENT_ID')
      
      # @TODO: If chrom_index_df is empty, throw an error
      
      "
        Get data from multiple chromatograms chromatogram
        - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        - data_type is one of 0 = mz, 1 = int, 2 = rt
        - data contains the raw (blob) data for a single data array
    "
      MazamaCoreUtils::logger.trace( "[mstools::getChromatogramDataPoints_] ** sqMass: Extracting Chromatogram data from database **")
      stmt <- "SELECT CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA FROM DATA WHERE CHROMATOGRAM_ID IN ("
      for ( myid in as.matrix(chrom_index_df[,1]) ){
        stmt <- paste( stmt,  myid, ",", sep='' )
      }
      stmt <- substr(stmt,1,nchar(stmt)-1)
      stmt <- paste(stmt, ')', sep='')
      MazamaCoreUtils::logger.trace( sprintf("[mstools::getChromatogramDataPoints_] Querying Database: %s", stmt) )
      data <- dplyr::collect( dplyr::tbl(sqmass_db, dbplyr::sql(stmt)) )
      
      data <- merge(data, chrom_index_df, by="CHROMATOGRAM_ID")
      
    } else if ( !is.null(mzPntrs)  ){ ## Had & id_type=='chromatogramIndex here before, but I think this is obsolete now.
      MazamaCoreUtils::logger.trace( "[mstools::getChromatogramDataPoints_] ** sqMass: Extracting Chromatogram data from cached pointers data.table **")
      if ( class(mzPntrs)=="list" ){
        if ( isListObj(mzPntrs, "mz") ){
          data <- mzPntrs$mz
        } else {
          MazamaCoreUtils::logger.error( sprintf("[mstools::getChromatogramDataPoints_] The passed cached mzPntrs object did not contain the raw chromatogram data in the list object.\nThe contents of the passed object contains %s.\n", names(mzPntrs)) )
        }
      } else {
        if ( dim(mzPntrs)[1] > 0 ){
        data <- mzPntrs
        } else {
          MazamaCoreUtils::logger.error( sprintf("[mstools::getChromatogramDataPoints_] The passed cached mzPntrs did not contain any data. Dimension of rows passed data.frame %s.\n", dim(mzPntrs)[1]) )
        }
      }
      if ( id_type=='chromatogramIndex' ){
        data %>%
          dplyr::filter( CHROMATOGRAM_ID %in% unlist(frag_ids) ) -> data
      } else {
        data %>%
          dplyr::filter( FRAGMENT_ID %in% unlist(frag_ids) ) -> data
      }
      
    }
    
    chrom <- list()
    for ( row in seq(1, nrow(data)) ){
      # row = 2
      ## Initialize an empty python list to store results in
      result <- reticulate::r_to_py( list() )
      ## Subset data for current row of data
      data_row <- data[ row, ]
      ## Check which method current row data was compressed with
      if ( data_row$COMPRESSION==5 ){
        pymsnumpress$decodeLinear(data = pybuiltins$bytearray( pyzlib$decompress( pybuiltins$bytes( reticulate::r_to_py( data_row$DATA[[1]] )))), 
                                  result = result )
        result <- reticulate::py_to_r( result )
      }
      if ( data_row$COMPRESSION==6 ){
        pymsnumpress$decodeSlof(data = pybuiltins$bytearray( pyzlib$decompress( pybuiltins$bytes( reticulate::r_to_py( data_row$DATA[[1]] )))), 
                                result = result )
        result <- reticulate::py_to_r( result )
      }
      
      if ( length(result) == 0 ){
        result = numeric()
      }
      if ( data_row$DATA_TYPE == 1 ){
        # Should Intensity be smoothed
        if ( !is.null(SgolayFiltOrd) & !is.null(SgolayFiltLen) ) result <- signal::sgolayfilt(result, p = SgolayFiltOrd, n = SgolayFiltLen)
        chrom[[ data_row$FRAGMENT_ID ]][[ tryCatch( expr={eval(parse(text=name_intensity))}, error=function(e) {name_intensity} ) ]] = result
      } else if ( data_row$DATA_TYPE == 2){
        chrom[[ data_row$FRAGMENT_ID ]][[name_time]] = result
      } else {
        MazamaCoreUtils::logger.error( sprintf("[mstools::getChromatogramDataPoints_] Only expected RT or Intensity data for chromatogram. Expected DATA_TYPE to be 1 or 2, instead got %s", data_row$DATA_TYPE) )
      }
    }
    if ( is.null(mzPntrs) ) {
      # Disconnect from database
      MazamaCoreUtils::logger.trace( sprintf("[mstools::getChromatogramDataPoints_] Disconnecting From Database: %s", filename) )
      DBI::dbDisconnect(sqmass_db)
    }
    return(chrom)
    
  } else if ( tolower(fileType)=='mzml' ){
    if ( F ){
      filename <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/mzml/chludwig_K150309_013_SW_0.chrom.mzML"
      filename <- '/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/DrawAlignR/inst/extdata/mzml/chludwig_K150309_013_SW_0.chrom.mzML'
      filename <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/mzML_Chroms_Decomp/chludwig_K150309_013_SW_0_osw_chrom.mzML"
    }
    # Read in an mzML chromatogram --------------------------------------------
    MazamaCoreUtils::logger.info(paste('[mstools::getChromatogramDataPoints_] Reading in chromatogram of ', crayon::blue$bold$underline('mzML type.\n'), sep=''))
    if ( is.null(mzPntrs) ){
      MazamaCoreUtils::logger.info( "[mstools::getChromatogramDataPoints_] ** mzR: Loading mzML chromatogram into mz_object **")
      # Create an mzR object that stores all header information, and use ProteoWizard api to access data from MzML file
      mz_object <- mzR::openMSfile(filename, backend = "pwiz", verbose = T)
      # Get header information for chromtagograms
      chromHead <- mzR::chromatogramHeader(mz_object)
    } else {
      mz_object <- mzPntrs$mz
      chromHead <- mzPntrs$chromHead
    }
    MazamaCoreUtils::logger.info( "[mstools::getChromatogramDataPoints_] ** mzR: Extracting chromatogram indices **")
    # Extract all the indices of chromatograms that match the transition names of the ones found in the TargetPetides file
    chromatogramIndices <- chromHead$chromatogramIndex[ match(frag_ids[[1]], chromHead$chromatogramId)  ]
    MazamaCoreUtils::logger.info( "[mstools::getChromatogramDataPoints_] ** mzR: Extracting chromatographic data **")
    # Check how many chromatogramIndices are present to extract
    if ( length(chromatogramIndices)==1 & !is.na(chromatogramIndices) ){
      rawChrom <- list(mzR::chromatograms(mz_object, chromatogramIndices))
    } else if ( length(chromatogramIndices)>1 ) {
      rawChrom <- mzR::chromatograms(mz_object, chromatogramIndices)
    } else {
      MazamaCoreUtils::logger.error( paste(crayon::red$bold$underline('[mstools::getChromatogramDataPoints_] There was no Chromatogramphic data for the following fragment(s): ', base::paste(frag_ids[[1]], collapse = ', ')), sep='') )
      chrom <- list(); rawChrom_idx <- 1
      for ( fragment_id in frag_ids[[1]] ){
        chrom[[ fragment_id ]] <- list(RT=NaN,
                                       Int=NaN)
        rawChrom_idx <- rawChrom_idx + 1
      }
      return( chrom )
    }
    chrom <- list(); rawChrom_idx <- 1
    for ( fragment_id in frag_ids[[1]] ){
      chrom[[ fragment_id ]] <- list(RT=rawChrom[[rawChrom_idx]][,1],
                                     Int=rawChrom[[rawChrom_idx]][,2])
      rawChrom_idx <- rawChrom_idx + 1
    }
    rm(mz_object)
    return(chrom)
  } else {
    MazamaCoreUtils::logger.error( paste(crayon::red$bold$underline(fileType, '[mstools::getChromatogramDataPoints_] FileType is not supported!!\n'), sep='') )
  }
} ## End Function
