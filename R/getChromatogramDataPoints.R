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
#' @return A list of fragment ids containing 2 arrays for Retention time and Intensity
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
getChromatogramDataPoints_ <- function( filename, frag_ids ){
  ## Get File Extension Type
  fileType <- gsub( '.*\\.', '', filename)
  ## Extract Chromatogram Data
  if ( tolower(fileType)=='sqmass' ){
    # Read in an sqMass Chromatogram ------------------------------------------
    
    if ( F ){
      cat('Reading in chromatogram of ', crayon::blue$bold$underline('sqmass type.\n', sep=''))
      
      
      filename <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/M_pools_re-run//lgillet_L160920_012-Manchester_dirty_phospho_-_Pool_M6_-_SW/lgillet_L160920_012-Manchester_dirty_phospho_-_Pool_M6_-_SW.mzML.gz_osw_chrom.sqMass"
      frag_ids <- list()
      frag_ids[[1]] <- c("7788", "7789", "7790", "7791", "7792", "7793")
      
      
      # chrom <- getChromatogramsbyIndice_( filename, frag_ids )
      
      # Connect to database
      sqmass_db <- DBI::dbConnect( RSQLite::SQLite(), filename )
      
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
      chrom_index_df <- dplyr::collect( dplyr::tbl(sqmass_db, dplyr::sql(sql_query)) )
      
      # @TODO: If chrom_index_df is empty, throw an error
      
      
      "
        Get data from multiple chromatograms chromatogram
        - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        - data_type is one of 0 = mz, 1 = int, 2 = rt
        - data contains the raw (blob) data for a single data array
    "
      
      stmt <- "SELECT CHROMATOGRAM_ID, COMPRESSION, DATA_TYPE, DATA FROM DATA WHERE CHROMATOGRAM_ID IN ("
      for ( myid in as.matrix(chrom_index_df[,1]) ){
        stmt <- paste( stmt,  myid, ",", sep='' )
      }
      stmt <- substr(stmt,1,nchar(stmt)-1)
      stmt <- paste(stmt, ')', sep='')
      
      data <- dplyr::collect( dplyr::tbl(sqmass_db, dplyr::sql(stmt)) )
      
      rt_array <- list()
      intensity_array <- list()
      
      for ( row in seq(1, nrow(data)) ){
        row = 1
        result <- numeric()
        
        data_row <- data[ row, ]
        
        if ( data_row$COMPRESSION==5 ){
          # install.packages("remotes")
          # remotes::install_github("statwonk/Rcompression")
          
          tmp <- base::charToRaw( Rcompression::uncompress( data_row$DATA[[1]] ) )
          
          
          
          if ( length( tmp ) > 0 ){
            mstools::decodeLinear( tmp, result)
          }
          
        }
        
      }
      
      
      
      # Disconnect from database
      DBI::dbDisconnect(sqmass_db)
      
      
      return(chrom)
    } else {
      cat( crayon::red$bold$underline(fileType, ' FileType is not supported yet, coming soon!!\n'), sep='')
    }
    
    
  } else if ( tolower(fileType)=='mzml' ){
    # Read in an mzML chromatogram --------------------------------------------
    cat('Reading in chromatogram of ', crayon::blue$bold$underline('mzML type.\n', sep=''))
    # Create an mzR object that stores all header information, and use ProteoWizard api to access data from MzML file
    mz_object <- mzR::openMSfile(filename, backend = "pwiz", verbose = T)
    # Get header information for chromtagograms
    chromHead <- mzR::chromatogramHeader(mz_object)
    # Extract all the indices of chromatograms that match the transition names of the ones found in the TargetPetides file
    chromatogramIndices <- chromHead$chromatogramIndex[ match(frag_ids[[1]], chromHead$chromatogramId)  ]
    # Check how many chromatogramIndices are present to extract
    if ( length(chromatogramIndices)==1 ){
      rawChrom <- list(mzR::chromatograms(mz_object, chromatogramIndices))
    } else if ( length(chromatogramIndices)>1 ) {
      rawChrom <- mzR::chromatograms(mz_object, chromatogramIndices)
    } else {
      cat( crayon::red$bold$underline('There was no Chromatogramphic data for the following fragment(s): ', base::paste(frag_ids[[1]], collapse = ', ')), sep='')
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
    cat( crayon::red$bold$underline(fileType, ' FileType is not supported!!\n'), sep='')
  }
} ## End Function