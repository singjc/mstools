#' Get a list of data.tables for chromatrogram id mappings
#' @param input A shiny input variable that contains Working Directory Information
#' @param global A list variable containing paths to chromatogram files
#' @return (A list of mzRpwiz)
#' 
#' @importFrom tictoc tic toc
#' 
#' @export
getsqMassPntrs <- function( dataPath, runs, nameCutPattern = "(.*)(/)(.*)", chrom_ext=".chrom.sqMass"  ){
  
  sql_query <- sprintf("SELECT
CHROMATOGRAM.NATIVE_ID AS chromatogramId,
CHROMATOGRAM.ID AS chromatogramIndex
FROM CHROMATOGRAM
")
  
  chromFiles <- list.files(dataPath, pattern = ".sqMass", recursive = T, full.names = T)
  names(chromFiles) <- basename(chromFiles)
  
  ## Get filenames from osw files and check if names are consistent between osw and mzML files. ######
  filenames <- getRunNames( dataPath, oswMerged=TRUE, nameCutPattern = nameCutPattern, chrom_ext = chrom_ext )
  filenames <- filenames[filenames$runs %in% runs,]
  
  tictoc::tic('Pre-Loading mzML Chromatogram Files onto disk')
  mzPntrs <- list()
  for ( chromatogram_input_index_num in seq(1, length(filenames$runs)) ){
    tryCatch(
      expr = {
        tictoc::tic()
        run <- rownames(filenames)[ chromatogram_input_index_num ]
        current_filename <- filenames$runs[ chromatogram_input_index_num ]
        # message(sprintf("\rCacheing mzML for %s of %s runs", run, length(filenames$runs)))
        ## Get path for current chromatogram file
        chromatogram_file_i <-  chromFiles[ grepl(current_filename, names(chromFiles)) ][[1]]
        # Establish connection to sqlite database chromatogram file
        conn <- DBI::dbConnect( RSQLite::SQLite(), chromatogram_file_i )
        ## Get table of chromatogram incidces and respective transtion ids
        chromHead <- dplyr::collect( dplyr::tbl(conn, dbplyr::sql(sql_query)) )
        ## Store run id and mz object into master list
        mzPntrs[[run]] <- list()
        ## TODO At somepoint make this store the chrom data for sqmass maybe
        # mzPntrs[[run]]$mz <- mz 
        ## Store file path
        # mzPntrs[[run]]$mz <- chromatogram_file_i
        mzPntrs[[run]]$mz <- dplyr::collect( dplyr::tbl(conn, dbplyr::sql( "SELECT * FROM DATA"  )) )
        mzPntrs[[run]]$chromHead <- chromHead
        ## Append chromHead to sqMass DATA table
        mzPntrs[[run]]$mz <- merge( mzPntrs[[run]]$chromHead , mzPntrs[[run]]$mz, by.x="chromatogramIndex", by.y="CHROMATOGRAM_ID", all=T)
        colnames(mzPntrs[[run]]$mz)[which(grepl("chromatogramIndex", colnames(mzPntrs[[run]]$mz)))] <- "CHROMATOGRAM_ID"
        colnames(mzPntrs[[run]]$mz)[which(grepl("chromatogramId", colnames(mzPntrs[[run]]$mz)))] <- "FRAGMENT_ID"
        ## Append filename
        mzPntrs[[run]]$filename <- filenames$filename[ chromatogram_input_index_num ]
        ## Append Run Name
        mzPntrs[[run]]$run_name <- filenames$runs[ chromatogram_input_index_num ]
        ## Append Run Id
        mzPntrs[[run]]$run_id <- run
        
        ## Disconnect form database
        DBI::dbDisconnect(conn)
        ## End timer
        exec_time <- tictoc::toc(quiet = T)
        message(sprintf("\rCacheing sqMass for %s of %s runs: Elapsed Time = %s sec", run, length(filenames$runs), round(exec_time$toc - exec_time$tic, 3) ))
      },
      error = function(e){
        message(sprintf("[getsqMassChromIdMapping] There was an issue cacheing %s, skipping...: %s\n", current_filename, e$message))
      }
    ) # End tryCatch
  }
  tictoc::toc()
  
  return( mzPntrs )
}