## Install python dependency
#' @import reticulate
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
find_python <- function(){
  ## Check if logging has been initialized
  if( MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  
  MazamaCoreUtils::logger.info(  "** Finding Python **"  )
  # "/home/justincsing/anaconda3/envs/mstools_test/"
  ## Check for available versions of python
  pythons_available <- reticulate::py_discover_config()
  ## Just use first isntance of python. @ TODO, might need to find a better way to do this
  use_first_python <- pythons_available$python_versions[1]
  MazamaCoreUtils::logger.warn( "Python Found: %s", use_first_python )
  Sys.setenv( RETICULATE_PYTHON=use_first_python )
  reticulate::use_python( use_first_python )
  
}

#' @import reticulate
#' @importFrom MazamaCoreUtils logger.isInitialized logger.info logger.error logger.warn 
install_python_dependencies <- function(method = "auto", conda = "auto") {
  
  ## Check if logging has been initialized
  if( MazamaCoreUtils::logger.isInitialized() ){
    log_setup()
  }
  MazamaCoreUtils::logger.info( "** Checking for Required Python Modules ** ")
  
  ## Check for available versions of python
  pythons_available <- reticulate::py_discover_config()
  ## Just use first isntance of python. @ TODO, might need to find a better way to do this
  use_first_python <- pythons_available$python_versions[1]
  
  # MazamaCoreUtils::logger.info("Installing modules for python: %s", Sys.getenv("RETICULATE_PYTHON"))
  
  pymsnumpress_available=FALSE
  cython_available=FALSE
  sqlite3_available=FALSE
  pandas_available=FALSE
  if ( !reticulate::py_module_available("PyMSNumpress") ){
    if ( !reticulate::py_module_available("cython") ){
      MazamaCoreUtils::logger.info( "** -> Installing Python Module: cython" )
      reticulate::py_install( "cython", pip=TRUE )
    }
    MazamaCoreUtils::logger.info( "** -> Installing Python Module: pymsnumpress" )
    reticulate::py_install( "PyMSNumpress", pip=TRUE )
    pymsnumpress_available=TRUE
  }
  if ( !reticulate::py_module_available("sqlite3") ){
    MazamaCoreUtils::logger.info( "** -> Installing Python Module: sqlite3" ) 
    reticulate::py_install("sqlite3", pip=TRUE )
    sqlite3_available=TRUE
  }
  if ( !reticulate::py_module_available("pandas") ){
    MazamaCoreUtils::logger.info( "** -> Installing Python Module: pandas" )
    reticulate::py_install("pandas", pip=TRUE )
    pandas_available=TRUE
  }
}

#' @import reticulate
python_dummy_modules <- function() {
  # Load the module and create dummy objects from it, all of which are NULL
  pymsnumpress <- reticulate::import( "PyMSNumpress", delay_load = TRUE, convert = FALSE )
  for (obj in c("decodeLinear", "decodeSlof") ) {
    assign(obj, NULL)
  }
  # Clean up
  rm(pymsnumpress)
  
  # Load the module and create dummy objects from it, all of which are NULL
  pyzlib <- reticulate::import( "zlib", delay_load = TRUE, convert = FALSE )
  for (obj in c("compress","decompress")) {
    assign(obj, NULL)
  }
  # Clean up
  rm(pyzlib)
  
  # Load the module and create dummy objects from it, all of which are NULL
  pybuiltins <- reticulate::import_builtins( convert = TRUE )
  for (obj in c("bytes","bytearray")) {
    assign(obj, NULL)
  }
  # Clean up
  rm(pybuiltins)
}

#' @import reticulate
.onload <- function( libname, pkgname ){
  pymsnumpress <- NULL
  pyzlib <- NULL
  pybuiltins <- NULL
  ## Use super assignment to update global reference to import
  #' @export
  pymsnumpress <<- reticulate::import( "PyMSNumpress", delay_load = TRUE, convert = FALSE )
  #' @export
  pyzlib <<- reticulate::import("zlib", delay_load = TRUE, convert = FALSE )
  #' @export
  pybuiltins <<- reticulate::import_builtins( convert = TRUE )
  
  # assignInMyNamespace(...) is meant for namespace manipulation
  # for (obj in c("decodeLinear", "decodeSlof") ) {
  #   assignInMyNamespace(obj, pymsnumpress[[obj]])
  # }
  # for (obj in c("compress","decompress")) {
  #   assignInMyNamespace(obj, pyzlib[[obj]])
  # }
  # for (obj in c("compress","decompress")) {
  #   assignInMyNamespace(obj, pybuiltins[[obj]])
  # }
  
}
