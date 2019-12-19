## Install python dependency

find_python <- function(){
  
  # "/home/justincsing/anaconda3/envs/mstools_test/"
  ## Check for available versions of python
  pythons_available <- reticulate::py_discover_config()
  ## Just use first isntance of python. @ TODO, might need to find a better way to do this
  use_python <- pythons_available$python_versions[1]
  reticulate::use_python( use_python )
  
}

install_python_dependencies <- function(method = "auto", conda = "auto") {
  pymsnumpress_available=FALSE
  sqlite3_available=FALSE
  pandas_available=FALSE
  if ( !reticulate::py_module_available("PyMSNumpress") ){
    system( paste( use_python, "-m pip install pymsnumpress", sep=" "))
    pymsnumpress_available=TRUE
  }
  if ( !reticulate::py_module_available("sqlite3") ){
    reticulate::py_install("sqlite3", method = method, conda = conda)
    sqlite3_available=TRUE
  }
  if ( !reticulate::py_module_available("pandas") ){
    reticulate::py_install("pandas", method = method, conda = conda)
    pandas_available=TRUE
  }
}

python_dummy_modules <- function() {
  # Load the module and create dummy objects from it, all of which are NULL
  pymsnumpress <- reticulate::import( "PyMSNumpress", delay_load = FALSE, convert = FALSE )
  for (obj in c("decodeLinear", "decodeSlof") ) {
    assign(obj, NULL)
  }
  # Clean up
  rm(pymsnumpress)
  
  # Load the module and create dummy objects from it, all of which are NULL
  pyzlib <- reticulate::import( "zlib", delay_load = FALSE, convert = FALSE )
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
