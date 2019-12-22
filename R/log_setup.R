#install.packages("MazamaCoreUtils")
#// **********************************************************************************************
#//                         log_setup.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' 
#' @title Logging system setup
#' @description Setup logging system to write out error, trace, warnings and general information to file
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
log_setup <- function(){
# Setup Logging System
tmpDir <- getwd()
MazamaCoreUtils::logger.setup(
  traceLog = file.path( tmpDir,"mstools-trace.log" ),
  errorLog = file.path( tmpDir, "mstools-error.log" ),
  warnLog = file.path( tmpDir, "mstools-warn.log" ),
  infoLog = file.path( tmpDir, "mstools-info.log" )
)
}