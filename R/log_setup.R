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
#' @import MazamaCoreUtils
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @import MazamaCoreUtils
log_setup <- function(){
# Setup Logging System
tmpDir <- getwd()
cat('log directory: ', tmpDir, '\n', sep = ' ')
time_stamp <- format(Sys.time(), "%a_%b_%d_%Hh%Mm%Ss_%Y")
MazamaCoreUtils::logger.setup(
  traceLog = file.path( tmpDir, sprintf("mstools-trace_%s.log", time_stamp) ),
  errorLog = file.path( tmpDir, sprintf("mstools-error_%s.log", time_stamp) ),
  warnLog = file.path( tmpDir, sprintf("mstools-warn_%s.log", time_stamp) ),
  infoLog = file.path( tmpDir, sprintf("mstools-info_%s.log", time_stamp) )
)
MazamaCoreUtils::logger.setLevel( level='INFO' )
}