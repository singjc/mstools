#' Find if an object exists in a list
#' Searches name structure for name search parameter
#' @param x A list object
#' @param name Name of sub-components in list
#' @return logical value if name is present in list
#' @export
isListObj <- function(x, name) {
  pos <- match(name, names(x))
  if (!is.na(pos)){ return(TRUE) } else { return(FALSE) }
  for (el in x) {
    if (class(el) == "list") {
      out <- Recall(el, name)
      if (!is.null(out)) return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

#' Find and Extract an obj from a structured list
#' @param x A list object
#' @param name Name of sub-components in list
#' @return returns the named object in list
#' @export
getListObj <- function(x, name) {
  pos <- match(name, names(x))
  if (!is.na(pos)) return(x[[pos]])
  for (el in x) {
    if (class(el) == "list") {
      out <- Recall(el, name)
      if (!is.null(out)) return(out)
    }
  }
}

#' check_sqlite_table
#' @param conn Connection to database
#' @param table Character vector to test for, if present in database
#' @param msg A character to pre-append to stop error message. (Optional)
#' @return Logical value, TRUE if table is present
#' 
#' @importFrom DBI dbExistsTable
#' @export
check_sqlite_table <- function( conn, table, msg="" ) {
  if( !DBI::dbExistsTable( conn, table ) ){
    out.msg <- sprintf("%s An Error occured! There was no %s Table found in %s.\nCheck to see if you are using the right file, or if the file is corrupted.\n", msg, table, conn@dbname)
    stop( out.msg, call.=FALSE )
  }  
}

#' Check is SCORE_IPF is in database
#' @param oswFile A character vector pointing to path of osw file
#' @return returns a logical value if SCORE_IPF is present or not
#' @export
Score_IPF_Present <- function( oswFile ){
  db <- DBI::dbConnect( RSQLite::SQLite(), oswFile )
  if ( DBI::dbExistsTable( db, "SCORE_IPF" ) ){
    use_ipf_score <- TRUE
  } else {
    use_ipf_score <- FALSE
  }
  DBI::dbDisconnect( db )
  return( use_ipf_score )
}

#' Convert Variable to character including NULL
#' @param x An object to coerce to character
#' @return returns a character vector
#' @export
as_character_null <- function( x ){
  
  if ( is.null(x) ){
    return( 'NULL' )
  } else {
    return( as.character( x ) )
  }
  
}

#' Convert list object to printable character vectory
#' @param list_obj A list object to coerce to character
#' @param collapse A character vector to collapse characters on
#' @return returns a character vector
#' @export
listTostring <- function( list_obj, collapse = '\n' ){
  paste( paste(names(list_obj), list_obj, sep=' = '), collapse = collapse )
}