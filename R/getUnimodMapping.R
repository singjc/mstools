#// **********************************************************************************************
#//                         getUnimodMapping.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Generate a UniMod mapping table
#' @description This function can be used to generate a data.table containing
#' UniMod ID to Modification code name mapping
#' 
#' @param unimod_xml A character path indicating the path to a unimod.xml path
#' @return A data.table containing UniMod ID to Modification Code Name mapping
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
#' @import dplyr
#' @importFrom dplyr %>%
#' @import xml2
#' @import purrr
generateUniModMapping_ <- function( unimod_xml=NULL ){
  # library(xml2)
  # library(purrr)
  if ( is.null(unimod_xml) ){
    ## Path to unimod xml database
    unimod_xml <- paste0(getwd(), "/R/sysdata/unimod_tables.xml")
  }
  ## Read in unimod xml
  unimod_dt <- xml2::read_xml( unimod_xml )
  ## Establish list of namespace in file
  ns <- xml2::xml_ns(unimod_dt)
  ## Extract modifications row information from unimod xml
  tmp <- xml2::xml_find_all(unimod_dt, "//d1:modifications_row", ns)
  ## Extract record id and modifications code name and map to datatable
  purrr::map(tmp, xml2::xml_attrs) %>%
    purrr::map_df( as.list ) %>%
    dplyr::select( record_id, ex_code_name ) -> unimod_mapping
  # save(unimod_mapping, file = './data/unimod_mapping.RData')
  ## Return a datable with UniMod ID and modification code name mapping
  return( unimod_mapping )
}

## Dummy Documentation
#' @title Conversion from UniMod ID to modification name
#' @usage unimodTocodename( mod_seq, out='sequence' )
#' @param mod_seq A character vector containing of a peptide sequence, containing UniMod IDs to be converted to code names
#' @name unimodTocodename
#' @examples 
#' mod_seq <- "EGHAQNPM(UniMod:35)EPS(UniMod:21)VPQLS(UniMod:21)LM(UniMod:35)DVK"
#' ## Will return peptide sequence with modifciation code name
#' unimodTocodename( mod_seq )
#' ## Out
#' ## [1] "EGHAQNPM(Oxidation)EPS(Phospho)VPQLS(Phospho)LM(Oxidation)DVK"
#' @import plyr
#' @importFrom dplyr %>%
#' @import stringr
NULL

#' @export
unimodTocodename <- function( mod_seq, out='sequence'){
  DEBUG = F
  if ( DEBUG ){
    mod_seq <- "EGHAQNPMEPSVPQLSLMDVK"
    mod_seq <- "EGHAQNPMEPSVPQLS(UniMod:21)LMDVK"
    mod_seq <- "EGHAQNPM(UniMod:35)EPS(UniMod:21)VPQLS(UniMod:21)LM(UniMod:35)DVK"
  }
  ## Get Modifications
  modification_labels <- base::regmatches(mod_seq, gregexpr("\\(.*?\\)", mod_seq))[[1]]
  ## Check if there were any modifications
  if ( length(modification_labels) >= 1 ){
  ## get UniMod IDs
  query_unimod_ids <- gsub( "\\(UniMod:(\\d+)\\)", "\\1", (modification_labels) )
  ## convert query UniMod IDs to modification code name
  out_code_names <- plyr::mapvalues(query_unimod_ids, from = mstools::unimod_mapping$record_id, to = mstools::unimod_mapping$ex_code_name, warn_missing = FALSE)
  names(out_code_names) <- gsub('\\(|\\)','',modification_labels)
  ## Convery to code name
  mod_seq %>%
    stringr::str_replace_all( out_code_names ) -> converted_mod_seq
  } else {
    out_code_names <- NULL
    converted_mod_seq <- mod_seq
  }
  if ( out=='sequence' ){
    return( converted_mod_seq )
  } else {
    return( out_code_names )
  }
}

## Dummy Documentation
#' @title Conversion from modification codename to UniMod ID
#' @usage codenameTounimod( mod_seq, out='sequence' )
#' @param mod_seq A character vector containing of a peptide sequence, containing modification code names to be converted to UniMod IDs
#' @name codenameTounimod 
#' @examples 
#' mod_seq <- "EGHAQNPM(Oxidation)EPS(Phospho)VPQLS(Phospho)LM(Oxidation)DVK"
#' ## Will return peptide sequence with UniMod ID
#' codenameTounimod ( mod_seq )
#' ## Out
#' ## [1] "EGHAQNPM(UniMod:35)EPS(UniMod:21)VPQLS(UniMod:21)LM(UniMod:35)DVK"
#' @import plyr
#' @importFrom dplyr %>%
#' @import stringr
NULL

#' @export
codenameTounimod <- function( mod_seq, out='sequence'){
  DEBUG = F
  if ( DEBUG ){
    mod_seq <- "EGHAQNPMEPSVPQLSLMDVK"
    mod_seq <- "EGHAQNPMEPSVPQLS(Phospho)LMDVK"
    mod_seq <- "EGHAQNPM(Oxidation)EPS(Phospho)VPQLS(Phospho)LM(Oxidation)DVK"
    mod_seq <- "ANS(Phospho)SPTTNIDHLK(Label:13C(6)15N(2))"
  }
  ## Get Modifications
  ## Need to make this more robust, somehow
  modification_labels <- base::regmatches(mod_seq, gregexpr("\\(\\w+\\)|\\(Label:13C\\(\\d+\\)15N\\(\\d+\\)\\)", mod_seq))[[1]]
  ## Check if there were any modifications
  if ( length(modification_labels) >= 1 ){
    ## convert query UniMod IDs to modification code name
    out_code_names <- plyr::mapvalues( gsub('^\\(|\\)$','',modification_labels), from = mstools::unimod_mapping$ex_code_name, to = paste0("UniMod:", mstools::unimod_mapping$record_id), warn_missing = FALSE)
    names(out_code_names) <- gsub('^\\(|\\)$','',modification_labels)
    ## Convery to code name
    mod_seq %>%
      stringr::str_replace_all( out_code_names ) -> converted_mod_seq
    converted_mod_seq <- gsub('\\Label:13C\\(\\d+\\)15N\\(\\d+\\)', out_code_names[grepl('Label.*', names(out_code_names))], converted_mod_seq )
  } else {
    out_code_names <- NULL
    converted_mod_seq <- mod_seq
  }
  if ( out=='sequence' ){
    return( converted_mod_seq )
  } else {
    return( out_code_names )
  }
}