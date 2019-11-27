#// **********************************************************************************************
#//                         getModificationPosition.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing


#' @export
#' @title Get Modificiation position from modified sequence. Currently only support phosphorylation.
#' @description This function can be used to get the modification position for a given modified sequence
#' Currently only support phosphorylation/UniMod:21. TODO: Make it useable for other cases
#' 
#' @param mod_seq A character vector of the modified position.
#' @param character_index A Logical value indicating to return the exact location of ammino acid modified, or the index of modified text in string.
#' @return A numeric value for modification position
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' 
getModificationPosition_ <- function( mod_seq, character_index=F ){ 
  if (character_index==F){
    modification_labels <- regmatches(mod_seq, gregexpr("\\(.*?\\)", mod_seq))[[1]]
    mod_seq_Phospho_only <- gsub( paste(gsub('\\)','\\\\)',gsub('\\(','\\\\(',modification_labels[!(grepl('\\(UniMod:21\\)|\\(Phospho\\)', modification_labels))])), collapse = '|'), '', mod_seq )
    pos_mod <- gregexpr('\\(UniMod:21\\)|\\(Phospho\\)', mod_seq_Phospho_only)
    pos_mod <- as.numeric(pos_mod[[1]])-1
  } else {
    pos_mod <- gregexpr('\\(UniMod:21\\)|\\(Phospho\\)', mod_seq)
    pos_mod <- as.numeric(pos_mod[[1]])-1
  }
  
  ## If there is no modification, pos_mod is negative, then return 0 for pos_mod
  if ( pos_mod < 0 ){
    pos_mod <- 0
  }
  return(pos_mod) 
  
}