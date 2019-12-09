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
  DEBUG = F
  if ( DEBUG ){
    
    mod_seq <- "EGHAQNPMEPSVPQLS(UniMod:21)LMDVK"
    mod_seq <- "EGHAQNPMEPS(UniMod:21)VPQLS(UniMod:21)LM(UniMod:35)DVK"
  }
  
  if (character_index==F){
    modification_labels <- regmatches(mod_seq, gregexpr("\\(.*?\\)", mod_seq))[[1]]
    mod_seq_Phospho_only <- gsub( paste(gsub('\\)','\\\\)',gsub('\\(','\\\\(',modification_labels[!(grepl('\\(UniMod:21\\)|\\(Phospho\\)', modification_labels))])), collapse = '|'), '', mod_seq )
    
    modification_index_list <- list()
    for ( mod in unique(modification_labels) ) {
      
      ## Remove other modifications from the sequence that are not being accessed to get a more accurate position index
      mods_to_remove_from_sequence <- modification_labels[ !(modification_labels %in% mod) ]
      if ( length(mods_to_remove_from_sequence)>0 ){
        ## Remove the other modifications not being assessed.
        current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mods_to_remove_from_sequence )), "", mod_seq )
      } else {
        current_mod_sequence <- mod_seq
      }
      ## Replace current Modification string with an identifier to get amino acid index
      current_mod_sequence <- gsub( gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mod )), "@", current_mod_sequence )
      ## Get positions of current modification based on identifier
      pos_mod <- gregexpr( '@', current_mod_sequence )
      ## Get actual amino acid index
      pos_mod <- as.numeric(pos_mod[[1]]) - 1
      ## Store position in list
      modification_index_list[[mod]] <- pos_mod
    }
    return( modification_index_list )    
    
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