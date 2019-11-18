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
  return(pos_mod) 
  
}