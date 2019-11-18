getSiteDeterminingIonInformation_ <- function( uni_mod, len_unmod ) {
  ###############################
  ## Get Site Determining Ions ##
  ###############################
  # Get modification locationgs on sequence
  mod_positions <- sapply(uni_mod, function(seq){ getModificationPosition_(seq)})
  all_mod_positions <- mod_positions
  if (is.matrix(all_mod_positions)){
    names(mod_positions) <- c(colnames(all_mod_positions)[1], colnames(all_mod_positions)[2])
  } else {
    names(mod_positions) <- c(names(all_mod_positions)[1])  
  }
  ### @ TODO: Need to add a check for uni_mod that are not isoforms of each other.
  
  uni_mod_list <- list()
  
  first_mod <- mod_positions[ which( (mod_positions==max(mod_positions)) ) ]
  
  y_site_determining_end = len_unmod - first_mod
  
  b_site_determining_start = len_unmod - y_site_determining_end
  
  uni_mod_list[['y']][['Modification']] <- names(first_mod)
  uni_mod_list[['y']][['Position']] <- as.numeric(first_mod)
  uni_mod_list[['y']][['site_determining_start']] <- len_unmod
  uni_mod_list[['y']][['site_determining_end']] <- y_site_determining_end
  
  uni_mod_list[['b']][['Modification']] <- names(first_mod)
  uni_mod_list[['b']][['Position']] <- as.numeric(first_mod)
  uni_mod_list[['b']][['site_determining_start']] <- b_site_determining_start
  uni_mod_list[['b']][['site_determining_end']] <- 1
  
  
  # uni_mod_list[[names(first_mod)]][['Modification']] <- names(first_mod)
  # uni_mod_list[[names(first_mod)]][['Position']] <- as.numeric(first_mod)
  # uni_mod_list[[names(first_mod)]][['y_site_determining_start']] <- y_site_determining_start
  # uni_mod_list[[names(first_mod)]][['y_site_determining_end']] <- y_site_determining_end
  # uni_mod_list[[names(second_mod)]][['b_site_determining_start']] <- b_site_determining_start
  # uni_mod_list[[names(second_mod)]][['b_site_determining_end']] <- b_site_determining_end
  

  
  return( uni_mod_list )
}