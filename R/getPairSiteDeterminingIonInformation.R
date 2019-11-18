getPairSiteDeterminingIonInformation_ <- function( uni_mod, len_unmod ) {
  ###############################
  ## Get Site Determining Ions ##
  ###############################
  # Get modification locationgs on sequence
  mod_positions <- sapply(uni_mod, function(seq){ getModificationPosition_(seq)})
  all_mod_positions <- mod_positions
  mod_positions <- as.numeric(matrix(c(setdiff(all_mod_positions[[1]],all_mod_positions[[2]]), setdiff(all_mod_positions[[2]],all_mod_positions[[1]])), nrow=1, ncol=2))
  if ( all(is.na(mod_positions)) ){ cat(red(paste(uni_mod, collapse=', '), ' do not differ in the modification position, these are most likely not isoforms of eachother!!\n', sep='')); return( 'skip' )}
  if (is.matrix(all_mod_positions)){
    names(mod_positions) <- c(colnames(all_mod_positions)[1], colnames(all_mod_positions)[2])
  } else {
    names(mod_positions) <- c(names(all_mod_positions)[1], names(all_mod_positions)[2])  
  }
  ### @ TODO: Need to add a check for uni_mod that are not isoforms of each other.
  
  uni_mod_list <- list()
  
  first_mod <- mod_positions[ which( !(mod_positions==max(mod_positions)) ) ]
  second_mod <- mod_positions[ which( (mod_positions==max(mod_positions)) ) ]
  
  y_site_determining_end = len_unmod - first_mod
  y_site_determining_start =  len_unmod - second_mod + 1 # + 1 to not include a.a. mod is on
  
  b_site_determining_end = len_unmod - y_site_determining_start
  b_site_determining_start = len_unmod - y_site_determining_end
  
  uni_mod_list[[names(first_mod)]][['Modification']] <- names(first_mod)
  uni_mod_list[[names(first_mod)]][['Position']] <- as.numeric(first_mod)
  uni_mod_list[[names(first_mod)]][['Type']] <- 'y'
  uni_mod_list[[names(first_mod)]][['site_determining_start']] <- y_site_determining_start
  uni_mod_list[[names(first_mod)]][['site_determining_end']] <- y_site_determining_end
  
  uni_mod_list[[names(second_mod)]][['Modification']] <- names(second_mod)
  uni_mod_list[[names(second_mod)]][['Position']] <- as.numeric(second_mod)
  uni_mod_list[[names(second_mod)]][['Type']] <- 'b'
  uni_mod_list[[names(second_mod)]][['site_determining_start']] <- b_site_determining_start
  uni_mod_list[[names(second_mod)]][['site_determining_end']] <- b_site_determining_end
  
  
  # uni_mod_list[[names(first_mod)]][['Modification']] <- names(first_mod)
  # uni_mod_list[[names(first_mod)]][['Position']] <- as.numeric(first_mod)
  # uni_mod_list[[names(first_mod)]][['y_site_determining_start']] <- y_site_determining_start
  # uni_mod_list[[names(first_mod)]][['y_site_determining_end']] <- y_site_determining_end
  # uni_mod_list[[names(second_mod)]][['b_site_determining_start']] <- b_site_determining_start
  # uni_mod_list[[names(second_mod)]][['b_site_determining_end']] <- b_site_determining_end
  
  
  
  return( uni_mod_list )
}