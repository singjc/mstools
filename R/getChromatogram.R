plotChromatogram <- function(data, run, precursor_id, return_out='plot', filter_transitions.include=NULL, filter_transitions.exclude=NULL, filter_RT=NULL, verbose=TRUE){
  'Definition:
  This function is used to plot an XIC (extracted ion chromatogram) given a list containing MS runs, that contains osw peak scoring information and
  retention time and intensity values for each transition.
  
  Inputs:
  data - should be a nested list of MS runs, see below:
  data 
  ----> Run1 
  ----------> osw.out # This is a data.table that is read in from openswathworkflows-pyprophet.tsv output file. This contains scored peak rank groups
  ----------> osw.out.filtered # This is a data.table that is a filtered version of osw.out. Filtered based on Peak-rank group=1, Decoy=0 and m-score<0.01
  ----------> frag_id_list # This is a list of precursors with 6 transition group entries. See below..
  ------------------------> [[1]]  "0" "1" "2" "3" "4" "5" # Precursor id whould be a number index matching in list
  ------------------------> [[2]]  "10" "11" "6"  "7"  "8"  "9"
  ----------> chroms_list # This is a list of precursors containing vectors of RT and Intensity
  ------------------------> 1 # This is the precursor id
  ---------------------------> [[1]] 3865.4 3868.9 3872.5 3876.0 3879.5 ... # This is the RT array
  ---------------------------> [[2]] 48.00060 24.00051  0.00000 48.00060 24.00051 ... # This is the Intensity array
  run - should be a string value indicating the run index to extract in the data list. ^i.e. "Run1"
  precursor_id - should be a numeric value indicating the precursor id to extract. ^i.e. 1
  
  Output:
  g - is a ggplot object
  '
  ## Logging function for verbose output
  log_ <- function( msg, verbose ){
    if( isTRUE(verbose) ){
      cat( msg )
    }
  }
  if (is.character(data)){
    data <- readRDS( data )
  }
  ## @TODO, need to fix boolean indexing to make it more robut for different cases.
  # peptide_subset <- data[[run]]$osw.out.filtered[grepl(paste('^',precursor_id,'$',sep=''),data[[run]]$osw.out.filtered$transition_group_id),]
  peptide_subset <- data[[run]]$osw.out.filtered[ data[[run]]$osw.out.filtered$transition_group_id %in% precursor_id,]
  # get transition ids from openswath output file
  peptide_tranisition_ids <- (mclapply(strsplit(peptide_subset$aggr_Fragment_Annotation, ';'), function(x){gsub('_[y,b]\\d+_\\d+$','',x)})[[1]])
  
  
  
  
  chrom_list_indice <- grep( paste(paste('^', peptide_tranisition_ids, '$',sep=''), collapse = '|'), names(data[[run]][['chroms_list']]) )
  
  
  
  data_points <- data[[run]][['chroms_list']][chrom_list_indice]
  product_mz_val <- data[[run]][['productmz_list']][chrom_list_indice]
  for (i in seq(1:length(data_points))){
    names(data_points[[i]]) <- c('RT','Int')
    data_points[[i]]$Int <- sgolayfilt(data_points[[i]]$Int, p = 4, n = 9)
    data_points[[i]]$ProdMZ <- c((rep(product_mz_val[[i]], length(data_points[[i]]$Int))))
  }
  # data_points$ProdMZ <- product_mz_val
  
  df <- bind_rows(mclapply(data_points, data.frame), .id='Transition')
  
  if (!(is.null(filter_transitions.include))){
    log_( paste('-- Size of df before transition filtering: ', dim(df)[1], '\n', sep=''), verbose )
    df %>%
      filter( (Transition %in% filter_transitions.include) ) -> df
    log_( paste('---- Size of df after transition filtering: ', dim(df)[1], '\n', sep=''), verbose )
  }
  if (!(is.null(filter_transitions.exclude))){
    log_( paste('-- Size of df before transition filtering: ', dim(df)[1], '\n', sep=''), verbose )
    df %>%
      filter( !(Transition %in% filter_transitions.exclude) ) -> df
    log_( paste('---- Size of df after transition filtering: ', dim(df)[1], '\n', sep=''), verbose )
  }
  if (!(is.null(filter_RT))){
    log_( paste('-- Size of df before RT filtering: ', dim(df)[1], '\n', sep=''), verbose )
    if (length(filter_RT)==2){
      df %>%
        filter( RT > filter_RT[1] ) %>%
        filter( RT < filter_RT[2] ) -> df
    } else if (length(filter_RT)==1){
      df %>%
        filter( RT > ( peptide_subset$leftWidth - filter_RT) ) %>%
        filter( RT < ( peptide_subset$rightWidth + filter_RT)  ) -> df
    } else{
      log_( paste(red('You entered an invalid value for "filter_RT argument"!!!')), verbose )
    }
    log_( paste('---- Size of df after RT filtering: ', dim(df)[1], '\n', sep=''), verbose )
  }
  
  g <- ggplot(df, aes(RT, Int, col=Transition)) + 
    geom_line(show.legend = FALSE) + 
    geom_vline(xintercept = peptide_subset$RT, color='red' ) +
    geom_vline(xintercept = peptide_subset$leftWidth, color='red', linetype='dotted' ) +
    geom_vline(xintercept = peptide_subset$rightWidth, color='red', linetype='dotted' ) +
    ggtitle(  peptide_subset$FullPeptideName ) +
    labs(subtitle = paste('Run: ', run, ' | Precursor: ', sapply(str_split(precursor_id, '_'), "[[", 2), ' | Charge: ', peptide_subset$Charge, ' | m/z: ', peptide_subset$mz, ' | RT: ', peptide_subset$RT, ' | m-score: ', peptide_subset$m_score, ' | ipf_pep: ', peptide_subset$ipf_pep, sep='')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 14)) +
    theme_bw()
  # print(g)
  if(return_out=='plot'){
    return(g)
  } else if (return_out=='data'){
    return(df)
  } else {
    cat(red('You did not enter the correct value for return_out!! Options {"plot","data"}'), '\nYour input: ', underline(return_out), '\n',sep='')
  }
}