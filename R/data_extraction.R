data_extraction_ <- function( file ) {
  osw.results <- list()
  run_id <- gsub('_SW.mzXML.gz.tsv','',gsub('yanliu_I170114_0\\d+_','',basename(file)))
  cat(blue('o'), ' Working on: ', run_id,'\n',sep='')
  
  osw.results[[run_id]][['osw.out']] <- fread(file)
  # Filter OpenSwath Results
  osw.results[[run_id]][['osw.out']] %>%
    filter( decoy==0 ) %>%
    filter( peak_group_rank==1 ) %>%
    filter( m_score<0.01 ) -> osw.results[[run_id]][['osw.out.filtered']]
  
  osw.results[[run_id]][['frag_id_list']] <- (strsplit(osw.results[[run_id]][['osw.out.filtered']]$aggr_Fragment_Annotation, ';'))
  osw.results[[run_id]][['frag_id_list']] <- mclapply(osw.results[[run_id]][['frag_id_list']], function(x){gsub('_[y,b]\\d+_\\d+$','',x)})
  sqMass_file <- grep(paste('_',run_id,'_',sep=''),list.files(path=paste(getwd(), sep=''), pattern="*.sqMass$", full.names = T, recursive = T),value=T)
  tic('--- Chromatogram Extraction took: ')
  chrom_extract_error <- tryCatch({
    osw.results[[run_id]][['chroms_list']] <- getChromatogramsbyIndice_( sqMass_file, (osw.results[[run_id]][['frag_id_list']]) )
    osw.results[[run_id]][['productmz_list']] <- getChromatogramsMZbyIndice_( sqMass_file, (osw.results[[run_id]][['frag_id_list']]) )
  }, error=function(e){
    e
    cat(red('There was an issue that occured while extracting chromatograms from: '), underline(sqMass_file), '\n', sep='')
  }
  )
  if(inherits(chrom_extract_error, "error")) next
  # index=1
  # for ( precursor in osw.results[[run_id]][['osw.out.filtered']]$transition_group_id ){
  #   # progress(index,max.value=dim(osw.results[[run_id]][['osw.out.filtered']])[1], progress.bar = T)
  #   cat('Index: ',index,', Precursor ', precursor,' of ', dim(osw.results[[run_id]]$osw.out.filtered)[1], '\n', sep='')
  #   osw.results[[run_id]][['chroms_list']][[as.character(precursor)]] <- (getChromatogramsbyIndice_( sqMass_file, (osw.results[[run_id]][['frag_id_list']][[index]]) ))
  #   index=index+1
  # }
  
  # osw.results[[run_id]][['chroms_list']] <- mclapply(seq(1,dim(osw.results[[run_id]]$osw.out.filtered)[1],1), 
  #          function(x){
  #            tmp_list <- list()
  #            tmp_list[[as.character(osw.results[[run_id]][['osw.out.filtered']]$transition_group_id[x])]] <- getChromatogramsbyIndice_( sqMass_file, (osw.results[[run_id]][['frag_id_list']][[x]]) )
  #   
  # }, mc.cores = detectCores()-1)
  
  toc()
  return(osw.results)
}