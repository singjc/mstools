getXIC_ <- function( graphic_obj, 
                     df_lib, mod, 
                     in_sqMass, 
                     transition_type='none', 
                     intersecting_mz=NULL, 
                     uni_mod_list=NULL, 
                     max_RT, 
                     min_RT, 
                     max_Int, 
                     in_osw=NULL, 
                     smooth_chromatogram=TRUE,
                     doFacetZoom=FALSE, 
                     FacetFcnCall=NULL,
                     top_trans_mod_list=NULL, 
                     Isoform_Target_Charge=NULL, 
                     RT_pkgrps=NULL,
                     show_manual_annotation=NULL,
                     plotIdentifying.Unique=NULL, 
                     plotIdentifying.Shared=NULL, 
                     plotIdentifying.Against=NULL,
                     show_n_transitions=NULL,
                     show_legend=T,
                     verbosity=0
){
  
  # Check Dataframe, if empty return an object and stop function call
  checkDataframe <- function( df, return_item, msg ){
    if ( dim(df)[1]==0 ){
      cat(bold(red(msg)))
      print(df)
      return( TRUE )
    } else {
      return( FALSE )
    }
  }
  # Check numeric values if not NULL
  checkNumeric <- function( numeric_obj ){
    if( !is.null( numeric_obj )){
      return( formatC(numeric_obj, format = "e", digits = 3) )
    } else {
      return( NULL )
    }
  }
  ## Filter library data for specific information on transition type selected
  if ( transition_type=='precursor' ){
    cat( Verbose(threshold = verbosity), '--> Extracting Precursor Transition...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( TYPE=="" ) %>%
      dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge )-> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for precursor transition in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting_intersection' ){
    cat( Verbose(threshold = verbosity), '--> Extracting Intersecting Detecting Transitions...\n')
    if ( !is.null( intersecting_mz ) ){
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==1 ) %>%
        dplyr::filter( PRODUCT_MZ %in% intersecting_mz ) -> df_lib_filtered
    } else {
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==1 ) -> df_lib_filtered
    }
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for common detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting_unique' ){
    cat( Verbose(threshold = verbosity), '--> Extracting Unique Detecting Transitions...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( DETECTING==1 ) %>%
      dplyr::filter( !(PRODUCT_MZ %in% intersecting_mz) )-> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for unique detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting' ){
    cat( Verbose(threshold = verbosity), '--> Extracting Unique Detecting Transitions...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( DETECTING==1 ) -> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='identifying' ){
    cat( Verbose(threshold = verbosity), '--> Extracting Identifying Transitions...\n')
    if ( !(is.null(top_trans_mod_list)) ){
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==0 ) -> df_lib_filtered
    } else {
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==0 ) -> df_lib_filtered
      # dplyr::filter( (TYPE==uni_mod_list[[1]][['Type']] & (ORDINAL %in% c(uni_mod_list[[1]][['site_determining_start']]:uni_mod_list[[1]][['site_determining_end']])) |
      #                   (TYPE==uni_mod_list[[2]][['Type']] & (ORDINAL %in% c(uni_mod_list[[2]][['site_determining_start']]:uni_mod_list[[2]][['site_determining_end']])) ))) %>%
      # dplyr::filter( !(grepl('*H3O4P1*', TRAML_ID)) ) %>%
      # # dplyr::filter( gsub('.*\\{|\\}.*','',TRAML_ID) %in% paste('',gsub('UniMod:21','Phospho',mod),'',sep='') )
      # dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), gsub('.*\\{|\\}.*','',TRAML_ID)) ) -> df_lib_filtered
    }
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for identifying transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else {
    ## OSW Information
    if ( !is.null( in_osw ) ){
      cat( Verbose(threshold = verbosity), '--> Extracting OpenSwathResults Info...\n')
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==1 ) -> df_lib_filtered
      run_name <- gsub('_osw_chrom[.]sqMass$', '', basename(in_sqMass))
      run <- gsub('_SW*|_SW_0|(*_-_SW[.]mzML[.]gz)', '', gsub('yanliu_I170114_\\d+_|chludwig_K150309_|lgillet_L\\d+_\\d+-Manchester_dirty_phospho_-_', '', run_name))
      mod_position <- getModificationPosition_(mod, character_index = T)+1
      mod_form_rename <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', mod)))))
      ## Filter for precursor id that matches target charge state
      df_lib_filtered %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        select( PRECURSOR_ID ) %>%
        unique() %>% as.matrix() %>% as.numeric() -> target_charge_precursor
      
      ## @TODPO: Need to make this more robust for later
      if( mod==mod_form_rename ){
        # Actual Modification Name Convention
        osw_df <- mstools::getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id='', mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
        # # Original OSW Peptide Names 
        # osw_pep_names <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', osw_df$FullPeptideName)))))
        # # Keep only Rows that correspond to the correct Assay
        # osw_df %>% dplyr::filter( osw_pep_names == osw_df$ipf_FullPeptideName )
      } else {
        # UniMod Convention
        #osw_df <- mstools::getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id=c(mod,mod_form_rename), mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
        osw_df <- mstools::getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id=c(mod,mod_form_rename), mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
      }
      ## Check if openswath dataframe is empty
      if ( checkDataframe( osw_df, graphic_obj, msg='There was no data found in OpenSwath Dataframe\n' ) ){ 
        graphic_obj <- graphic_obj +
          ggtitle(  mod ) +
          labs(subtitle = paste('Run: ', run, 
                                ' | Precursor: ', df_lib_filtered$PRECURSOR_ID, 
                                ' | Peptide: ', df_lib_filtered$PEPTIDE_ID, 
                                ' | Charge: ', df_lib_filtered$PRECURSOR_CHARGE, sep=''))
        return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) 
      }
      if ( !is.null(Isoform_Target_Charge) ){
        osw_df %>%
          dplyr::filter( Charge==Isoform_Target_Charge ) -> osw_df
      }
      if ( !is.null( unlist(osw_df$ipf_pep) ) ){
        # Remove rows with NULL value in ipf_pep
        osw_df %>%
          dplyr::filter( !is.null(ipf_pep) ) %>%
          dplyr::filter( !is.nan(ipf_pep) ) -> osw_df
        
        osw_df %>%
          dplyr::filter( ipf_pep==min(ipf_pep) ) -> osw_df_filtered #### No longer filter by peak group
        
        if ( dim(osw_df_filtered)[1] > 1 ){
          osw_df %>%
            dplyr::filter( peak_group_rank==min(peak_group_rank) ) -> osw_df_filtered #### No longer filter by peak group
        }
        
      } else {
        osw_df %>%
          dplyr::filter( peak_group_rank==1 ) -> osw_df_filtered #### No longer filter by peak group
      }
      m_score <- checkNumeric( osw_df_filtered$ms2_m_score[[1]] )
      ipf_pep <- checkNumeric( osw_df_filtered$ipf_pep[[1]] )
      ipf_m_score <- checkNumeric( osw_df_filtered$m_score[[1]] )
      ms2_pkgrp_rank <- osw_df_filtered$peak_group_rank[[1]]
      
      graphic_obj <- graphic_obj +
        geom_vline(xintercept = osw_df_filtered$RT, color='red', size = 1.5 ) +
        geom_vline(xintercept = osw_df_filtered$leftWidth, color='red', linetype='dotted', size=1 ) +
        geom_vline(xintercept = osw_df_filtered$rightWidth, color='red', linetype='dotted', size=1 ) + # default size 0.65
        ggtitle(  mod ) +
        labs(subtitle = paste(
          'Run: ', run,
          ' | Precursor: ', df_lib_filtered$PRECURSOR_ID,
          ' | Peptide: ', df_lib_filtered$PEPTIDE_ID,
          ' | Charge: ', osw_df_filtered$Charge,
          ' | m/z: ', osw_df_filtered$mz, 
          ' | RT(s): ', osw_df_filtered$RT, 
          ' | RT(m) ', round((osw_df_filtered$RT/60), digits = 2), 
          '\nlib_RT: ', round(osw_df_filtered$assay_iRT, digits = 4),
          ' | ms2-m-score: ', m_score,
          ' | ipf-m-score: ', ipf_m_score,
          ' | ipf_pep: ', ipf_pep,
          ' | ms2_pkgrp_rank: ', ms2_pkgrp_rank,
          sep='')) 
      ### Plot Other peak rank groups
      if( !is.null(RT_pkgrps) ){
        osw_RT_pkgrps <- mstools::getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id=c(mod,mod_form_rename), mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
        
        osw_RT_pkgrps %>%
          dplyr::filter( RT %in% RT_pkgrps ) %>%
          dplyr::filter( ipf_FullPeptideName==mod_form_rename ) %>%
          dplyr::filter( RT != osw_df_filtered$RT ) %>%
          dplyr::select( RT, leftWidth, rightWidth, peak_group_rank, ms2_m_score, ipf_pep, m_score  ) -> osw_RT_pkgrps_filtered
        if ( dim(osw_RT_pkgrps_filtered)[1]!=0 ){
          # Define unique set of colors to annotate different peak rank groups
          jBrewColors <- brewer.pal(n = dim(osw_RT_pkgrps_filtered)[1], name = "Dark2")
          # Y intcrements
          y_increment = 0
          for( RT_idx in seq(1, dim(osw_RT_pkgrps_filtered)[1],1) ){
            # if ( any(osw_RT_pkgrps_filtered$RT[RT_idx] %in% c(5691.37, 5737.26)) ){ cat('Skipping: ', osw_RT_pkgrps_filtered$RT[RT_idx], '\n', sep=''); next }
            
            ## get scores and format them
            rank <- osw_RT_pkgrps_filtered$peak_group_rank[RT_idx]
            ms2_m_score <- formatC(osw_RT_pkgrps_filtered$ms2_m_score[RT_idx], format = "e", digits = 3)
            ipf_pep <- formatC(osw_RT_pkgrps_filtered$ipf_pep[RT_idx], format = "e", digits = 3)
            ipf_m_score <- formatC(osw_RT_pkgrps_filtered$m_score[RT_idx], format = "e", digits = 3)
            
            point_dataframe <- data.frame(RT=(osw_RT_pkgrps_filtered$RT[RT_idx]),
                                          RT_m = round((osw_RT_pkgrps_filtered$RT[RT_idx])/60, digits = 2),
                                          # y=((max(max_Int)/ 4 )-y_increment),
                                          rank = rank,
                                          ms2_m_score = ms2_m_score,
                                          ipf_pep = ipf_pep,
                                          ipf_m_score = ipf_m_score
                                          # # y=((1000)-y_increment),
                                          # label=paste('Rank:',rank,'\n',
                                          #             osw_RT_pkgrps_filtered$RT[RT_idx], '\n',
                                          #             'ms2-m-score:', ms2_m_score, '\n',
                                          #             'ipf_pep:', ipf_pep, '\n',
                                          #             'ipf-m-score:', ipf_m_score, 
                                          #             sep=' ')
            )
            
            graphic_obj <- graphic_obj +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$RT[RT_idx], color=jBrewColors[RT_idx], alpha=0.65, size = 1.5 ) +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$leftWidth[RT_idx], color=jBrewColors[RT_idx], linetype='dotted', alpha=0.85, size=1 ) +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$rightWidth[RT_idx], color=jBrewColors[RT_idx], linetype='dotted', alpha=0.85, size=1 ) 
            # geom_rect( aes(xmin=osw_RT_pkgrps_filtered$leftWidth[RT_idx], xmax=osw_RT_pkgrps_filtered$rightWidth[RT_idx], ymin=0, ymax=Inf), col=jBrewColors[RT_idx], fill=jBrewColors[RT_idx], alpha=0.5) +
            # geom_label(data=point_dataframe, aes(x=RT, y=y,label=label), alpha=0.7, fill=jBrewColors[RT_idx], size=3)
            y_increment = y_increment + 500
            
            #****************************************************
            # Check to see if master annotation table exits
            #****************************************************
            if ( !exists("master_annotation_table") ){
              master_annotation_table <- point_dataframe
            } else {
              master_annotation_table <- rbind(master_annotation_table, point_dataframe)
            }
            
          }
          
          theme_col <- gridExtra::ttheme_default(core = list(fg_params = list( col = c( jBrewColors )),
                                                             bg_params = list( col = NA)),
                                                 rowhead = list(bg_params = list(col = NA)),
                                                 colhead = list(bg_params = list(col = NA))
          )       
          
          #*****************************************************
          # Make tableGrob object with annotation information
          # and get height of table
          #*****************************************************
          annotationGrob <- tableGrob(master_annotation_table, theme = theme_col, row = NULL) 
          th <- sum(annotationGrob$heights)
          
        }
      }
     
      #***************************************
      # Plot Mannual Annotation Coordinates
      #***************************************
      if ( !is.null(show_manual_annotation) ){
        graphic_obj <- graphic_obj + 
          geom_rect(data = data.frame(xmin = show_manual_annotation[1],
                                      xmax = show_manual_annotation[2],
                                      ymin = -Inf,
                                      ymax = Inf),
                    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                    fill = "blue", alpha = 0.1)
        
      }
      
      if ( doFacetZoom==TRUE ){
        if ( !is.null(FacetFcnCall) ){
          graphic_obj <- graphic_obj + FacetFcnCall
        } else {
          if ( max(max_Int) > 1000 ){
            graphic_obj <- graphic_obj +
              # facet_zoom(ylim = c(0, (1000) ))
              facet_zoom(ylim = c(0, (max(max_Int)/ 4 ) ))
            # facet_zoom(ylim = c(0, (max(max_Int)/ (max(max_Int)-mean(max_Int)) ) ))
          }
        }
      }
      if ( !is.null(top_trans_mod_list) ){
        graphic_obj <- graphic_obj +
          theme(plot.title = element_text(hjust = 0.5), 
                plot.subtitle = element_text(hjust = 0.5, size = 10),
                legend.text = element_text(size = 2),
                legend.key.size = unit(0.5, "cm")) +
          # theme( panel.background = element_rect( fill='lightblue', color='lightblue', size=0.5, linetype='solid'),
          #       panel.border = element_rect( fill=NA, color='black', size=0.5, linetype='solid'),
          #       panel.grid.major = element_line( color='white', size=0.5, linetype='solid'),
          #       panel.grid.minor = element_line( color='white', size=0.25, linetype='solid') ) 
          theme_bw()
        
        if ( show_legend==FALSE ){
          graphic_obj <- graphic_obj + theme(legend.position="none") 
        }
        
      } else {
        graphic_obj <- graphic_obj +
          theme(plot.title = element_text(hjust = 0.5), 
                plot.subtitle = element_text(hjust = 0.5, size = 10)) +
          # theme(panel.background = element_rect( fill='#c1e1ec', color='#c1e1ec', size=0.5, linetype='solid'),
          #       panel.border = element_rect( fill=NA, color='black', size=0.5, linetype='solid'),
          #       panel.grid.major = element_line( color='white', size=0.5, linetype='solid'),
          #       panel.grid.minor = element_line( color='white', size=0.25, linetype='solid') ) 
          # guides( Transition=FALSE ) +
          theme_bw()
        if ( show_legend==FALSE ){
          graphic_obj <- graphic_obj + theme(legend.position="none") 
        }
      }
      
      if ( exists("annotationGrob") ){
        graphic_obj <- ggpubr::as_ggplot( arrangeGrob( graphic_obj, annotationGrob, nrow = 2, heights = unit.c(unit(1, "null"), th )) )
      }
      return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
    } else {
      cat(red( bold(underline(transition_type)), ' is not a supported argument for transition_type!!!\n'), sep='')
      return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
    }
  }
  
  ## Get Transition IDs for chromatogram data extraction
  cat( Verbose(threshold = verbosity), '----> Getting Transition IDs and Extracting Chromatogram Data...\n')
  frag_ids <- list()
  if( length(as.character( df_lib_filtered$TRANSITION_ID ))>1 ){
    frag_ids[[1]] <- as.character( df_lib_filtered$TRANSITION_ID )
  } else {
    frag_ids[[1]] <- list(as.character( df_lib_filtered$TRANSITION_ID ))
  }
  chrom <- getChromatogramsbyIndice_( in_sqMass, frag_ids )
  if (smooth_chromatogram==TRUE){
    cat( Verbose(threshold = verbosity), '----> Smoothing Chromatogram Data...\n')
  }
  ## Smooth Intensity values to make Chromatogram look nice
  for (i in seq(1:length(chrom))){
    names(chrom[[i]]) <- c('RT','Int')
    if (smooth_chromatogram==TRUE){
      chrom[[i]]$Int <- signal::sgolayfilt(chrom[[i]]$Int, p = 4, n = 9)
    }
  }
  
  df_plot <- bind_rows(mclapply(chrom, data.frame), .id='Transition')
  
  if ( transition_type=='precursor' ){
    df_lib_filtered %>% 
      dplyr::filter( TRANSITION_ID %in% unique(df_plot$Transition) ) %>%
      select( TRANSITION_ID, PRECURSOR_CHARGE, PRECURSOR_MZ ) -> transition_info
    transition_ids <- paste( paste('0 - ',transition_info$TRANSITION_ID,sep=''), paste(transition_info$PRECURSOR_CHARGE,'+',sep=''), transition_info$PRECURSOR_MZ, sep='_')
    tmp <- sapply(seq(1,length(df_plot$Transition)), function(i){ transition_ids[grepl(paste('0 - ',df_plot$Transition[i],'_*',sep=''), transition_ids)] } )
  } else {
    # if ( transition_type=='identifying' ){
    #   if( !(is.null(top_trans_mod_list)) ){
    #     cat( '----> Extracting Top Transitions with low PEPs..\n')
    #     df_plot %>%
    #       dplyr::filter( Transition %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1])[1:10] ) -> df_plot
    #   }
    # }
    df_lib_filtered %>% 
      dplyr::filter(  TRANSITION_ID %in% unique(df_plot$Transition) ) %>%
      select( TRANSITION_ID, CHARGE, TYPE, ORDINAL, PRODUCT_MZ ) -> transition_info
    transition_ids <- paste( transition_info$TRANSITION_ID, paste(transition_info$CHARGE,'+',sep=''), paste(transition_info$TYPE, transition_info$ORDINAL,sep=''), transition_info$PRODUCT_MZ, sep='_')
    tmp <- sapply(seq(1,length(df_plot$Transition)), function(i){ transition_ids[grepl(paste('^',df_plot$Transition[i],'_*',sep=''), transition_ids)] } )
  }
  df_plot$TRANSITION_ID <- df_plot$Transition
  df_plot$Transition <- tmp
  
  # if ( !is.null(min_RT) & !is.null(max_RT) ){
  #   df_plot %>%
  #     dplyr::filter( RT > min_RT ) %>%
  #     dplyr::filter( RT < max_RT ) -> df_plot
  # }
  # print(unique(df_plot$Transition))
  ## Check Max Intensity
  if ( transition_type!='precursor' & dim(df_plot)[1]>0 ){
    # if ( max(df_plot$Int) > max_Int ){ max_Int <- max(df_plot$Int) }
    max_Int <- c(max_Int, max(df_plot$Int[is.finite(df_plot$Int)]))
  }
  # ## Show a legend @TODO: Make this a top level argument.
  # if ( !is.null(top_trans_mod_list) ){
  #   show_legend = T
  # } else {
  #   show_legend = F
  # }
  
  ## Plotting Action
  if ( transition_type=='precursor' ){
    graphic_obj <- graphic_obj + 
      geom_line( data=df_plot, aes(RT, Int, group=Transition, alpha=0.65, text=paste('Transition: ', Transition, sep='')), show.legend = show_legend, linetype='solid', col='black' )  +
      scale_alpha_identity(name="Precursor", guide='legend', labels=unique(df_plot$Transition))
  } else if ( transition_type=='detecting_intersection' ){
    graphic_obj <- graphic_obj + 
      geom_line(data=df_plot, aes(RT, Int, group=Transition, fill=Transition, text=paste('Transition: ', Transition, sep='')), show.legend = show_legend, alpha=0.95, col='gray')  +
      # scale_fill_manual(name="Detecting",values=rep('gray', length(unique(df_plot$Transition)))) 
      scale_fill_identity(name="Detecting", guide='legend', labels=unique(df_plot$Transition))
  } else if ( transition_type=='detecting_unique' ){
    graphic_obj <- graphic_obj + 
      geom_line(data=df_plot, aes(RT, Int, group=Transition, text=paste('Transition: ', Transition, sep='')), alpha=0.75, col='black', linetype='dotted', show.legend = show_legend) 
    # theme( legend.position = 'bottom' )
  } else if ( transition_type=='detecting' ){
    graphic_obj <- graphic_obj + 
      geom_line(data=df_plot, aes(RT, Int, group=Transition, fill=Transition, text=paste('Transition: ', Transition, sep='')), alpha=0.75, col='gray', linetype='solid', show.legend = show_legend) +
      scale_fill_manual(name="Detecting",values=rep('gray', length(unique(df_plot$Transition))) )
  } else{
    
    df_plot <- merge( df_plot, select( df_lib_filtered[df_lib_filtered$TRANSITION_ID %in% unique(df_plot$TRANSITION_ID),], c(TRANSITION_ID, TRAML_ID) ), by='TRANSITION_ID' )
    
    uni_forms <- unique (unlist( strsplit( unique(gsub(".*\\{|\\}.*", "", df_plot$TRAML_ID)), split='\\|') ) )
    
    # uni_forms <- c("EGHAQNPM(UniMod:35)EPS(UniMod:21)VPQLSLM(UniMod:35)DVK", "EGHAQNPM(UniMod:35)EPSVPQLS(UniMod:21)LM(UniMod:35)DVK")
    
    peptide_modification_positions <- lapply( uni_forms, function(pep){ mstools::getModificationPosition_( pep ) } )
    names(peptide_modification_positions) <- uni_forms
    
    do.call( rbind, lapply(peptide_modification_positions, function(x){ x[match(names(peptide_modification_positions[[1]]), names(x))]}) ) %>% ## Turn list into dataframe, and ensure items in list are matching by name
      as.data.frame() %>%
      tibble::rownames_to_column( "peptides" ) %>%
      mutate( unimod_peptides = mstools::codenameTounimod(peptides) ) %>%
      mutate( target_peptide = unimod_peptides==mod ) -> peptides_information_df
    
    ## Get modification columns
    modification_col_indices <- grep('modification_.*', colnames(peptides_information_df))
    
    ## convert naked peptide length and modification position columns to nemeric vectors
    lapply(modification_col_indices, function(col){
      set( peptides_information_df, NULL, col, as.numeric(peptides_information_df[[col]]) )
    })
    class(peptides_information_df$naked_peptide_length) <- "numeric"
    
    if ( length(modification_col_indices)>1 ){
    modifications_position_df <- peptides_information_df[,modification_col_indices] 
    
    target_modification <- apply(modifications_position_df, 2, function( positions ){
      if( length(unique(positions)) > 1){
        return( TRUE )
      } else {
        return( FALSE )
      }
    })
    
    target_modification <- names(target_modification)[target_modification]
    } else {
      target_modification <- grep('modification_.*', colnames(peptides_information_df), value = TRUE)
    }
    
    peptides_information_df %>%
      dplyr::mutate( 
        modification_position = ifelse( !!rlang::sym(target_modification)==max( !!rlang::sym(target_modification) ), 'most_right',
                                        ifelse( !!rlang::sym(target_modification)==min( !!rlang::sym(target_modification) ), 'most_left',
                                        ifelse( !!rlang::sym(target_modification)!=max( !!rlang::sym(target_modification) | !!rlang::sym(target_modification)!=max( !!rlang::sym(target_modification) ) ), 'inbetween', 'nan') 
                                      ) )
        ) -> peptides_information_df
    
    ### Get series starting from up to modification amino acid
    
    peptides_information_df %>%
      dplyr::mutate(
        
        ion_series_start = ifelse( modification_position=='most_right' &  target_peptide, 
                                   dplyr::filter( dplyr::select(peptides_information_df, !!rlang::sym(target_modification)), peptides_information_df$unimod_peptides!=mod )+1, 
                                   dplyr::filter( dplyr::select(peptides_information_df, !!rlang::sym(target_modification)), peptides_information_df$unimod_peptides==mod )-1
          
        ),
        ion_series_end = ifelse( modification_position=='most_right' &  target_peptide, 
                                 dplyr::filter( dplyr::select(peptides_information_df, !!rlang::sym(target_modification)), peptides_information_df$unimod_peptides==mod ), 
                                 dplyr::filter( dplyr::select(peptides_information_df, !!rlang::sym(target_modification)), peptides_information_df$unimod_peptides!=mod )
                                 
        ),
        use_ion_series = ifelse( modification_position=='most_right', 'y', 'b')
        
        
      ) -> peptides_information_df
    class(peptides_information_df$ion_series_start) <- "numeric"
    class(peptides_information_df$ion_series_end) <- "numeric"
    
    if ( !(is.null(top_trans_mod_list)) ){
      
      #### Unique Transitions
      if ( plotIdentifying.Unique==TRUE ){
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
                               gsub('.*\\{|\\}.*','',TRAML_ID)) &
                           (!grepl(glob2rx('*\\|*'), 
                                   gsub('.*\\{|\\}.*','',TRAML_ID))) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1]) ) -> tmp_plot
        
        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx("*_\\{ESTAEPDSLS(Phospho)R(Label:13C(6)15N(4))\\}_*"), TRAML_ID) ) -> tmp_plot
        # 
        # 
        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx("*_\\{ESTAEPDS(Phospho)LSR(Label:13C(6)15N(4))\\}_*"), TRAML_ID) ) -> tmp_plot
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        
        
        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) +
            guides(col=guide_legend(title="Identifying")) #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
          #       
          
          # graphic_obj + geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) + facet_zoom(ylim = c(0, 5000))
          
        }
      }
      #### Shared Transitions
      if ( plotIdentifying.Shared==TRUE ){
        # if ( transition_type=='identifying' ){
        #   if( !(is.null(top_trans_mod_list)) ){
        #     df_plot <- df_plot_org
        #     df_plot %>%
        #       dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1])[1:10] ) -> df_plot
        #   }
        # }
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
                               gsub('.*\\{|\\}.*','',TRAML_ID)) &
                           grepl(glob2rx('*\\|*'), 
                                 gsub('.*\\{|\\}.*','',TRAML_ID)) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<0.6]) ) -> tmp_plot
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='F1', alpha=0.5, show.legend = show_legend)  +
            guides(col=guide_legend(title="Identifying")) #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
      ##### Transitions against
      if ( plotIdentifying.Against==TRUE ){
        df_plot %>%
          dplyr::filter( !grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
                                gsub('.*\\{|\\}.*','',TRAML_ID)) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[ top_trans_mod_list[[mod]]$transition_pep<1 ]) ) -> tmp_plot
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='dashed', show.legend = show_legend) +
            guides(col=guide_legend(title="Identifying"))  #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
    } else {
      #### Unique Identifying Transitions
      if ( plotIdentifying.Unique==TRUE ){
        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
        #                        gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
        #                    (!grepl(glob2rx('*\\|*'), 
        #                            gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID)))) ) -> tmp_plot

        # 
        # mods_present <- names(uni_mod_list)
        # alternate_mod <- mods_present[ !( mods_present %in% mod)]
        # alternate_mod_index <- unlist(lapply(alternate_mod, getModificationPosition_))
        # current_mod_index <- getModificationPosition_(mod)
        # 
        # identification_traml_ids <- str_split( gsub('.*\\{|\\}.*', '', df_plot$TRAML_ID), '\\|' )
        # current_mod_transitions <- unlist(lapply( identification_traml_ids, function( traml_id ){
        #   mods_for_ident_transition <- unique(unlist(lapply(traml_id, getModificationPosition_)))
        #   if ( any( current_mod_index == mods_for_ident_transition ) ){
        #     if ( any( (mods_for_ident_transition %in% alternate_mod_index) ) ){
        #       return( FALSE )
        #     } else {
        #       return( TRUE )
        #     }
        #   } else {
        #     return( FALSE )
        #   } 
        # } ) )
        # 
        # 
        # df_plot %>%
        #   dplyr::filter( current_mod_transitions ) -> df_plot
        # tmp_plot <- df_plot
        
        
        df_plot %>%
          dplyr::filter( 
            grepl( glob2rx( paste( paste0( peptides_information_df$use_ion_series[peptides_information_df$target_peptide], 
                                    base::seq(
                                      from=peptides_information_df$ion_series_end[peptides_information_df$target_peptide], 
                                      to=peptides_information_df$ion_series_start[peptides_information_df$target_peptide], 
                                      by=ifelse(peptides_information_df$use_ion_series[peptides_information_df$target_peptide]=='y',-1,1) 
                                      )  
                                    ), collapse = '|' ) 
                          ), gsub('.*\\{[a-zA-Z0-9\\|\\(\\)]+\\}_\\d+.\\d+_\\d+.\\d+_-\\d+.\\d+_([abcxyz0-9]+).*', '\\1', df_plot$TRAML_ID) 
              
            )
            ) -> tmp_plot
        
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        if ( dim(tmp_plot)[1] > 0){
          
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) +
            guides(col=guide_legend(title="Identifying"))  #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
      ##### Shared Identifying Transitions
      if ( plotIdentifying.Shared==TRUE ){
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
                               gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
                           grepl(glob2rx('*\\|*'), 
                                 gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) ) -> tmp_plot
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        if ( dim(tmp_plot)[1] > 0){
          
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='F1', show.legend = show_legend) +
            guides(col=guide_legend(title="Identifying"))  #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
      
      ##### Transitions Against
      if ( plotIdentifying.Against==TRUE ){
        df_plot %>%
          dplyr::filter( !grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')), 
                                gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',TRAML_ID)) ))  -> tmp_plot
        
        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }
        
        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int
        
        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot
        
        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj + 
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition, text=paste('Transition: ', Transition, sep='')), linetype='dashed', show.legend = show_legend) +
            guides(col=guide_legend(title="Identifying"))  #+ 
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
    }
  }
  
  return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
  
}

