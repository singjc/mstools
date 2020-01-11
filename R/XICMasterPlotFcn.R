#// **********************************************************************************************
#//                         XICMasterPlotFcn.R
#// **********************************************************************************************
#//
#// 
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing

#' @export
#' @title Master Plotting function for plotting XIC's
#' @description This function can be used to draw XIC's by calling getXIC 
#' 
#' @param dup_peps A character vector of a peptide sequence(s). 
#' @param uni_mod A character vector of a modified peptide sequence. (Default: NULL) If using this argument, dup_peps must be a single peptide
#' @param sqMass_files A list of character vectors. Full paths to chromatogram file(s).
#' @param in_lib A character vector. Full path to a pqp assay library.
#' @param in_osw A character vector. Full path to an osw results file.
#' @param plotPrecursor A logical value. True will plot precursor chromatogram
#' @param plotIntersectingDetecting A logcail value. True will plot intersecting detecting transitions if comparing two peptidoforms.
#' @param plotUniqueDetecting A logical value. True will plot unique detecting transitions if comparing two peptidoforms.
#' @param plotIdentifying A logical value. True will plot identifying transitions.
#' @param plotIdentifying.Unique A logical value. True will plot unique identifying transitions.
#' @param plotIdentifying.Shared A logical value. True will plot shared identifying transitions.
#' @param plotIdentifying.Against A logical value. True will plot against identifying transitions.
#' @param smooth_chromatogram A list containing p (numeric) for the polynomial order for sgolay, and n (numeric) the bandwidth for sgolay. (Defualt: list(p=4, n=9)
#' @param doFacetZoom A logical valie. Should the plot be zoomed in. The default zooming operation is based on the max int divided by 4. (Default: FALSE)
#' @param FacetFcnCall A facet_zoom function with user defined parameters. i.e. FacetFcnCall = facet_zoom(xlim = c(7475, 7620), ylim = c(0, 4000) ). (Default: NULL)
#' @param doPlot A logical value. TRUE will perform steps to save the plot as a pdf.
#' @param Charge_State A numeric value. The target charge state. (Default: NULL)
#' @param N_sample A numeric value. The number of peptidoforms to sample, if more than one. (Default: 1)
#' @param idx_draw_these A vector of numeric values. The indices for the peptidoforms you wish to exaine.
#' @param store_plots_subdir A character vector. The location to store plots.
#' @param printPlot A logical value. TRUE will print plot in RStudio display.
#' @param use_top_trans_pep A logical value. TRUE will rank the transitions based on the posterior error probabilities.
#' @param show_n_transitions A numeric value. Show n number of transitions
#' @param show_all_pkgrprnk A logical value. Show all feature peak-group ranks. Usually 5. (Default: 5)
#' @param show_manual_annotation A dataframe with leftWidth and rightWidth retention time boundary values of a manually annotated peak. Will draw a transparent blue shaded rectangle indicating manual annotation. I.e data.frame(leftWidth=300, rightWidth=330)
#' @param show_legend A logical value. Display legend information for transition id, m/z and charge. (Default: TRUE)
#' 
#' @return A ggplot-grobs table of a XIC
#' 
#' @author Justin Sing \url{https://github.com/singjc}
#' @import tictoc
#' @import crayon
#' @import parallel
#' @import ggplot2
#' @import dplyr
#' @importFrom dplyr %>%
#' @import stringr
XICMasterPlotFcn_ <- function( dup_peps, 
                               uni_mod=NULL, 
                               sqMass_files,  in_lib, in_osw, 
                               plotPrecursor=T,
                               plotIntersectingDetecting=T,
                               plotUniqueDetecting=F,
                               plotIdentifying=F,
                               plotIdentifying.Unique=F,
                               plotIdentifying.Shared=F,
                               plotIdentifying.Against=F,
                               smooth_chromatogram=list(p = 4, n = 9), 
                               doFacetZoom=F,
                               FacetFcnCall=NULL,
                               doPlot=T,
                               Charge_State=NULL,
                               N_sample=1,
                               idx_draw_these=NULL,
                               store_plots_subdir = '/Results/',
                               printPlot=F,
                               use_top_trans_pep=F,
                               show_n_transitions=NULL,
                               show_all_pkgrprnk=T,
                               show_manual_annotation=NULL,
                               show_legend=T
                               ){
  
  # Get XICs for Modified Peptides  ---------------------------------------------------------------
  # 
  
  # verbose <- Verbose(threshold = 0)
  # cat(verbose, red('bleh'))
  
  if ( F ){
    file_idx <- 1
    dup_peps <- "SQGGDGYYGR"
    pep <- dup_peps
    uni_mod <- "SQGGDGYY(UniMod:21)GR"
    Charge_State <- 2
    sqMass_files <-  "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/sqMass/lgillet_L160915_014-Manchester_dirty_phospho_-_Pool_U7_-_SW.mzML.gz_osw_chrom.sqMass"
    in_lib <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/lib/ConsLib_5600_dpp_UK_peakview_curated_optimized_decoys.pqp"
    in_osw <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Christian_Doerig_Dataset/results/U_pools_sgolay_24/pyprophet_parametric/merged_runs_group_id_MS1MS2_intergration.osw"
    plotPrecursor=T;
    plotIntersectingDetecting=T;
    plotUniqueDetecting=F;
    plotIdentifying=T;
    plotIdentifying.Unique=T;
    plotIdentifying.Shared=F;
    plotIdentifying.Against=F;
    smooth_chromatogram=list(p = 4, n = 9) 
    doFacetZoom=T
    FacetFcnCall=NULL
    doPlot=T
    RT_padding=100
    N_sample=1
    idx_draw_these=NULL
    store_plots_subdir = '/Results/'
    printPlot=F
    use_top_trans_pep=F
    show_n_transitions=NULL
    show_all_pkgrprnk=T
    show_manual_annotation=NULL
    show_legend=T
  }
  
  tictoc::tic( paste('XIC plotting for ', length(dup_peps), ' peptides took: ', sep=' '))
  pep_counter = 0
  for ( pep in dup_peps){
    pep_counter = pep_counter + 1
    cat( crayon::blue('#', pep, ' | (',pep_counter, ' of ', length(dup_peps), ')\n', sep='') )
    # record_list <- list()
    # record_i dx <- 1
    # for ( file_idx in seq(1,length(sqMass_files),1) ){
    record_list <- parallel::mclapply( seq(1,length(sqMass_files),1), function(file_idx){
      in_sqMass <- sqMass_files[file_idx]
      run_name <- gsub('_osw_chrom[.]sqMass$|[.]chrom.mzML$|[.]chrom.sqMass$', '', basename(in_sqMass))
      run <- gsub('_SW*|_SW_0|(*_-_SW[.]mzML[.]gz|[.]chrom[.]sqMass)', '', gsub('yanliu_I170114_\\d+_|chludwig_K150309_|lgillet_L\\d+_\\d+-Manchester_dirty_phospho_-_', '', run_name))
      
      cat(crayon::blue('@ Run: ', run),'\n', sep='')
      
      cat('  o Working on: ', pep, '\n', sep='')
      plot_chrom_error <- tryCatch({
        
        #####################################
        ## TIME LIMIT FOR SCRIPT EXECUTION ##
        #####################################
        ## Time limit in seconds
        time_limit <- 4000 
        
        setTimeLimit(cpu = time_limit, elapsed = time_limit, transient = TRUE)
        on.exit({
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        })
        
        
        
        cat('   ~ Getting peptide library data.. for ', pep, '\n', sep='')
        # Retrieve library data for specific peptide
        df_lib <- mstools::getPepLibData_( in_lib, peptide_id=pep )
        
        cat('   ~ Getting OpenSwath data.. for ', pep, '\n', sep='')
        # Load OSW Merged df
        osw_df <- mstools::getOSWData_( in_osw, run_name, precursor_id='', peptide_id=pep, mod_residue_position='', peak_group_rank_filter=T, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=F )
        if ( dim(osw_df)[1]==0 ){ cat(crayon::red(pep, ' was not found as a peak_rank_group=1 in osw file!!!, skipping...\n'),sep=''); return(list()) }
        
        # Get unique number of modifications
        uni_mod_precursor_id <- unique(paste( df_lib$MODIFIED_SEQUENCE, df_lib$PRECURSOR_ID, df_lib$PRECURSOR_CHARGE, sep='_'))
        
        # Get the unique desired modifications found in openswath/pyprophet results
        desired_uni_mods <- grep(paste(paste(osw_df$id_precursor, osw_df$Charge, sep='_'), collapse = '|'), uni_mod_precursor_id, value=T)
        
        # Get the charge states of the uni modifications found in OSW results
        current_uni_mod_charges <- as.numeric(gsub('.*_', '', desired_uni_mods))
        # if ( length(current_uni_mod_charges)==2 & (current_uni_mod_charges[1]!=current_uni_mod_charges[2]) ){ cat(crayon::red('The two isoforms for:', underline(pep), 'are of different precursor charge. Skipping...\n'), sep=' '); return(list()) }
        
        # Get a list of unique modifications
        #### @TODO: make this more versatile for cases with more than 2 isoforms
        if (is.null(uni_mod)){
          uni_mod <- gsub('_\\d+_\\d+$', '', desired_uni_mods[grep( paste('_\\d+_',mstools::Mode(as.numeric(gsub('.*_', '', desired_uni_mods))), sep=''), desired_uni_mods)])
          uni_mod <- sort(uni_mod)
          # uni_mod <- uni_mod[1:2]
          if(length(uni_mod)>2 | N_sample==1){
            cat(crayon::red('There are more than 2 Modification forms for', crayon::underline(pep),'\n Currently cannot handle more than 2 peptidoforms...\nPetidoforms:\n', paste(uni_mod,collapse='\n'),'\n\n', sep=''))
            
            cat('Will randomly sample 2 of the ', length(uni_mod), ' peptidoforms to process...\n')
            n_mod_sites <- str_count( uni_mod, '\\(UniMod:21\\)|\\(Phospho\\)' ) + (str_count( uni_mod, '\\(UniMod:35\\)|\\(Oxidation\\)' )*3) + (str_count( uni_mod, '\\(UniMod:4\\)|\\(Carbamidomethyl\\)' )*6)
            n_mod_sites_common <- mstools::Mode(n_mod_sites)
            cat( crayon::magenta$underline('There is(are) ', crayon::bold(length(n_mod_sites_common)), ' group(s) of isoforms...\n'))
            if ( !is.null(idx_draw_these) ){
              uni_isoform_group_list <- list()
              uni_isoform_group_list[[1]] <- uni_mod[idx_draw_these]
            } else {
              if ( length(uni_mod)==1 ){ n_sample=1 } else { n_sample=N_sample }
              # if ( length(unique(n_mod_sites_common))!=(length(uni_mod))/2 ){ cat(crayon::red('These are not isoforms of each other... Skipping... \n')); return(list()) }
              uni_isoform_group_list <- lapply(n_mod_sites_common, function(isoform_group){ sample(uni_mod[grepl(paste('^',isoform_group,'$',sep=''), n_mod_sites)], n_sample) })
              uni_isoform_group_list[[1]] <- sort(uni_isoform_group_list[[1]])
            }
            
          } else {
            n_mod_sites <- str_count( uni_mod, '\\(UniMod:21\\)|\\(Phospho\\)' ) + (str_count( uni_mod, '\\(UniMod:35\\)|\\(Oxidation\\)' )*3) + (str_count( uni_mod, '\\(UniMod:4\\)|\\(Carbamidomethyl\\)' )*6)
            n_mod_sites_common <- mstools::Mode(n_mod_sites)
            if ( length(uni_mod)==1 ){ mod=uni_mod; max_Int=0; return( drawNakedPeptide_(df_lib=df_lib, mod=mod, pep=pep, in_sqMass=in_sqMass, plotPrecursor=plotPrecursor, plotIntersectingDetecting=plotIntersectingDetecting,  plotIdentifying=plotIdentifying, plotUniqueDetecting=plotUniqueDetecting, plotIdentifying.Unique=plotIdentifying.Unique, plotIdentifying.Shared=plotIdentifying.Shared, plotIdentifying.Against=plotIdentifying.Against, intersecting_mz=NULL, uni_mod_list=NULL, max_Int=max_Int, in_osw=in_osw, smooth_chromatogram=smooth_chromatogram, doFacetZoom=doFacetZoom, top_trans_mod_list=NULL, show_all_pkgrprnk=show_all_pkgrprnk, FacetFcnCall=FacetFcnCall, show_legend = show_legend ) ) } 
            if ( length(uni_mod)==1 ){ n_sample=1 } else { n_sample=2 } # @TODO: Need to figure something out for this and the line above...
            if ( length(unique(n_mod_sites_common))!=(length(uni_mod))/2 ){ cat(crayon::red('These are not isoforms of each other... Skipping... \n')); return(list()) }
            uni_isoform_group_list <- lapply(n_mod_sites_common, function(isoform_group){ sample(uni_mod[grepl(paste('^',isoform_group,'$',sep=''), n_mod_sites)], n_sample) })
            
          }
        }else {
          uni_isoform_group_list <- uni_mod
        }
        # modification_types_found <- unique(unlist(str_extract_all(uni_mod_list, '(UniMod:\\d+)')))
        len_unmod <- nchar( unique(df_lib$UNMODIFIED_SEQUENCE) )
        final_g_list <- list()
        for ( uni_isoform_group_list_idx in seq(1,length(uni_isoform_group_list),1) ){
          skipped_bool=FALSE
          uni_mod <- uni_isoform_group_list[[uni_isoform_group_list_idx]]
          uni_mod <- stringr::str_replace_all(uni_mod, "\t|\n", "")
          osw_df %>%
            dplyr::filter( FullPeptideName %in% uni_mod) %>%
            select( Charge ) %>%
            as.matrix() %>%
            mstools::Mode() -> Isoform_Target_Charge
          
          if ( length(Isoform_Target_Charge)>1 ){
            
            if ( is.null(Charge_State) ){
              cat('* There are ', length(Isoform_Target_Charge), ' charge states.. Will Randomly sample one of the charge states...\n', sep='')
              Isoform_Target_Charge <- sample(Isoform_Target_Charge,1) # @ CHANGEEE
            } else {
              Isoform_Target_Charge <- Charge_State
            }
            cat('* Will analyze peptidoform with charge state: ', (Isoform_Target_Charge), '\n', sep='')
          }
          
          # Filter df_lib based on only the uni modifications with specified charge state found in OSW results
          df_lib %>%
            dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) -> df_lib
          
          cat('Peptidoforms of Charge State ', Isoform_Target_Charge, ' to Analyze..\n', paste(uni_mod, collapse = '\n'), '\n', sep='')
          if ( length(uni_mod)==1 & N_sample!=1 ){ cat(crayon::red('There is only 1 form... Skipping...\n')); skipped_bool=TRUE; next }
          
          
          # Display other peak group rank features
          if ( show_all_pkgrprnk==T ){
            osw_df_all <- getOSWData_( in_osw, run_name, precursor_id='', peptide_id=pep, mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=F )
            
            osw_df_all %>%
              dplyr::filter( FullPeptideName %in% uni_mod) %>%
              dplyr::filter( Charge %in% Isoform_Target_Charge )-> osw_df_all_filtered
            
            RT_Table <- table( osw_df_all_filtered$RT )
            
            # RT_pkgrps <- as.numeric(names(RT_Table)[RT_Table==2])
            RT_pkgrps <- as.numeric(names(RT_Table))
            if ( length(RT_pkgrps)==0 ){
              cat(crayon::red('WARNING: There were no common RT pkgrps found, will plot all pkgrps...\n'))
              RT_pkgrps <- as.numeric(names(RT_Table))
            }
            rm(osw_df_all, osw_df_all_filtered, RT_Table)
            
          } else {
            RT_pkgrps <- NULL
          }
          
          ###############################
          ## Get Site Determining Ions ##
          ###############################
          # cat(paste('   ~ Getting site Determining Ions information for each peptidoform\n', sep=''))
          # if ( N_sample!=1 ){
          #   print(uni_mod)
          #   print(len_unmod) 
          #   print(N_sample)
          #   print('getPairSiteDeterminingIonInformation_')
          #   uni_mod_list <- getPairSiteDeterminingIonInformation_( uni_mod, len_unmod )
          # } else {
          #   print(uni_mod)
          #   print(len_unmod)
          #   print(N_sample)
          #   print('getSiteDeterminingIonInformation_')
          #   uni_mod_list <- getSiteDeterminingIonInformation_( uni_mod, len_unmod )
          # }
          # ## Check for errors in uni_mod_list
          # if ( suppressWarnings( uni_mod_list=='skip' ) ){ cat('skipping...\n'); next }
          uni_mod_list <- NULL
          
          if ( plotIntersectingDetecting==T & plotUniqueDetecting==T ){
            # df_lib %>%
            #   dplyr::filter( MODIFIED_SEQUENCE==uni_mod[1] ) %>%
            #   dplyr::filter( DETECTING==1 ) -> df1
            # 
            # df_lib %>%
            #   dplyr::filter( MODIFIED_SEQUENCE==uni_mod[2] ) %>%
            #   dplyr::filter( DETECTING==1 ) -> df2
            # 
            # intersecting_mz <- intersect(df1$PRODUCT_MZ, df2$PRODUCT_MZ)
            
            product_mz <- unlist( lapply(uni_mod, function( mod_sequence ){ 
              df_lib %>%
                dplyr::filter( MODIFIED_SEQUENCE==mod_sequence ) %>%
                dplyr::filter( DETECTING==1 ) %>%
                dplyr::select( PRODUCT_MZ ) %>%
                as.matrix() %>%
                as.numeric()
            }) )
            
            intersecting_mz <- product_mz[table(product_mz)>1]
          } else {
            intersecting_mz <- NULL
          }
          
          #######################################
          ## Get Top Transitions with good PEP ##
          #######################################
          if ( use_top_trans_pep==T ){
            tictoc::tic(paste('   ~ Getting Top Transitions that scored a low PEP', sep=''))
            Transition_Scores <- getTransitionScores_( in_osw, run_name, precursor_id='', peptide_id=pep )
            
            if( is.null(unique(unlist(Transition_Scores$transition_id))) ){ cat(crayon::red('These is no Transition ID information... Skipping... \n')); next  }
            
            top_trans_mod_list <- lapply(uni_mod, function(mod, df_lib, Transition_Scores){
              df_lib %>%
                dplyr::filter( MODIFIED_SEQUENCE == mod ) -> mod_df_lib
              
              mod_df_lib %>%
                dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
                dplyr::filter( DETECTING==0 ) -> mod_df_lib_filtered
              
              Transition_Scores %>%
                dplyr::filter( FullPeptideName == mod) %>%
                dplyr::filter( !is.nan(transition_id) ) %>%
                dplyr::filter( transition_id %in% mod_df_lib_filtered$TRANSITION_ID ) %>%
                distinct() %>%
                arrange( transition_pep ) -> mod_transition_scores
            }, df_lib, Transition_Scores)
            names(top_trans_mod_list) <- uni_mod
            tictoc::toc()
          } else {
            top_trans_mod_list <- NULL
          }
          cat('   ~ Starting Plotting Action\n', sep='')
          plot_list <- list()
          max_Int <- 0
          uni_mod <- sort(uni_mod)
          for( mod in uni_mod ){ # names(uni_mod_list) ){
            cat( crayon::green('   --- Peptidoform: ', mod), '\n', sep='')
            ###########################
            ##     PLOT PRECURSOR    ##
            ###########################
            if ( plotPrecursor==T ){
              
              g <- ggplot2::ggplot()
              g <- mstools::getXIC( graphic_obj = g, 
                                    df_lib = df_lib, 
                                    mod = mod, 
                                    Isoform_Target_Charge = Isoform_Target_Charge,
                                    chromatogram_file = in_sqMass,  
                                    transition_type = 'precursor', 
                                    uni_mod_list = uni_mod_list, 
                                    max_Int = max_Int, 
                                    in_osw=NULL, 
                                    smooth_chromatogram=smooth_chromatogram, 
                                    doFacetZoom=F, 
                                    top_trans_mod_list=NULL, 
                                    show_n_transitions=show_n_transitions )
              max_Int <- g$max_Int
              g <- g$graphic_obj
            } else {
              g <- ggplot2::ggplot()
            }
            
            #################################
            ##   DETECTING TRANSITIONS     ##
            #################################
            
            ## INTERSECTING
            if ( plotIntersectingDetecting==T ){
              g <- mstools::getXIC( graphic_obj = g, 
                                    df_lib = df_lib, 
                                    mod = mod, 
                                    Isoform_Target_Charge = Isoform_Target_Charge,
                                    chromatogram_file = in_sqMass, 
                                    transition_type = 'detecting', 
                                    uni_mod_list = uni_mod_list, 
                                    max_Int = max_Int, 
                                    in_osw=NULL, 
                                    smooth_chromatogram=smooth_chromatogram, 
                                    doFacetZoom=F, 
                                    top_trans_mod_list=NULL, 
                                    show_n_transitions=show_n_transitions )
              max_Int <- g$max_Int
              g <- g$graphic_obj
            }
            ## UNIQUE
            if ( plotUniqueDetecting==T ){
              g <- mstools::getXIC( graphic_obj = g, 
                                    df_lib = df_lib, 
                                    mod = mod, 
                                    Isoform_Target_Charge = Isoform_Target_Charge,
                                    chromatogram_file = in_sqMass, 
                                    transition_type='detecting_unique', 
                                    uni_mod_list = uni_mod_list, 
                                    max_Int = max_Int, 
                                    in_osw=NULL, 
                                    smooth_chromatogram=smooth_chromatogram, 
                                    doFacetZoom=F, 
                                    top_trans_mod_list=NULL, 
                                    show_n_transitions=show_n_transitions )
              max_Int <- g$max_Int
              g <- g$graphic_obj
            }
            ###################################
            ##    IDENTIFYING TRANSITIONS   ###
            ###################################
            if ( plotIdentifying==T ){
              g <- mstools:: getXIC( graphic_obj = g, 
                                    df_lib = df_lib, 
                                    mod = mod, 
                                    Isoform_Target_Charge = Isoform_Target_Charge,
                                    chromatogram_file = in_sqMass,  
                                    transition_type='identifying', 
                                    uni_mod_list = uni_mod_list, 
                                    max_Int = max_Int, 
                                    in_osw=NULL, 
                                    smooth_chromatogram=smooth_chromatogram, 
                                    doFacetZoom=F, 
                                    top_trans_mod_list=NULL, 
                                    show_n_transitions=show_n_transitions, 
                                    plotIdentifying.Unique=plotIdentifying.Unique, 
                                    plotIdentifying.Shared=plotIdentifying.Shared, 
                                    plotIdentifying.Against=plotIdentifying.Against )
              max_Int <- g$max_Int
              g <- g$graphic_obj
              
            } else {
              cat(crayon::red('-- Identifying Transitions were not found for: ', underline(mod)), '\n', sep='')
            }
            
            ###################################
            ##     ADD OSW RESULTS INFO     ###
            ###################################
            g <- mstools::getXIC( graphic_obj = g, 
                                  df_lib = df_lib, 
                                  mod = mod, 
                                  Isoform_Target_Charge = Isoform_Target_Charge,
                                  chromatogram_file = in_sqMass, 
                                  transition_type='none', 
                                  uni_mod_list = uni_mod_list, 
                                  max_Int = max_Int, 
                                  in_osw = in_osw, 
                                  doFacetZoom=doFacetZoom, 
                                  top_trans_mod_list=NULL, 
                                  RT_pkgrps=RT_pkgrps, 
                                  show_manual_annotation=show_manual_annotation, 
                                  FacetFcnCall=FacetFcnCall, 
                                  show_legend = show_legend  )
            max_Int <- g$max_Int
            g <- g$graphic_obj
            
            plot_list[[mod]] <- g
          } 
          graphics.off()
          final_g <- (gridExtra::arrangeGrob(grobs=plot_list, nrow=length(uni_mod)))
          final_g_list[[uni_isoform_group_list_idx]] <- final_g
          final_g_list <- plot_list # TODO: Check if this is right 
          # grid.draw(final_g)
          
        } # uni_isoform_group_list_idx
        
        return( final_g_list )
        # record_list[[record_idx]] <- final_g
        # record_idx <- record_idx + 1
      }, error=function(e){
        if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
          # we reached timeout, apply some alternative method or do something else
          cat(crayon::red('Reached CPU Timelimit for ', crayon::underline(pep), ' from run: '), crayon::underline(run), '\n', sep='')
        } else {
          # error not related to timeout
          cat(crayon::red('There was an issue trying to process ', crayon::underline(pep), ' from run: '), crayon::underline(run), '\n', sep='')
          stop(e)
        }
      }
      )
      
    }, mc.cores = parallel::detectCores()-3 ) ## mclapply
    # }# for loop
    # print(record_list)
    # leh <- 0
    # bleh = leh + 2
    if( length(record_list)>0){
      draw_chrom_error <- tryCatch({
        if( doPlot==T & length(record_list[[1]])!=0 ){
          # cat( 'Length of record_list: ', length(record_list), '\n', sep='' )
          # print( record_list 
          # print('stop')
          # print((record_list))
          ### RECORD
          store_record <- list()
          store_idx <- 1
          for ( grph_idx in seq(1,length(record_list),1) ){
            # cat('Graph_Idx: ', grph_idx, '\n', sep ='')
            if ( !is.null(record_list[[grph_idx]]) ){
              # cat( 'eval( !is.null(record_list[[grph_idx]]) ) -> ', !is.null(record_list[[grph_idx]]), '\n', sep ='')
              if ( length(record_list[[grph_idx]])>1 ){
                # cat( ' if eval( length(record_list[[grph_idx]])>1 ) -> ', length(record_list[[grph_idx]])>1, '\n', sep='')
                ## There could be multiple isoform groups for 1 peptide in 1 run
                for ( sub_record_idx in seq(1,length(record_list[[grph_idx]]),1) ){
                  # print('list draw')
                  grid::grid.draw(record_list[[grph_idx]][[sub_record_idx]])
                  store_record[[store_idx]] <- recordPlot()
                  store_idx <- store_idx + 1
                  graphics.off()
                }
              } else {
                # cat( 'else eval( length(record_list[[grph_idx]])>1 ) -> ', length(record_list[[grph_idx]])>1, '\n', sep='')
                # print('no list draw')
                grid::grid.draw(record_list[[grph_idx]][[1]])
                store_record[[store_idx]] <- recordPlot()
                store_idx <- store_idx + 1
                graphics.off()
              }
            }
          }
          
          ## Check to see if save directory exits
          dir.create(file.path(getwd(), store_plots_subdir), showWarnings = FALSE )
          
          message( sprintf("Plots will be saved in: %s", file.path(getwd(), store_plots_subdir)) )
          
          # write and Append plots to a pdf file
          datatime_log <- sub(':', 'm', sub(':', 'h', gsub(' ', '_', as.character(as.POSIXct(Sys.time())))))
          pdf(paste(getwd(),store_plots_subdir,pep,'_ChargeState_', Charge_State, '_Precursor_', plotPrecursor, '_Detecting_', plotIntersectingDetecting, '_IdenifyingUnique_', plotIdentifying.Unique, '_IdentifyingShared_', plotIdentifying.Shared, '_IdentifyingAgainst_', plotIdentifying.Against, '_', datatime_log, '.pdf',sep=''), onefile = T, width=120, height=120, paper='a4r', pointsize=0.5 )
          for (myplot in store_record){
            if (class(myplot)=='recordedplot'){
              replayPlot(myplot)
            }
            
          }
          graphics.off()
          # rm(record_list)
          
          if ( printPlot==T ){
            for (myplot in store_record){
              if (class(myplot)=='recordedplot'){
                replayPlot(myplot)
              }
              
            }
          }
        } #doPlot END BLOCK
      }, error=function(e){
        # error not related to timeout
        cat(crayon::red('There was an issue trying to draw ', crayon::underline(pep)), '\n', sep='')
      })
    }# length of record plot list chec
  }
  tictoc::toc()
  return( record_list[[1]][[1]] )
}

