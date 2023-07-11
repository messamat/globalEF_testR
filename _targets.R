source("R/packages.R")
source("R/functions.R")

rootdir = rprojroot::find_root(rprojroot::has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

path_resgdb = file.path(resdir, "processing_outputs.gdb")
path_riveratlas = file.path(resdir, "RiverATLAS_v10tab.csv") #Too heavy to have as a target
path_riveratlas_clzlabels = file.path(datdir, 'HydroATLAS', 'HydroATLAS_v10_Legends.xlsx')

#tar_option_set(packages = c("biglm", "tidyverse"))
list(
  tar_target(
    efp_riverjoin,
    sf::st_read(dsn = path_resgdb, layer = 'EFpoints_20230424_clean_riverjoin')
  ),
  tar_target(
    path_efp_metadata,
    file.path(datdir, 'GEFIS_test_data', 'Master Data Table_20230424.xlsx'),
    format = 'file'
  ),
  tar_target(
    path_efp_mod,
    file.path(resdir, 'EFpoints_20230424_clean_globalEF.csv'),
    format = 'file'
  ),
  tar_target(
    path_riveratlas_meta,
    file.path(datdir, "HydroATLAS_metadata_MLMv11.xlsx"),
    format = 'file'
  )
  ,
  tar_target(
    clz_labels,
    as.data.table(read_xlsx(path_riveratlas_clzlabels, sheet='clz_cl'))
  )
  ,
  
  tar_target(
    eftab,
     read_excel(path=path_efp_metadata, sheet='Data') %>%
      format_eftab %>%
      merge(efp_riverjoin,
            by=c("Point_db", "UID_Mathis")) %>%
      comp_derivedvar
  ),

  tar_target(
    riveratlas_varsdt,
    selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
                          in_eftab = eftab)
  ),

  tar_target(
    rivernetwork,
    rformat_network(in_predvars = riveratlas_varsdt,
                    inp_riveratlasmeta = path_riveratlas_meta,
                    inp_riveratlas = path_riveratlas
    )
  ),

  tar_target(
    envhist,
    layout_ggenvhist(in_rivernetwork = rivernetwork,
                     in_sitedt = eftab,
                     in_predvars = riveratlas_varsdt,
                     out_path = file.path(resdir, 'figures', 
                                          paste0('envhist_', 
                                                 format(Sys.Date(), "%Y%m%d"),
                                                 '.pdf')
                     ))
  ),
  
  tar_target(
    countrytab,
    country_summary(
      in_tab = eftab,
      out_path = file.path(resdir, 'figures', 
                           paste0('country_method_histo_', 
                                  format(Sys.Date(), "%Y%m%d"),
                                  '.pdf')
      )
    )
  ),

  tar_target(
    efmap, #To redo with bigger sizes for France
    map_ef(in_efp = eftab,
           out_path = file.path(resdir, 'figures', 
                                paste0('efmap_', 
                                       format(Sys.Date(), "%Y%m%d"),
                                       '.png')
           ))
  ),

  tar_target(
    efp_efmod_join,
    join_efp_to_efmod(eftab,
                      path_efp_mod)

  ),
  
  tar_target(
    general_stats,
    list(
      N_efsites = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #Number of unique sites with ef data
                        nrow(unique(.SD)), 
                        .SDcols=c("UID_Mathis", "Point_db")],
      N_efa = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #Number of unique ef assessments
                    length(unique(EFUID))],
      N_countries = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #By country
                          length(unique(Country))],
      N_countries = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #By country
                          nrow(unique(.SD[, c("UID_Mathis", "Point_db")]))/
                            eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #Number of unique sites with ef data
                                  nrow(unique(.SD)), 
                                  .SDcols=c("UID_Mathis", "Point_db")], 
                          by='Country'],
      N_bysource = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)),  #By source
                         .N,
                         by="Point_db"],
      N_bymethod = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)), #By method
                         .N,
                         by="efname_ref"],
      N_bymethod_type = eftab[(!is.na(efvol_ref) | !is.na(efper_ref)), #By method type
                              .N,
                              by="eftype_ref"],
      minmax_MAF = efp_efmod_join[(!is.na(efvol_ref) | !is.na(efper_ref)), #By method type
                                  range(mar_ref, na.rm=T)],
      unique_sites_chars = efp_efmod_join[(!is.na(efvol_ref) | !is.na(efper_ref)) & 
                                   !(duplicated(efp_efmod_join, 
                                                by=c("UID_Mathis", "Point_db"))) & 
                                     Country != 'France', 
                                 list(median_MAF = median(mar_ref, na.rm=T),
                                      median_DA = median(up_area_skm_15s),
                                      N_DA_u100skm = .SD[up_area_skm_15s<100, .N, by=Country])
      ]
    )
  ),  
  
  tar_target(
    hydrology_comparison,
    compare_hydrology(efp_efmod_join,
                      in_clz_labels = clz_labels)
  ),
  
  tar_target(
    EFestimate_comparison,
    compare_EFestimate(in_efp_efmod_join = efp_efmod_join,
                       in_eftab_ensemble = hydrology_comparison$eftab_ensemblemod)
  ),
  tar_target(
    export_hydro_comparison,
    command = {
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_obspred_discharge.png')),
             hydrology_comparison$plot_obspred_discharge,
             width=8,
             height=8)
      ;
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_APEarea_discharge.png')),
             hydrology_comparison$plot_APEarea_discharge,
             width=8,
             height=8)
      ;
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_obspred_downscaledqtot.png')),
             hydrology_comparison$plot_obspred_downscaledqtot,
             width=8,
             height=8)
      ;
      ggsave(
        file.path(resdir, 'figures', 'hydrology_comparison',
                  paste0('plot_APEarea_downscaledqtot.png')),
        hydrology_comparison$plot_APEarea_downscaledqtot,
        width=8,
        height=8);
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_APEaridity_downscaledqtot.png')),
             hydrology_comparison$plot_APEaridity_downscaledqtot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_APEurban_downscaledqtot.png')),
             hydrology_comparison$plot_APEurban_downscaledqtot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_APEcrop_downscaledqtot.png')),
             hydrology_comparison$plot_APEcrop_downscaledqtot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('plot_APEreservoir_downscaledqtot.png')),
             hydrology_comparison$plot_APEreservoir_downscaledqtot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('qstats_dis_DA_plot.png')),
             hydrology_comparison$qstats_dis_DA_plot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('qstats_qtot_all_plot.png')),
             hydrology_comparison$qstats_qtot_all_plot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('qstats_qtot_DA_plot.png')),
             hydrology_comparison$qstats_qtot_DA_plot,
             width=8,
             height=8
      );
      ggsave(file.path(resdir, 'figures', 'hydrology_comparison',
                       paste0('qstats_qtot_country_plot.png')),
             hydrology_comparison$qstats_qtot_country_plot,
             width=8,
             height=8
      )
    }
    
    
    
    
    # EFestimate_comparison$plot_efper_ecpresent_ref
    # EFestimate_comparison$efvol_stats_all
    # EFestimate_comparison$efvol_stats_country
    # EFestimate_comparison$efvol_stats_eftype
    # EFestimate_comparison$efvol_stats_ensemble_all
    # EFestimate_comparison$efvol_stats_ensemble_country
    # EFestimate_comparison$efvol_stats_ensemble_eftype
    # EFestimate_comparison$efvol_stats_ensemble_eftypecountry
    # EFestimate_comparison$EFcompare_ensemble_best_vol_scatter
    # EFestimate_comparison$EFcompare_ensemble_best_vol_country_scatter
    # EFestimate_comparison$EFcompare_ensemble_best_vol_eftype_scatter
    # EFestimate_comparison$statsplot_all_vol_smape
    # EFestimate_comparison$statsplot_country_pbiassmape
    # EFestimate_comparison$statsplot_eftype_pbiassmape
    # EFestimate_comparison$statsplot_eftypecountry_pbiassmape
    # EFestimate_comparison$efper_stats_all
    # EFestimate_comparison$efper_stats_ensemble_all
    # EFestimate_comparison$efper_stats_ensemble_country
    # EFestimate_comparison$efper_stats_ensemble_eftype
    # EFestimate_comparison$EFcompare_ensemble_best_per_scatter
    # EFestimate_comparison$EFcompare_ensemble_best_per_country_scatter
    # EFestimate_comparison$EFcompare_ensemble_best_per_eftype_scatter
    # EFestimate_comparison$EFMARdiff_ensemble_best_per_country
    # EFestimate_comparison$EFDAdiff_ensemble_best_per_country 
    # EFestimate_comparison$EFGAIdiff_ensemble_best_per_country
    # EFestimate_comparison$EFGAIdiff_ensemble_best_per_eftype
    # EFestimate_comparison$statsplot_all_per_smape 
    
  )
)
