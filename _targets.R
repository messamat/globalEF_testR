source("R/packages.R")
source("R/functions.R")

rootdir = rprojroot::find_root(rprojroot::has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

path_resgdb = file.path(resdir, "processing_outputs.gdb")
path_riveratlas = file.path(resdir, "RiverATLAS_v10tab.csv") #Too heavy to have as a target
# path_grdcp = file.path(resdir, 'GRDCstations_predbasic800.gpkg', 'GRDCstations_predbasic800')
# path_GRDCgaugedir = file.path(datdir, 'GRDCdat_day')

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
    eftab,
     read_excel(path=path_efp_metadata, sheet='Data') %>%
      format_eftab %>%
      merge(efp_riverjoin,
            by=c("Point_db", "UID_Mathis")) %>%
      comp_derivedvar
  )
,

  # tar_target(
  #   riveratlas_varsdt,
  #   selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
  #                         in_eftab = eftab)
  # ),
  # 
  # tar_target(
  #   rivernetwork,
  #   rformat_network(in_predvars = riveratlas_varsdt,
  #                   inp_riveratlasmeta = path_riveratlas_meta,
  #                   inp_riveratlas = path_riveratlas
  #   )
  # ),
  # 
  # tar_target(
  #   envhist,
  #   layout_ggenvhist(in_rivernetwork = rivernetwork,
  #                    in_sitedt = eftab,
  #                    in_predvars = riveratlas_varsdt)
  # ),
  # 
  # tar_target(
  #   countrytab,
  #   country_summary(
  #     in_tab = eftab
  #   )
  # ),
  # 
  # tar_target(
  #   efmap, #To redo with bigger sizes for France
  #   map_ef(in_efp = eftab)
  # ),

  tar_target(
    efp_efmod_join,
    join_efp_to_efmod(eftab,
                      path_efp_mod)

  )
  #,
  # 
  # tar_target(
  #   eftab_gefis_grdc,
  #   link_EFtoGRDC(in_eftab_gefis = eftab,
  #                 inp_riveratlas = path_riveratlas,
  #                 inp_grdcp = path_grdcp)
  # ),
  # 
  # tar_target(
  #   efgrdc_compare,
  #   compare_EFtoGRDC(eftab_gefis_grdc)
  # ),
  # 
  # tar_target(
  #   hydrology_comparison,
  #   compare_hydrology(efp_efmod_join)
  # ),
  # 
  # tar_target(
  #   EMC_comparison,
  #   compare_EMC(eftab,
  #               in_riveratlas_varsdt = riveratlas_varsdt)
  # ),
  # 
  # tar_target(
  #   EFestimate_comparison,
  #   compare_EFestimate(eftab)
  # ),
  # 
  # tar_target(
  #   mask_analysis,
  #   check_masking(eftab)
  # )
)
