source("R/packages.R")
source("R/functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

path_riveratlas = file.path(resdir, "RiverATLAS_v10tab.csv") #Too heavy to have as a target
path_grdcp = file.path(resdir, 'GRDCstations_predbasic800.gpkg', 'GRDCstations_predbasic800')
path_GRDCgaugedir = file.path(datdir, 'GRDCdat_day')

#tar_option_set(packages = c("biglm", "tidyverse"))
list(
  tar_target(
    path_efp,
    file.path(resdir, "Master_20211104_QAQCed_riverjoin.csv"),
    format = 'file'
  ),
  
  tar_target(
    path_riveratlas_meta,
    file.path(datdir, "HydroATLAS_metadata_MLMv11.xlsx"),
    format = 'file'
  ),
  
  tar_target(
    eftab,
    fread(path_efp, encoding = 'UTF-8') %>%
      comp_derivedvar %>%
      format_eftab
  )
  ,

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
                     in_predvars = riveratlas_varsdt)
  ),

  tar_target(
    countrytab,
    country_summary(
      in_tab = eftab
    )
  ),

  tar_target(
    efmap,
    map_ef(in_efp = eftab)
  ),

  tar_target(
    eftab_gefis_grdc,
    link_EFtoGRDC(in_eftab_gefis = eftab,
                  inp_riveratlas = path_riveratlas,
                  inp_grdcp = path_grdcp)
  ),

  tar_target(
    efgrdc_compare,
    compare_EFtoGRDC(eftab_gefis_grdc)
  ),
  
  tar_target(
    hydrology_comparison,
    compare_hydrology(eftab)
  ),
  
  tar_target(
    EMC_comparison,
    compare_EMC(eftab,
                in_riveratlas_varsdt = riveratlas_varsdt)
  ),
  
  tar_target(
    EFestimate_comparison,
    compare_EFestimate(eftab)
  ),
  
  tar_target(
    mask_analysis,
    check_masking(eftab)
  )
)
