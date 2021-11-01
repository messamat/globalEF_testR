source("R/packages.R")
source("R/functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

path_efp = file.path(resdir, "processing_outputs.gdb/EFpoints_cleanjoin")
path_grdcp = file.path(resdir, 'GRDCstations_predbasic800.gpkg', 'GRDCstations_predbasic800')
path_GRDCgaugedir = file.path(datdir, 'GRDCdat_day')

#tar_option_set(packages = c("biglm", "tidyverse"))
list(

  tar_target(
    path_riveratlas,
    file.path(resdir, "RiverATLAS_v10tab.csv"),
    format = 'file'
  ),
  
  tar_target(
    path_eftab,
    file.path(datdir, "GEFIS_test_data/Master Data Table_20211020_format.xlsx"),
    format = 'file'
  ),
  
  tar_target(
    path_gefistats,
    file.path(resdir, "efsites_gefisstats.csv"),
    format = 'file'
  ),
  
  tar_target(
    path_riveratlas_meta,
    file.path(datdir, "HydroATLAS_metadata_MLMv11.xlsx"),
    format = 'file'
  ),
  
  tar_target(
    efp,
    as.data.table(
      sf::st_read(dsn = dirname(path_efp),
                  layer = basename(path_efp))
    ) %>%
      comp_derivedvar
  ),
  
  tar_target(
    eftab,
    readformat_eftab(path_eftab)
  ),
  
  tar_target(
    eftab_gefis,
    joineftab_gefis(in_eftab = efp,  ###################UPDATE
                    inp_gefistats = path_gefistats,
                    idcol = 'no')
  ),
  
  tar_target(
    riveratlas_varsdt,
    selectformat_predvars(inp_riveratlas_meta = path_riveratlas_meta,
                          in_eftab = efp) ###################UPDATE
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
                     in_sitedt = efp, ###################UPDATE
                     in_predvars = riveratlas_varsdt)
  ),
  
  tar_target(
    eftab_gefis_grdc,
    link_EFtoGRDC(in_eftab_gefis = eftab_gefis,
                  inp_riveratlas = path_riveratlas,
                  inp_grdcp = path_grdcp)
    
  ),
  
  tar_target(
    grdc_paths,
    read_GRDCgauged_paths(inp_GRDCgaugedir = path_GRDCgaugedir,
                          in_gaugep = eftab_gefis_grdc,
                          exclude_missing = T)
  )
)