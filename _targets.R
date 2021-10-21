source("R/packages.R")
source("R/functions.R")

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

#tar_option_set(packages = c("biglm", "tidyverse"))
list(
  tar_target(
    inp_efp,
    file.path(resdir, "processing_outputs.gdb/EFpoints_cleanjoin")
  ),
  
  tar_target(
    inp_eftab,
    file.path(datdir, "GEFIS_test_data/Master Data Table_20211020_format.xlsx"),
    format = 'file'
  ),
  
  tar_target(
    inp_gefistats,
    file.path(resdir, "efsites_gefisstats.csv"),
    format = 'file'
  ),
  
  tar_target(
    efp,
    as.data.table(
      sf::st_read(dsn = dirname(inp_efp),
                  layer = basename(inp_efp))
    )
  ),
  
  tar_target(
    eftab,
    readformat_eftab(inp_eftab)
  ),
  
  tar_target(
    eftab_gefis,
    joineftab_gefis(in_eftab = efp,  ###################UPDATE
                    inp_gefistats = inp_gefistats,
                    idcol = 'no')
  )
)