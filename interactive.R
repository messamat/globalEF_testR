source('_targets.R')

tar_make()
tar_manifest()
tar_glimpse()
tar_visnetwork()
 #Sys.setenv(RENV_DOWNLOAD_METHOD = "libcurl")