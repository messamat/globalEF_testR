source('R/packages.R')

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

#Master database at October 28th freeze
#Formatting of the database from the original format
db3 <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, "GEFIS_test_data/Master Data Table_20211104.xlsx"),
    sheet = 'Data')
) 

#Standardize "Original" coordinates
setnames(db3, old=c('Latitude', 'Longitude'), 
         new=c('Latitude_original', 'Longitude_original'))
db3[, `:=`(
  latitude_parzer = parse_lat(Latitude_original),
  longitude_parzer = parse_lon(Longitude_original))
]

check_lat <- db3[, .(Latitude_original, latitude_parzer)]
check_lon <- db3[, .(Longitude_original, longitude_parzer)]

fwrite(db3[is.na(no),], 
       file.path(resdir, 'Master_20211104_parzered_notIWMI.csv'))
