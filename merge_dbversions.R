source('R/packages.R')

rootdir = rprojroot::find_root(has_dir('src'))
datdir = file.path(rootdir, 'data')
resdir = file.path(rootdir, 'results')

#Read in different versions of the database ------------------------------------
#Original dataset from Chandima
db1 <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, 
                     'GEFIS_test_data/Database of EF assessments_IWMI_202109.xlsx')),
  sheet = 'EF estimates')

setnames(db1, old = c('EFA Type', 'Ecological Class'),
         new = c("E-flow Method/Model Type", "Ecological Condition(s) in Present Day"))

#Second dataset given to and processed by Chandima
db2 <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, "GEFIS_test_data/Master Data Table_20211020_format.xlsx"),
    sheet = 'Data')
) 

#Master database at October 28th freeze
#Formatting of the database from the original format
db3 <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, "GEFIS_test_data/Master Data Table_20211031.xlsx"),
    sheet = 'Data')
) 

#Merge databases ------------------------------------------------------------------
# db1 and db2
db1_2_bind <- cbind(db1[-16, 'no'], db2[2:(nrow(db1)),]) %>% #Row 16 is not in db2
  rbind(db2[nrow(db1):nrow(db2)+1], fill=T) %>% #Add the rest of db2
  .[!(Country %in% c('South Africa (from Retha)' , #Clean up the rest of db2
                     'Other African countries (from Retha)',
                     'Country',
                     NA)) &
      !(Country == 'Brazil' & is.na(Scale)),]
  

#db1+db2 and db3
db_all <- cbind(db1_2_bind[,"no", with=F], db3[1:nrow(db1_2_bind),])  %>%
  rbind(db1[16, c('no', 'Country', 'Latitude', 'Longitude',  #Add row 16 from db1
                  'E-flow Method/Model Type', 'Ecological Condition(s) in Present Day'), 
            with=F],
        fill = T) %>%
  .[no == 16, `:=`(
    River = 'Snowy River',
    `E-flow Location Name/No.` = 'Jindabyne Dam',
    `E-flow Requirement as Volume (106m3 per year)` = 44375,
    `Comments` = 'Snowy River at Jindabyne Dam flow was 3.2MCM/day prior to Dam construction. Present flow is 0.025MCM/day',
    `Link/Attachment for E-flows Assessment Report` = 'www.dnr.nsw.gov.au/water/snowy_benchmarking.shtml',
    `Sources of Additional Information` = 'p.dissanayake@cgiar.org'
  )] %>%
  rbind(db3[(nrow(db1_2_bind)+1):nrow(db3),], fill=T)
  
  
#Check database for cleaning
names(db_all)
unique(db_all$Country)
unique(db_all$Scale)
unique(db_all$`E-flow Method/Model Type`)
unique(db_all$`E-flow Method/Model Name`)
unique(db_all$`System to Determine Ecological Condition`)
unique(db_all$`Ecological Condition(s) in Present Day`)

#Replace Habitat simulation. in EF method name with habitat modelling for EF type?

#Standardize "Original" coordinates
setnames(db_all, old=c('Latitude', 'Longitude'), 
         new=c('Latitude_original', 'Longitude_original'))
db_all[, `:=`(
  latitude_parzer = parse_lat(Latitude_original),
  longitude_parzer = parse_lon(Longitude_original))
]

check_lat <- db_all[, .(Latitude_original, latitude_parzer)]
check_lon <- db_all[, .(Longitude_original, longitude_parzer)]

fwrite(db_all[is.na(no),], 
       file.path(resdir, 'Master_20211031_parzered_notIWMI.csv'))
