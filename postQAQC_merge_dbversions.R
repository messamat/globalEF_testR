source('R/packages.R')

rootdir <- rprojroot::find_root(has_dir('src'))
datdir <- file.path(rootdir, 'data')
resdir <- file.path(rootdir, 'results')
process_gdb <- file.path(resdir, 'processing_outputs.gdb') 


#Read cleaned spatial data
EFp_clean <- as.data.table(st_read(dsn = process_gdb, 
                                   layer = 'EFpoints_20211104_clean')) %>%
  setnames('Id', 'ID')

#Original dataset from Chandima ------------------------------------
db1 <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, 
                     'GEFIS_test_data/Database of EF assessments_IWMI_202109.xlsx')),
  sheet = 'EF estimates')

setnames(db1, old = c('EFA Type', 'Ecological Class'),
         new = c("E-flow Method/Model Type", "Ecological Condition(s) in Present Day"))

db1_p <- merge(db1,
               merge(db1, EFp_clean, by='no', all.y=F)[
                 , c('River/Site Name', 'Latitude', 'Longitude', 
                     'Point_shift_mathis', 
                     'Comment_mathis', 'POINT_X', 'POINT_Y')],
               by=c('River/Site Name', 'Latitude', 'Longitude')
)

#Second dataset given to and processed by Chandima
db2_india <- fread(file.path(datdir, "Formatted_data_Chandima_20211102", "India.csv"))
db2_lesotho <- fread(file.path(datdir, "Formatted_data_Chandima_20211102", "Lesotho.csv"))
db2_miscaf <- fread(file.path(datdir, "Formatted_data_Chandima_20211102", "SA_Botswana_Zimbabwe_Mozambique.csv"))
db2_sa <- fread(file.path(datdir, "Formatted_data_Chandima_20211102", "Soth_africa_subset_3.csv"))

db2 <- rbindlist(list(db2_india, db2_lesotho, db2_miscaf, db2_sa), fill=T)

db2_p <- merge(db2, EFp_clean, by=c('ID', 'Country')) %>%
  setnames('River_', 'River')


#dataset from Oct 28th 2021 freeze
dbmaster <- as.data.table(
  readxl::read_xlsx(
    path = file.path(datdir, 
                     'GEFIS_test_data/Master Data Table_20211104.xlsx')),
  sheet = 'Data') %>%
  setnames("E-flow Location Name/No.", "E_flow_Location_Name_No.")


####################### Merge them all to master database ######################

db1_masterbind <- cbind(dbmaster[1:(nrow(db1)-1),][EFUID != 25,],
                        db1_p[order(no),
                              c('Point_shift_mathis', 'Comment_mathis',
                                'POINT_X', 'POINT_Y'), with=F][-16,] 
)


db2_masterbind <- merge(dbmaster, 
                        db2_p[, c("E_flow_Location_Name_No.", "River",
                                  'Point_shift_mathis', 
                                  'Comment_mathis', 'POINT_X', 'POINT_Y'), with=F],
                        by = c("E_flow_Location_Name_No.", "River"))


db3_masterbind <- merge(dbmaster, 
                        EFp_clean[,c('Point_shift_mathis', 'Comment_mathis', 
                                     'EFUID', 'POINT_X', 'POINT_Y')],
                        by='EFUID') 

all_masterbind <- rbind(db1_masterbind, db2_masterbind, db3_masterbind)

fwrite(all_masterbind, file.path(resdir, 'Master_20211104_QAQCed.csv'))
