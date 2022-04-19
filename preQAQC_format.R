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


#------------------------------- Format Mexico --------------------------------
#Original source of data for officially-determined eflows:
#https://www.gob.mx/conagua/documentos/programa-nacional-hidrico-pnh-2020-2024

outdir <- file.path(datdir, 'GEFIS_test_data', 'Data by country', 'Mexico')
download_mexico <- function(outdir) {
  onlinedat <- "https://www.mdpi.com/2071-1050/13/3/1240/s1"
  zip_path <- file.path(outdir, 'sustainability-1066121-supplementary.zip')
  if (!file.exists(zip_path)) {
    download.file(url = onlinedat, 
                  destfil = zip_path,
                  method = 'auto',
                  mode = 'wb')
  } else {
    print(paste(zip_path, "already exists... Skipping download."))
  }
  
  unzipped_list <- unzip(zip_path, exdir = dirname(zip_path))
  
  
  onlinebasins <- "http://201.116.60.46/DatosAbiertos/Cuencas_hidrolOgicas_que_cuentan_con_reserva.zip"
  zip_path_basins <- file.path(outdir, basename(onlinebasins))
  if (!(file.exists(zip_path_basins))) {
    download.file(url = onlinebasins, 
                  destfil = zip_path_basins,
                  method = 'auto',
                  mode = 'wb')
  } else {
    print(paste(zip_path_basins, "already exists... Skipping download."))
  }
  unzipped_list_basins <- unzip(zip_path_basins, exdir = dirname(zip_path_basins))
  
  return(list(
    eflows_list = file.path(dirname(unzipped_list)[[1]], 
                            'sustainability-1066121-supplementary.xlsx'),
    eflows_basins = grep('.*[.]shp$', unzipped_list_basins, value=T)
  ))
}
mexico_dat <- download_mexico(outdir)

#Formatting of the database from the original format
mexico_refdata_raw <- as.data.table(
  readxl::read_xlsx(
    path = mexico_dat[['eflows_list']],
    sheet = '%_EWR_met',
    trim_ws	= T,
    range = "A1:S288"
    )
)

mexico_refdata_formatted <- mexico_refdata_raw[!is.na(Method),] %>% #Remove lines without data
  .[,`:=`(
    Hydro_Region_name = enc2utf8(Hydro_Region_name), #Make sure that spanish names are well encoded
    `River basin` = enc2utf8(`River basin`), #Make sure that spanish names are well encoded
    Country = 'Mexico',
    Scale = 'Basin',
    River = NA,
    E_flows_Assessment_Report_Av_13 = NA, #To be completed with a standardized format
    Report_Content_as_Summary_Da_14 = "Summary",
    Link_Attachment_for_E_flows__15 = "https://doi.org/10.3390/su13031240",
    Raw_Data_Available_on_Reques_16 = "Yes",
    E_flow_Method_Model_Type = NA,
    eftype_ref = NA,
    efname_ref = NA,
    ecname_ref = "Ecological objective (Objetivos ambientales, NORMA MEXICANA NMX-AA-159-SCFI-2012)",
    Ecological_condition_comment = "See NMX-AA-159-SCFI-2012 8/118 and 19/118.",
    ecpresent_ref = NA,
    hydrotype_ref = "To be completed",
    mar_unit = "hm3/yr",
    National_Legislation_for_E_f_39 = "Yes",
    Name_s__of_Laws = "Ley de Aguas Nacionales",
    Supporting_Regulation_s__Exi_41 = "Yes",
    Name_s__of_Regulations = "NORMA MEXICANA NMX-AA-159-SCFI-2012 QUE ESTABLECE EL PROCEDIMIENTO PARA LA DETERMINACIÓN DEL CAUDAL ECOLÓGICO EN CUENCAS HIDROLÓGICAS",
    Sources_of_Additional_Inform_43 = "ssalinas@ecosur.mx"
  )] 


#Change column names
mexico_colnamedt <-  rbindlist(
  list(
    list('Hydro_Region_name', 'Basin'), 
    list('River basin', 'Sub_Basin'), 
    list('ID_R2018', 'E_flow_Location_Name_No_'),
    list('MAR_hm3/yr', 'mar_ref'), 
    list('Env_Obj_NMX2016', 'ecfuture_ref'),
    list('Ewr_Est_hm3/yr', 'efvol_ref'),
    list('EWR1_%MAR', 'efper_ref')
  )
) %>% setnames(c('old', 'new'))

setnames(mexico_refdata_formatted, mexico_colnamedt$old, mexico_colnamedt$new)

fwrite(mexico_refdata_formatted,
       file.path(resdir, 'mexico_refdata_preformatted.csv'), 
       bom = TRUE)


# E_flow_Location_Name_No_
# Latitude
# Longitude




