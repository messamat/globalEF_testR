#download_france_rhone -> manual processing -> format_points.py -> process_france_rhone.R

source('R/packages.R')
paste3 <- function(...,sep=", ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

rootdir <- rprojroot::find_root(has_dir('src'))
resdir <- file.path(rootdir, 'results')
rhonedir <- file.path(rootdir, 'data', 'GEFIS_test_data', 'Data by Country', 'Europe', 'France', 'Rhone')

ecostate <- fread(file.path(rhonedir, "etat_stations_filtrees.csv"), 
                  sep=";", encoding = "UTF-8")
#ECO column is the ecological state of the masse d'eau following the Water Framework Directive system

#Import EF points for the Rhone with minimal attributes
efpoints <- sf::st_read(file.path(resdir, "france_preprocessing.gdb"),
                        layer = "Rhone_EFpoints_reportsraw") 

#Import sub-basin polygons
ssbv_FRD <- sf::st_read(file.path(rhonedir, "ssbv_FRD.shp"))

#Create new points based on original coordinates and project to Lambert 93
efcoors_wgs84 <- st_drop_geometry(efpoints)[, c('UID_Mathis', 
                                                "Longitude_original", 
                                                "Latitude_original")]

efcoors_lambert <- st_as_sf(
  x = efcoors_wgs84,                         
  coords = c("Longitude_original", "Latitude_original"),
  crs = 4326) %>%
  sf::st_transform(efpoints, crs=2154) %>%
  st_join(ssbv_FRD[,"LIB_SSBV"], join = st_intersects) #Get basin name for each point

#Get full attributes for Rhone EF assessments
efdata <- readxl::read_xlsx(
  tail(
    file.path(rhonedir, 
              list.files(path=rhonedir, 
                         pattern="France_Rhone_catchment_extracted_.*")),
    n=1),
  
) %>%
  merge(efcoors_lambert, by='UID_Mathis', all.x=T, all.y=F) #Merge to projected points

efdata[, c("Longitude_original_Lambert", 
           "Latitude_original_Lambert")] <- st_coordinates(efdata$geometry)

#Only keep the most recent ecological state assessment
ecostate_last <- ecostate[order(annee), .SD[.N,], by=numero_station]

#Merge EF points with all ecological state assessments on the Masse d'eau 
#Multiple stations per Masse d'eau
efdata_ecostate <- merge(efdata, 
                         ecostate_last[, c("code_MDO", "x_lambert93", 
                                           "y_lambert93", "annee", "ECO")],  
                         by.x="Code_masse_deau", by.y="code_MDO",
                         all.x=T, all.y=T) %>%
  setDT
efdata_ecostate[!is.na(Longitude_original_Lambert) & 
                  !is.na(Latitude_original_Lambert),
                efecodist := sqrt((Longitude_original_Lambert - x_lambert93)^2 +
                                    (Latitude_original_Lambert - y_lambert93)^2)]


#Only keep the most recent ecological state assessment
efdata_ecostate_nearest <- efdata_ecostate[
  order(efecodist), .SD[1,], by=UID_Mathis] %>%
  .[, -c("x_lambert93", "y_lambert93", "geometry", "efecodist",
         "Longitude_original_Lambert", "Latitude_original_Lambert"
         ), with=F] %>%
  merge(efcoors_wgs84, by='UID_Mathis', all.x=T, all.y=F) %>%
  setnames(old=c("LIB_SSBV", "annee", 'ECO'), new=c("Sub-Basin", "ecostate_annee", "ecostate"))

#Remove records with NAs for e-flows
efdata_noNA <- efdata_ecostate_nearest[
  , lapply(.SD, function(x) gsub('^NA$', NA, x))] %>% #Replace NA strings with actual NAs
  .[(!is.na(Débit_minimum_biologique_m3s_DMO) | 
       !is.na(Débit_biologique_valeur_haute_m3s_DBh) |
       !is.na(Débit_biologique_valeur_basse_m3s_DBb) |
       !is.na(Débit_de_survie_proposée_m3s_DBs)),]

#Convert column types to have numeric ones
efdata_noNA[, 
            `:=`(Débit_minimum_biologique_m3s_DMO = as.numeric(Débit_minimum_biologique_m3s_DMO),
                 Débit_biologique_valeur_haute_m3s_DBh = as.numeric(Débit_biologique_valeur_haute_m3s_DBh),
                 Débit_biologique_valeur_basse_m3s_DBb = as.numeric(Débit_biologique_valeur_basse_m3s_DBb),
                 Débit_de_survie_proposée_m3s_DBs = as.numeric(Débit_de_survie_proposée_m3s_DBs)
                 )]

#Select e-flow value to use
#If there is a débit minimum biologique, use that one
efdata_noNA[!is.na(Débit_minimum_biologique_m3s_DMO),
            `:=`(eflow_selected = Débit_minimum_biologique_m3s_DMO,
                 eflow_selected_type = "Débit minimum biologique",
                 eflow_selected_methodtype = Méthode_type_Débit_minimum_biologique,
                 eflow_selected_methodname = Méthode_nom_Débit_minimum_biologique
                 )]

#If no débit minimum biologique, but a single "débit biologique" value, use that one
efdata_noNA[is.na(Débit_minimum_biologique_m3s_DMO) & 
              ((!is.na(Débit_biologique_valeur_haute_m3s_DBh) &
                   is.na(Débit_biologique_valeur_basse_m3s_DBb)) | 
                  (is.na(Débit_biologique_valeur_haute_m3s_DBh) &
                     !is.na(Débit_biologique_valeur_basse_m3s_DBb))),
            `:=`(eflow_selected = rowMeans(.SD, na.rm=T),
                 eflow_selected_type = "Débit biologique - single value",
                 eflow_selected_methodtype = fifelse(
                   (Méthode_type_DBh == Méthode_type_DBb) & 
                     !(is.na(Méthode_type_DBh) | is.na(Méthode_type_DBb)),  #Use %in% to deal with NAs
                   as.character(Méthode_type_DBh),
                   paste3(Méthode_type_DBh, Méthode_type_DBb, sep="")
                 ),
                 eflow_selected_methodname = fifelse(
                   Méthode_nom_DBh == Méthode_nom_DBb &
                     !(is.na(Méthode_nom_DBh) | is.na(Méthode_nom_DBb)), 
                   as.character(Méthode_nom_DBh),
                   paste3(Méthode_nom_DBh, Méthode_nom_DBb, sep="")
                 )
            ),
            .SDcols = c('Débit_biologique_valeur_haute_m3s_DBh',
                        'Débit_biologique_valeur_basse_m3s_DBb')]

#If no débit minimum biologique, but a low and high "débit biologique" value, get the average
efdata_noNA[is.na(Débit_minimum_biologique_m3s_DMO) & 
              (!is.na(Débit_biologique_valeur_haute_m3s_DBh) &
                 !is.na(Débit_biologique_valeur_basse_m3s_DBb)),
            `:=`(eflow_selected = rowMeans(.SD, na.rm=T),
                 eflow_selected_type = "Débit biologique - average of high and low values",
                 eflow_selected_methodtype = fifelse(
                   Méthode_type_DBh == Méthode_type_DBb & 
                     !(is.na(Méthode_type_DBh) | is.na(Méthode_type_DBb)), 
                   as.character(Méthode_type_DBh),
                   paste3(Méthode_type_DBh, Méthode_type_DBb, sep=' + ')
                 ),
                 eflow_selected_methodname = fifelse(
                   Méthode_nom_DBh == Méthode_nom_DBb &
                     !(is.na(Méthode_nom_DBh) | is.na(Méthode_nom_DBb)), 
                   as.character(Méthode_nom_DBh),
                   paste3(Méthode_nom_DBh, Méthode_nom_DBb, sep=' + ')
                 )
            ),
            .SDcols = c('Débit_biologique_valeur_haute_m3s_DBh',
                        'Débit_biologique_valeur_basse_m3s_DBb')]

#Don't use débit biologique de survie, it doesn't translate well to yearly or monthly e-flows

#Translate ecological status
#The notations mean the following:
#TBE: Très bon état - High
#BE: Bon état - Good
#MOY: Etat moyen - Moderate 
#MED: Etat médiocre - Poor
#MAUV: Etat mauvais - Bad
#IND: État indéterminé - NA/not determined
efdata_noNA[,
            ecostate := fcase(
              ecostate == 'TBE', 'High',
              ecostate == 'BE', 'Good',
              ecostate == 'MOY', 'Moderate',
              ecostate == 'MED', 'Poor',
              ecostate == 'MAUV', 'Bad',
              ecostate == 'IND', 'Not determined',
              default = NA
            )]


#Write records with eflow values
fwrite(efdata_noNA[!is.na(eflow_selected),], 
       file.path(rhonedir, 
                 paste0("France_Rhone_catchment_preformatted_", 
                        format(Sys.Date(), "%Y%m%d"),
                        ".csv")
       )
)

