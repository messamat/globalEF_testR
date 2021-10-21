#---------------------- Utility functions --------------------------------------




#---------------------- Workflow function --------------------------------------
#Clean up and format master table
readformat_eftab <- function(inp_eftab) {
  #Read tab
  eftab <- as.data.table(
    readxl::read_xlsx(path = inp_eftab,
                      sheet = 'Data')
  ) 
  
  #Change column names
  namedt <-  rbindlist(
    list(
      list('E-flow Method/Model Type', 'eftype_ref'), 
      list('E-flow Method/Model Name', 'efname_ref'), 
      list('System to Determine Ecological Condition', 'ecname_ref'), 
      list('Ecological Condition(s) in Present Day', 'ecpresent_ref'),
      list('Ecological Condition(s) Set as E-flow Objective for Future Management', 'ecfuture_ref'),
      list('Type of Hydrological Regime Used to Calculate E-flow', 'hydrotype_ref'),
      list('Natural/Naturalised Mean Annual Runoff at E-flow Location (106m3)', 'mar_ref'),
      list('E-flow Requirement as Volume (106m3 per year)', 'efvol_ref'),
      list('E-flow as % of Natural/Naturalised Mean Annual Runoff', 'efper_ref'),
      list('Evaluation of Strength of E-flow Evidence', 'efstrength'),
      list('Specific Confidence Rating Where Applied', 'efconfidence')
    )
  ) %>% setnames(c('old', 'new'))
  
  setnames(eftab, namedt$old, namedt$new)
  
  #Format cols
  charcols <- names(eftab)[eftab[,sapply(.SD, is.character)]]
  
  catcolstonormalize <- c('eftype_ref',
                          'hydrotype_ref',
                          'efstrength',
                          'efconfidence')
  
  unique(eftab$efvol_ref)
  
  #Clean it up
  eftab_format <- eftab %>%
    .[, (charcols) := lapply(.SD, function(x) str_trim(x, side = 'both')), 
      .SDcols = charcols] %>%
    .[, (catcolstonormalize) := lapply(.SD, function(x) str_to_sentence(x)), 
      .SDcols = catcolstonormalize] %>%
    .[!(Country %in% c(NA, 'SA/ Botswana/ Zimbabwe/ Mozambique', 
                       'Other African countries (from Retha)')),] %>% #Remove extraneous rows
    .[!(Latitude %in% c('Latitude', 'Coordinates')),] %>%
    .[eftype_ref == 'Intermdiate', eftype_ref := 'Intermediate'] %>%
    .[eftype_ref %in% c('Hydrological index'), eftype_ref := 'Hydrological'] %>%
    .[eftype_ref %in% c('Comprehensive', 'Comprehensive/holistic'),  eftype_ref := 'Holistic'] %>%
    .[eftype_ref %in% c('Rapid 1', 'Rapid 2', 'Rapid 3'), eftype_ref := 'Rapid'] %>%
    .[eftype_ref %in% c('Rapid 3- intermediate', 'Rapid/ intermediate'), eftype_ref := 'Rapid/intermediate'] %>% 
    .[efname_ref %in% c('Desktop Reserve Model', 'RDRM'), eftype_ref := 'DRM'] %>%
    .[grep('.*[(]not sure whether present day[)]', ecpresent_ref),
      ecpresent_ref := gsub('\\s[(]not sure whether present day[)]', '', ecpresent_ref)] %>%
    .[hydrotype_ref %in% c('Existing water withdrawal'), hydrotype_ref := 'Existing water withdrawals'] %>%
    .[hydrotype_ref %in% c('Reduced impact scenario'), hydrotype_ref := 'Reduced impact scenario'] %>%
    .[Country == 'Brasil', efvol_ref := gsub('[,]', '', efvol_ref)]
  
  #Convert mar_ref, efvol_ref, and efper_ref to numeric
  numcols <- c('mar_ref', 'efvol_ref', 'efper_ref')
  eftab_format[, (numcols) := lapply(.SD, as.numeric), .SDcols = numcols]
  
  #For Brazil (double check that all records are in this case):
  #Convert Mean long term flow (m3/s) to Natural/Naturalised Mean Annual Runoff at E-flow Location (106m3)
  #Format efvol_ref
  eftab_format[Country == 'Brasil', mar_ref := efvol_ref*2]

  #Create unique ID for each study and site ##################################################Update with new data
  eftab_format[, IDef := .I]
  eftab_format[, IDefsite := as.numeric(factor(do.call(paste0,.SD))), .SDcols=c('Latitude', 'Longitude')] #Update with POINT_X and POINT_Y
  
  return(eftab_format)
}

#Join points with Master Table - does not work. Not the same tables at all.
join_eftabp <- function(in_eftab, in_efp) {
  eftabp <- merge(in_eftab, in_efp, by=c('Latitude', 'Longitude'))
  
}

#Join and format GEFIS stats to master table
joineftab_gefis <- function(in_eftab, inp_gefistats, idcol) {
  #Read GEFIS stats
  gefistats <- fread(inp_gefistats) %>%
    .[, no := as.numeric(gsub('ws_', '', no))]
  jointab <- merge(in_eftab, gefistats, by=idcol)
  
  return(jointab)
}

#-------------------- Summary
#Simple categorical histograms and tables of variables
multhist <- function(in_tab, idcols, measurecols) {
  catdt <- in_tab[, measurecols, with=F] %>%
    .[, .SD, .SDcols = function(x) {is.character(x) | is.factor(x)}]
  
  numdt <- in_tab[, measurecols, with=F] %>%
    .[, .SD, .SDcols = is.numeric] 

  if (ncol(catdt) > 1) {
    catdt_melt <- melt(cbind(catdt, in_tab[, idcols, with=F]), 
                       idcols = idcols, measure.vars = names(catdt))
    
    p_cat <- ggplot(catdt_melt, aes(x=value)) + 
      geom_histogram(stat='count') +
      stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
      coord_flip(clip='off') +
      scale_x_discrete('') + 
      scale_y_continuous('Number of e-flow assessments') +
      facet_wrap(~variable, scales = 'free_y') +
      theme_classic() +
      theme(axis.title.y = element_text(vjust=-5),
            plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  }
  
  if (ncol(numdt) > 1) {
    numdt_melt <- melt(cbind(numdt, in_tab[, idcols, with=F]), 
                       idcols = idcols, measure.vars = names(numdt))
    
    p_num <- ggplot(numdt_melt, aes(x=value)) + 
      geom_histogram() +
      #stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
      coord_flip(clip='off') +
      scale_x_sqrt() + 
      scale_y_continuous('Number of e-flow assessments') +
      facet_wrap(~variable, scales = 'free_y') +
      theme_classic() +
      theme(axis.title.y = element_text(vjust=-5),
            plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  }

  
  return(list(cat=p_cat, num=p_num))  
}



efdb_hists <- function(in_eftab) {
  colstoplot <- c('eftype_ref',
                  'efname_ref',
                  'ecpresent_ref',
                  'ecfuture_ref',
                  'hydrotype_ref',
                  #'efstrength',
                  #'efconfidence'
                  'mar_ref',
                  'efvol_ref',
                  'efper_ref'
  )
  
  all_p <- multhist(in_tab = in_eftab, idcols = 'IDef', measurecols = colstoplot)
  
  return(all_p)
}

#Histogram of number of sites x country x EFA_type
country_summary <- function(in_tab) {
  in_tab[, countryN := .N, by=Country]
  
  country_hist <- ggplot(in_tab, aes(x=reorder(Country, countryN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
    coord_flip(clip='off') +
    scale_x_discrete('Country') + 
    scale_y_continuous('Number of unique sites') +
    scale_color_discrete('Type of e-flow assessment') +
    theme_classic() +
    theme(axis.title.y = element_text(vjust=-5),
          plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  
  return(country_hist)
}

#Histogram of number of sites x country x EFA_type
ecoregions_summary <- function(in_tab) {
  in_tab[, countryN := .N, by=Country]
  
  country_hist <- ggplot(in_tab, aes(x=reorder(Country, countryN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
    coord_flip(clip='off') +
    scale_x_discrete('Country') + 
    scale_y_continuous('Number of unique sites') +
    scale_color_discrete('Type of e-flow assessment') +
    theme_classic() +
    theme(axis.title.y = element_text(vjust=-5),
          plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  
  return(country_hist)
}


#Representativeness:
#Climate zones
#

#Match

#-- Plot studies
map_ef <- function(in_efp) {
  datap <- st_as_sf(x = in_efp[!is.na(POINT_X),],
                    coords = c('POINT_X', 'POINT_Y'),
                    crs = 4326)
  
  
#Intersect studies with climate zones
# 
# 
# #Intersect with climate
# tempfst <- file.path(resdir, 'datap_gens.fst')
# if (!file.exists(tempfst)) {
#   gens <- read_sf(gens_shp) %>%
#     st_transform(crs=4326)
#   datap_gens <- st_join(datap, gens) %>%
#     as.data.table %>%
#     cbind(st_coordinates(datap))
#   write.fst(datap_gens[,-"geometry", with=F], tempfst)
# } else {
#   datap_gens <- read.fst(tempfst) %>%
#     as.data.table
# }

  # cluster all points using a hierarchical clustering approach
  studyclus <- geodist(st_coordinates(datap),
                       measure='geodesic') %>%
    as.dist %>%
    hclust(method="complete")
  

  # %>%
  #   .[, Climate := factor(Climate, levels=.SD[order(GEnZname), unique(Climate)])]
  datap_gens[, .N, by=.(cluster)]
  
  clustcentro <- in_efp[, `:=`(cluster = cutree(studyclus, h=200000))] %>%#,
                                      # Climate = gsub('^[A-Z][.]\\s', '', GEnZname))
    .[, list(lon = mean(POINT_X), lat = mean(POINT_Y), N = .N), 
      by=.(cluster)] %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326)
  
  
  #Map it
  crs_wintri = "+proj=wintri +datum=WGS84 +no_defs +over"
  
  # gens_ras <- raster(gens_ras) %>%
  #   projectRaster(crs=crs_wintri)
  cluscentro_wintri <- st_transform(clustcentro, crs_wintri)
  
  wcountries <- rnaturalearth::ne_countries(
    scale = "medium", returnclass = "sf") %>%
    st_as_sf %>%
    st_transform(crs_wintri, use_gdal = FALSE)
  
  
  studiesmap <- ggplot() +
    geom_sf(data=wcountries, color='#bdbdbd', alpha=1/2) +
    geom_sf(data = cluscentro_wintri, 
            aes(size=N), alpha=0.75) + #color = alpha('#3182bd', 1/2)
    scale_size_continuous(name=str_wrap('Number of studies within 200 km', 20),
                          range=c(1, 5), breaks=c(1, 2, 5, 9)) +
    # scale_color_manual(name='', values=c('#08519c', '#2171b5', '#67a9cf', '#006d2c',
    #                                      '#41ab5d', '#fd8d3c', '#e31a1c',
    #                                      '#800026')) +
    coord_sf(datum = NA, expand=F) +
    theme_map() +
    theme(legend.position=c(0, 0.25),
          #legend.direction="horizontal",
          legend.text = element_text(size=9),
          legend.title = element_text(size=9),
          legend.key.height = unit(0.1,"in"),
          legend.background = element_blank())

}

#---------------------- Report functions ---------------------------------------
