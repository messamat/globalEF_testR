#---------------------- Utility functions --------------------------------------




#---------------------- Workflow function --------------------------------------
#Clean up master table
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
  
  eftab_format <- eftab %>%
    .[!(Country %in% c(NA, 'SA/ Botswana/ Zimbabwe/ Mozambique', 
                       'Other African countries (from Retha)')),] %>% #Remove extraneous rows
    .[!(Latitude %in% c('Latitude', 'Coordinates')),] %>%
    .[, eftype_ref := tolower(eftype_ref)] %>%
    .[eftype_ref == 'intermdiate', eftype_ref := 'intermediate']
  
  return(eftab_format)
}

#Join points with Master Table - does not work. Not the same tables at all.
join_eftabp <- function(in_eftab, in_efp) {
  eftabp <- merge(in_eftab, in_efp, by=c('Latitude', 'Longitude'))
  
}


#-------------------- Summary
#Check river distance
#Histogram of number of sites x country x EFA_type
country_summary <- function(in_tab) {
  in_tab[, countryN := .N, by=Country]
  
  ggplot(in_tab, aes(x=reorder(Country, countryN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
    coord_flip(clip='off') +
    scale_x_discrete('Country') + 
    scale_y_continuous('Number of unique sites') +
    theme_classic() +
    theme(axis.title.y = element_text(vjust=-5),
          plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  
}

#Number of freshwater ecoregions
#


#Ecological conditions
#Total EFR
#Level_of_i

#Representativeness:
#Climate zones
#

#Match

#---------------------- Report functions ---------------------------------------
