#---------------------- Utility functions --------------------------------------
#Simple histograms and tables for multiple variables
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

#------ fread_cols -----------------
#' fread columns
#'
#' Fast data.table-based reading of a subset of columns from a table
#'
#' @param file_name path of the table to be read
#' @param cols_tokeep character vector, names of columns to read
#'
#' @return a data.table with the \code{cols_tokeep} that are found in the table
#'
#' @examples
#' fread_cols(iris, c('Sepal.Length', 'Sepal.Width')
#'
#' @export
fread_cols <- function(file_name, cols_tokeep) {
  #Only read the first row from the file
  header <- fread(file_name, nrows = 1, header = FALSE)
  #Check which columns are in the table
  keptcols <- cols_tokeep[cols_tokeep %chin% unlist(header)]
  missingcols <- cols_tokeep[!(cols_tokeep %chin% unlist(header))]
  paste('Importing', file_name, 'with ', length(keptcols),
        'columns out of ', length(cols_tokeep), 'supplied column names')
  #Reading in table
  paste(missingcols, 'columns are not in the file...')
  
  dt <- fread(input=file_name, header=TRUE,
              select=keptcols, verbose=TRUE)
  return(dt)
}
#------ comp_ymean -------------
#' Compute annual mean
#'
#' Compute annual mean of a variable based on columns of monthly means.
#'
#' @param in_dt data table, contains columns of monthly values.
#' @param fieldex character, example column name containing monthly value.
#' It must include the month number in %m format (two digits e.g. 02 for February)
#' @param mstart integer; position of the first digit of the month in the column
#' of monthly values.
#' @param outcol character; name of new column to be created containing monthly mean
#'
#' @return Modifies in_dt in place. Add column of weighted mean across months
#' (based on number of days in months)
#'
#' @export
comp_ymean<- function(in_dt, fieldex, mstart, outcol) {
  refd <- data.table(month=c('01', '02', '03', '04', '05', '06',
                             '07', '08', '09', '10', '11', '12'),
                     dimonth = c(31, 28.25, 31, 30, 31, 30,
                                 31, 31, 30, 31, 30, 31))
  refd[, mcols := paste0(substr(fieldex, 1, mstart-1),
                         month)]
  in_dt2 <- copy(in_dt[, refd$mcols, with=F])
  in_dt2[, comp_ymeanID := .I]
  in_dt[, outcol] <- melt(in_dt2,
                          id.vars='comp_ymeanID',
                          variable.name = 'mcols') %>%
    .[refd, on='mcols'] %>%
    .[, weighted.mean(value, dimonth), by=comp_ymeanID] %>%
    .[,'V1',with=F]
}
#------ selectformat_predvars -----------------
#' Select + format predictor variables
#'
#' Select candidate predictor variables to use in subsequent modelling and create
#' a table of formatted names.
#'
#' @param inp_riveratlas_meta path to metadata table for river atlas variables
#' @param in_eftab data.table of formatted gauging station summary statistics and hydro-environmental attributes.
#' 
#' 
#' @return data.table of selected predictor variable codes, names, category, 
#' attribute, sources, references, etc. Most sources are from RiverATLAS 
#' technical documentation available at https://www.hydrosheds.org/page/hydroatlas.
#' 
#' @export
selectformat_predvars <- function(inp_riveratlas_meta, in_eftab) {
  #---- List predictor variables ----
  monthlydischarge_preds <- paste0('DIS_',
                                   str_pad(seq(1, 12), 2, pad=0),
                                   '_CMS')
  
  predcols<- c(
    monthlydischarge_preds,
    'UPLAND_SKM',
    'dis_m3_pyr',
    'dis_m3_pmn',
    'dis_m3_pmx',
    'dis_m3_pvar',
    'dis_m3_pvaryr',
    'run_mm_cyr',
    'runc_ix_cyr', #runoff coefficient (runoff/precipitation)
    'sdis_ms_uyr', #specific discharge
    'sdis_ms_umn',
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    'dor_pc_pva', #anthropogenic - degree of regulation
    'gwt_m_cav',
    'ele_pc_rel',
    'slp_dg_cav',
    'slp_dg_uav',
    'clz_cl_cmj',
    'tmp_dc_uyr',
    'tmp_dc_cyr',
    'tmp_dc_cmn',
    'tmp_dc_cmx',
    'pre_mm_cyr',
    'pre_mm_uyr',
    'snw_pc_uyr',
    'snw_pc_cyr',
    'snw_pc_cmx',
    'glc_cl_cmj',
    'glc_pc_c16',
    'glc_pc_u16',
    
    
    'pnv_cl_cmj',
    'wet_pc_cg1',
    'wet_pc_cg2',
    'wet_pc_ug1',
    'wet_pc_ug2',
    'wet_pc_c07',
    'wet_pc_c09',
    'wet_pc_u07',
    'wet_pc_u09',
    'for_pc_use',
    'for_pc_cse',
    'ire_pc_use', #anthropogenic - irrigated area extent
    'ire_pc_cse', #anthropogenic - irrigated area extent
    'gla_pc_use',
    'gla_pc_cse',
    'prm_pc_use',
    'prm_pc_cse',
    'swc_pc_uyr',
    'swc_pc_cyr',
    'swc_pc_cmn',
    'lit_cl_cmj',
    'kar_pc_use',
    'kar_pc_cse',
    'ppd_pk_cav', #anthropogenic - pop density
    'ppd_pk_uav', #anthropogenic - pop density
    'urb_pc_cse', #anthropogenic - urban cover
    'urb_pc_use', #anthropogenic - urban cover
    'hft_ix_c93', #anthropogenic - human footprint
    'hft_ix_u93', #anthropogenic - human footprint
    'hft_ix_c09', #anthropogenic - human footprint
    'hft_ix_u09', #anthropogenic - human footprint
    'hdi_ix_cav', #anthropogenic - dev index
    'crp_pc_cse',
    'crp_pc_use',
    'pst_pc_cse',
    'pst_pc_use',
    'rdd_mk_cav',
    'rdd_mk_uav',
    
    'cly_pc_cav',
    'cly_pc_uav',
    'slt_pc_cav',
    'slt_pc_uav',
    'snd_pc_cav',
    'snd_pc_uav',
    # 'pre_mm_c01',
    # 'pre_mm_c02',
    # 'pre_mm_c03',
    # 'pre_mm_c04',
    # 'pre_mm_c05',
    # 'pre_mm_c06',
    # 'pre_mm_c07',
    # 'pre_mm_c08',
    # 'pre_mm_c09',
    # 'pre_mm_c10',
    # 'pre_mm_c11',
    # 'pre_mm_c12',
    # 'pre_mm_u01',
    # 'pre_mm_u02',
    # 'pre_mm_u03',
    # 'pre_mm_u04',
    # 'pre_mm_u05',
    # 'pre_mm_u06',
    # 'pre_mm_u07',
    # 'pre_mm_u08',
    # 'pre_mm_u09',
    # 'pre_mm_u10',
    # 'pre_mm_u11',
    # 'pre_mm_u12',
    'aet_mm_cyr',
    'aet_mm_uyr',
    'pet_mm_cyr',
    'pet_mm_uyr',
    # 'cmi_ix_cyr',
    # 'cmi_ix_uyr',
    'cmi_ix_cmn',
    'cmi_ix_umn',
    # 'cmi_ix_cvar',
    # 'cmi_ix_uvar',
    # 'cmi_ix_c01',
    # 'cmi_ix_c02',
    # 'cmi_ix_c03',
    # 'cmi_ix_c04',
    # 'cmi_ix_c05',
    # 'cmi_ix_c06',
    # 'cmi_ix_c07',
    # 'cmi_ix_c08',
    # 'cmi_ix_c09',
    # 'cmi_ix_c10',
    # 'cmi_ix_c11',
    # 'cmi_ix_c12',
    # 'cmi_ix_u01',
    # 'cmi_ix_u02',
    # 'cmi_ix_u03',
    # 'cmi_ix_u04',
    # 'cmi_ix_u05',
    # 'cmi_ix_u06',
    # 'cmi_ix_u07',
    # 'cmi_ix_u08',
    # 'cmi_ix_u09',
    # 'cmi_ix_u10',
    # 'cmi_ix_u11',
    # 'cmi_ix_u12',
    'ari_ix_cav',
    'ari_ix_uav',
    'wloss_pc_cav',
    'wdryp_pc_cav',
    'wwetp_pc_cav',
    'whfrq_pc_cav',
    'wseas_pc_cav',
    'wperm_pc_cav',
    'wfresh_pc_cav',
    'wloss_pc_uav',
    'wdryp_pc_uav',
    'wwetp_pc_uav',
    'whfrq_pc_uav',
    'wseas_pc_uav',
    'wperm_pc_uav',
    'wfresh_pc_uav',
    'bio1_dc_cav',
    'bio2_dc_cav',
    'bio3_dc_cav',
    'bio4_dc_cav',
    'bio5_dc_cav',
    'bio6_dc_cav',
    'bio7_dc_cav',
    'bio8_dc_cav',
    'bio9_dc_cav',
    'bio10_dc_cav',
    'bio11_dc_cav',
    'bio12_mm_cav',
    'bio13_mm_cav',
    'bio14_mm_cav',
    'bio15_mm_cav',
    'bio16_mm_cav',
    'bio17_mm_cav',
    'bio18_mm_cav',
    'bio19_mm_cav',
    'bio1_dc_uav',
    'bio2_dc_uav',
    'bio3_dc_uav',
    'bio4_dc_uav',
    'bio5_dc_uav',
    'bio6_dc_uav',
    'bio7_dc_uav',
    'bio8_dc_uav',
    'bio9_dc_uav',
    'bio10_dc_uav',
    'bio11_dc_uav',
    'bio12_mm_uav',
    'bio13_mm_uav',
    'bio14_mm_uav',
    'bio15_mm_uav',
    'bio16_mm_uav',
    'bio17_mm_uav',
    'bio18_mm_uav',
    'bio19_mm_uav'
  )
  
  #Check that all columns are in dt
  message(paste(length(predcols[!(predcols %in% names(in_eftab))]),
                'variables are missing from formatted gauge dataset'))
  
  #---- Associate HydroATLAS column names with variables names ----
  
  #Get predictor variable names
  metaall <- readxl::read_xlsx(inp_riveratlas_meta,
                               sheet='Overall') %>%
    setDT
  
  metascale <- readxl::read_xlsx(inp_riveratlas_meta,
                                 sheet='scale') %>%
    setDT %>%
    setnames(c('Key','Spatial representation'),
             c('Keyscale', 'Spatial.representation'))
  
  metastat <- readxl::read_xlsx(inp_riveratlas_meta,
                                sheet='stat') %>%
    setDT %>%
    setnames(c('Key','Temporal or statistical aggregation or other association'),
             c('Keystat', 'Temporal.or.statistical.aggregation.or.other.association'))
  
  meta_format <- as.data.table(expand.grid(`Column(s)`=metaall$`Column(s)`,
                                           Keyscale=metascale$Keyscale,
                                           Keystat=metastat$Keystat)) %>%
    .[metaall, on='Column(s)'] %>%
    .[metascale, on = 'Keyscale'] %>%
    .[metastat, on = 'Keystat',
      allow.cartesian=TRUE]
  
  meta_format[, `:=`(
    unit = substr(`Column(s)`, 5, 6),
    varcode = paste0(gsub('[-]{3}', '', `Column(s)`),
                     Keyscale,
                     fifelse(grepl("[0-9]",Keystat),
                             str_pad(Keystat, 2, side='left', pad='0'),
                             Keystat)),
    varname = paste(Attribute,
                    Spatial.representation,
                    Temporal.or.statistical.aggregation.or.other.association))]
  
  #Add newly generated variables to meta_format (variable labels)
  addedvars <- data.table(varname=c('Precipitation catchment Annual min/max',
                                    'Discharge watershed Annual min/max',
                                    'Discharge watershed Annual min/average',
                                    'Elevation catchment average - watershed average',
                                    'Runoff coefficient catchment Annual average',
                                    'Specific discharge watershed Annual average',
                                    'Specific discharge watershed Annual min',
                                    paste0('Discharge watershed ', month.name),
                                    'Drainage area',
                                    'Groundwater table depth catchment average'),
                          varcode=c('pre_mm_cvar',
                                    'dis_m3_pvar', 'dis_m3_pvaryr',
                                    'ele_pc_rel',
                                    'runc_ix_cyr',
                                    'sdis_ms_uyr', 'sdis_ms_umn',
                                    monthlydischarge_preds,
                                    'UPLAND_SKM',
                                    'gwt_m_cav'
                          )
  )
  
  oldcolnames <- c('Spatial.representation',
                   'Temporal.or.statistical.aggregation.or.other.association',
                   'Source Data')
  newcolnames <- c('Spatial representation',
                   'Temporal/Statistical aggreg.',
                   'Source')
  
  predcols_dt <- merge(data.table(varcode=predcols),
                       rbind(meta_format, addedvars, fill=T),
                       by='varcode', all.x=T, all.y=F)   %>%
    setnames(oldcolnames, newcolnames) %>%
    setorder(Category, Attribute,
             `Spatial representation`, `Temporal/Statistical aggreg.`)
  
  #Format table
  predcols_dt[varcode=='UPLAND_SKM', `:=`(
    Category = 'Physiography',
    Attribute= 'Drainage Area',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='',
    Source = 'HydroSHEDS',
    Citation = 'Lehner & Grill 2013'
  )]
  
  predcols_dt[varcode=='dis_m3_pvar', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='dis_m3_pvaryr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Natural Discharge',
    `Spatial representation`='p',
    `Temporal/Statistical aggreg.`='mn/yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='runc_ix_cyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Runoff coefficient',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2, WorldClim v2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='sdis_ms_uyr', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='yr',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  predcols_dt[varcode=='sdis_ms_umn', `:=`(
    Category = 'Hydrology',
    Attribute= 'Specific discharge',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn',
    Source = 'WaterGAP v2.2',
    Citation = 'Döll et al. 2003'
  )]
  
  # predcols_dt[grepl('DIS_[0-9]{2}.*', varcode), `:=`(
  #   Category = 'Hydrology',
  #   Attribute= 'Natural Discharge',
  #   `Spatial representation`='p',
  #   `Temporal/Statistical aggreg.`= gsub('[A-Z_]', '', varcode),
  #   Source = 'WaterGAP v2.2',
  #   Citation = 'Döll et al. 2003'
  # )]
  
  predcols_dt[varcode=='ele_pc_rel', `:=`(
    Category = 'Physiography',
    Attribute= 'Elevation',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='(cav-uav)/uav',
    Source = 'EarthEnv-DEM90',
    Citation = 'Robinson et al. 2014'
  )]
  
  predcols_dt[varcode=='cmi_ix_cvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index catchment monthly mn/mx'
  )]
  
  predcols_dt[varcode=='cmi_ix_uvar', `:=`(
    Category = 'Climate',
    Attribute= 'Climate Moisture Index',
    `Spatial representation`='u',
    `Temporal/Statistical aggreg.`='mn/mx',
    Source = 'WorldClim v2 & Global-PET v2',
    Citation = 'Fick et al. 2017',
    varname = 'Climate moisture index watershed monthly mn/mx'
  )]
  
  predcols_dt[varcode=='gwt_m_cav', `:=`(
    Category = 'Hydrology',
    Attribute= 'Groundwater table depth',
    `Spatial representation`='c',
    `Temporal/Statistical aggreg.`='av',
    Source = 'Global Groundwater Map',
    Citation = 'Fan et al. 2013'
  )]
  
  #Remove duplicates (which were created with keystat meaning different things e.g; 09 meaning september, 2009, class 9)
  predcols_dtnodupli<- predcols_dt[!(
    (Category == "Climate" &
       grepl('Class.*', `Temporal/Statistical aggreg.`)) |
      (Category == "Landcover" &
         !grepl('(Class|Spatial).*', `Temporal/Statistical aggreg.`))
  ),]  %>%
    .[grepl('hft_ix_[cu]09', varcode),
      `:=`(`Temporal/Statistical aggreg.`='2009',
           varname = gsub('(?<=Human\\sFootprint\\s(watershed|catchment)).*',
                          ' 2009',
                          varname,
                          perl=T)
      )] %>%
    .[grepl('hft_ix_[cu]93', varcode),
      `:=`(`Temporal/Statistical aggreg.`='1993',
           varname = gsub('(?<=Human\\sFootprint\\s(watershed|catchment)).*',
                          ' 1993',
                          varname,
                          perl=T)
      )] %>%
    unique(by='varcode')
  
  return(predcols_dtnodupli)
}




#------ comp_derivedvar -----------------
#' Compute derived variables
#'
#' Format and compute derived a set of environmental variables from input
#' dataset that contains
#' \href{https://www.hydrosheds.org/page/hydroatlas}{HydroATLAS} variables for
#' the purpose of predicting intermittency.
#'
#' @param dt data.table which contains rows as records and HydroATLAS
#'   environmental variables as columns.
#' @param copy (logical) whether to make a copy of \code{dt} (if TRUE) or to
#'   modify it in place (if FALSE)
#'
#' @details This function only formats and compute derived variables deemed
#'   relevant for the intermittency analysis. \cr
#'   \cr
#'   Steps include: \cr
#'   1. Convert some -9999 that should be 0 after checking them
#'   2. Count the number of -9999 values per column
#'   3. Convert the rest of the -9999 values to NA
#'   4. Compute a set of derived variables based on existing variables in
#'   HydroATLAS
#'   (e.g. \code{pre_mm_cvar= fifelse(pre_mm_cmx==0, 0, pre_mm_cmn/pre_mm_cmx)})
#'   5. Correct scaling from RiverATLAS
#'   (e.g. Degree of regulation is out of 10,000 in HydroATLAS)
#'
#' ```
#' @source Linke, S., Lehner, B., Dallaire, C. O., Ariwi, J., Grill, G., Anand,
#'   M., ... & Tan, F. (2019). Global hydro-environmental sub-basin and river
#'   reach characteristics at high spatial resolution. Scientific Data, 6(1),
#'   1-15.
#'
#' @export
comp_derivedvar <- function(in_dt, copy=FALSE) {
  if (copy) {
    in_dt2 <- copy(in_dt)
  } else {
    in_dt2 <- in_dt
  }
  
  
  #---- Inspect and correct -9999 and NA values ----
  print('Inspect and correct -9999 values')
  #check <- riveratlas[snw_pc_cyr == -9999,] #One reach in the middle of the Pacific
  in_dt2[snw_pc_cyr == -9999, snw_pc_cyr:=0]
  in_dt2[snw_pc_cmx == -9999, snw_pc_cmx:=0]
  
  #Places with very high slope value (in Greenland) have -9999 - replace with high value
  #in_dt2[, max(slp_dg_uav)]
  in_dt2[slp_dg_cav == -9999, `:=`(slp_dg_cav = 750,
                                   slp_dg_uav = 650)]
  
  print('Number of NA values per column')
  colNAs<- in_dt2[, lapply(.SD, function(x) sum(is.na(x)))]
  print(colNAs)
  
  print('Number of -9999 values per column')
  col9999<- in_dt2[, lapply(.SD, function(x) sum(x==-9999))]
  print(col9999)
  
  #-9999 in cly_pc_cav, slt, and snd are places with no soil mask (urban areas, lakes, glaciers, etc.)
  
  sgcols <- c('cly_pc_cav','slt_pc_cav', 'snd_pc_cav',
              'cly_pc_uav','slt_pc_uav', 'snd_pc_uav')
  
  #Convert -9999 to NAs
  for (j in which(sapply(in_dt2,is.numeric))) { #Iterate through numeric column indices
    if (!(j %in% which(names(in_dt2) %in% c(sgcols, 'wet_cl_cmj')))) {
      set(in_dt2,which(in_dt2[[j]]==-9999),j, NA) #Set those to NA if -9999
    }
  }
  
  #Scale variables based on HydroATLAS v1.0 documentation and v1.0.9 processing
  in_dt2[, `:=`(
    ari_ix_cav = ari_ix_cav/1000,
    ari_ix_uav = ari_ix_uav/1000,
    dor_pc_pva = dor_pc_pva/10,
    lka_pc_cse = lka_pc_cse/10,
    lka_pc_use = lka_pc_use/10,
    gwt_m_cav = gwt_cm_cav/100
  )]
  
  #---- Compute derived predictor variables ----
  print('Compute derived predictor variables')
  
  in_dt2[, swc_pc_cmn := do.call(pmin, c(.SD, list(na.rm=TRUE))),
         .SDcols= paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0))] %>% #Get minimum monthly swc
    .[, `:=`(
      #min/max monthly watershed discharge
      dis_m3_pvar=fifelse(dis_m3_pmx==0, 1, dis_m3_pmn/dis_m3_pmx),
      #min monthly/average yearly watershed discharge
      dis_m3_pvaryr=fifelse(dis_m3_pyr==0, 1, dis_m3_pmn/dis_m3_pyr),
      #catchment average elv - watershec average elev
      ele_pc_rel = fifelse(ele_mt_uav==0, 0, (ele_mt_cav-ele_mt_uav)/ele_mt_uav),
      #runoff coefficient (runoff/precipitation)
      runc_ix_cyr = fifelse(pre_mm_cyr==0, 0, run_mm_cyr/pre_mm_cyr),
      #Specific discharge
      sdis_ms_uyr = dis_m3_pyr/UPLAND_SKM,
      sdis_ms_umn = dis_m3_pmn/UPLAND_SKM
    )]
  return(in_dt2)
}
#------ formatscales ------------
#' Format plot scales
#'
#' Utility function to format plot scales for density distribution plots of environmental variables.
#' 
#' @param in_df data.frame with all records for environmental variables. 
#' Used to determine the appropriate range of values for each variable. In this case, 
#' a data.table of the river network hydro-environmental attributes.
#' @param varstopplot vector of variable names that will be plots and for which 
#' to return a list of scales
#'  
#' @return list of x and y scale objects + cartesian coordinates for ggplot
#' 
#' 
#' @export
formatscales <- function(in_df, varstoplot) {
  scales_x <- list(
    ari_ix_uav = scale_x_continuous(expand=c(0,0)),
    cly_pc_uav = scale_x_continuous(labels=percent_format(scale=1), expand=c(0,0)),
    cmi_ix_uyr = scale_x_continuous(),
    dis_m3_pyr = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               labels=c(0, 10^2,
                                        10^(0:log10(max(in_df$dis_m3_pyr)))),
                               expand=c(0,0)),
    dor_pc_pva = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    for_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    gla_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    kar_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    lka_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    sdis_ms_uyr = scale_x_continuous(expand=c(0,0)),
    snw_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    tmp_dc_uyr = scale_x_continuous(expand=c(0,0)),
    hdi_ix_cav = scale_x_continuous(expand=c(0,0)),
    hft_ix_c93 = scale_x_continuous(expand=c(0,0)),
    ORD_STRA = scale_x_continuous(expand=c(0,0)),
    UPLAND_SKM = scale_x_log10(breaks=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               labels=c(1, 10^2,
                                        10^(0:log10(max(in_df$UPLAND_SKM)))),
                               expand=c(0,0)),
    gwt_m_cav = scale_x_sqrt(expand=c(0,0)),
    ire_pc_use = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0))
  ) %>%
    .[(names(.) %in% names(in_df)) & names(.) %in% varstoplot]
  #Only keep those variables that are actually in df and that we want to plot
  
  scales_y <- unlist(rep(list(scale_y_continuous(expand=c(0,0))),
                         labels = scientific_format(),
                         length(scales_x)),
                     recursive=F) %>%
    setNames(names(scales_x))
  
  scales_y[['dis_m3_pmn']] <- scale_y_sqrt(expand=c(0,0))
  scales_y[['glc_pc_u16']] <- scale_y_continuous(trans='log1p',
                                                 breaks=c(10, 1000, 100000, 10000000))
  
  coordcart <- lapply(varstoplot, function(var) {
    coord_cartesian(xlim=as.data.table(in_df)[, c(min(get(var), na.rm=T),
                                                  max(get(var), na.rm=T))])
  }) %>%
    setNames(varstoplot)
  
  coordcart[['clz_cl_cmj']] <-  coord_cartesian(
    xlim=c(1,max(in_df$clz_cl_cmj)))
  coordcart[['kar_pc_use']] <-  coord_cartesian(
    xlim=c(0, 100))
  coordcart[['pet_mm_uyr']] <-  coord_cartesian(
    xlim=c(0, max(in_df$pet_mm_uyr)))
  coordcart[['ORD_STRA']] <-  coord_cartesian(
    xlim=c(1, 10))
  coordcart[['ari_ix_uav']] <-  coord_cartesian(
    xlim=c(0, 100))
  
  return(list(scales_x=scales_x, scales_y=scales_y, coordcart=coordcart))
}
#------ ggenvhist -------------
#' Plot of environmental histogram
#'
#' Utility function to create an individual density plot of the distribution of a 
#' given environmental variables across gauges and the whole global river network.
#' 
#' @param vartoplot (column) variable for which to produce a density plot.
#' @param in_sitedt data.table of gauging stations' environmental attributes.
#' @param in_rivdt data.table of global river network's environmental attributes 
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' @param scalesenvhist list of scale objects to format plot. From \link{formatscales}.
#' 
#' @return ggplot with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
ggenvhist <- function(vartoplot, in_sitedt, in_rivdt, in_predvars,
                      scalesenvhist) {
  print(vartoplot)
  
  varname <- in_predvars[varcode==vartoplot, Attribute]
  #paste0(Attribute, ' ',Keyscale,Keystat,' (',unit,')')]
  
  if (vartoplot == "clz_cl_cmj") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=as.factor(clz_cl_cmj)]
    gclz <- in_sitedt[,.N/in_sitedt[,.N],by=as.factor(clz_cl_cmj)]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))
    
    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'dodge', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
    
  } else if (vartoplot == "glc_pc_u16") {
    rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
                       by=glc_pc_u16]
    gclz <- in_sitedt[,.N/in_sitedt[,.N],by=glc_pc_u16]
    bindclz <- rbind(rivclz, gclz, idcol='source')%>%
      setnames(c( 'source', vartoplot, 'density'))
    
    penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
      geom_bar(aes(fill=as.factor(source)), stat='identity',
               position = 'identity', alpha=1/2, width=.6) +
      scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
    #
    #     penvhist <- ggplot(in_sitedt, aes_string(x=vartoplot)) +
    #       geom_histogram(data=in_rivdt, aes(weight = LENGTH_KM),
    #                      fill='#2b8cbe', alpha=0.5, bins=101) +
    #       geom_histogram(fill='#dd3497', alpha=0.5, bins=101)
    
  } else {
    penvhist <- ggplot(in_sitedt, aes_string(x=vartoplot)) +
      geom_density(data=in_rivdt, aes(weight = LENGTH_KM),
                   fill='#2b8cbe', alpha=0.5) +
      geom_density(fill='#dd3497', alpha=0.5) +
      ylab('Density')
  }
  
  penvhist <- penvhist +
    scalesenvhist$scales_x[[vartoplot]] +
    #scalesenvhist$scales_y[[vartoplot]] +
    scalesenvhist$coordcart[[vartoplot]] +
    xlab(varname) +
    theme_classic() +
    theme(strip.background=element_rect(colour="white", fill='lightgray'),
          legend.position = 'none',
          axis.title.y = element_blank(),
          axis.title = element_text(size=12))
  
  # if (which(vartoplot %in% varstoplot_hist)!=length(varstoplot_hist)) {
  #   penvhist <- penvhist +
  #     theme(legend.position='none')
  # }
  
  return(ggplotGrob(penvhist))
}

#------ getqstats ------------------------------------------------
getqstats <- function(dt, x, y, rstudthresh= 3, log=FALSE) {
  if (log) {
    in_form = paste0('log10(', y, '+0.1)~log10(', x, '+0.1)')
  } else {
    in_form = paste0(y,'~',x)
  }
  
  mod <- lm(as.formula(in_form), data=dt)
  
  outrows <- setDT(olsrr::ols_prep_rstudlev_data(mod)$`levrstud`)[
    abs(rstudent)>rstudthresh & color == 'outlier & leverage',]
  if (nrow(outrows) >0) {
    dtsub <- dt[-(outrows$obs),]
  } else {
    dtsub <- dt
  }
  
  mod_nooutliers <- lm(as.formula(in_form), data=dtsub)
  
  outstats <- dt[, list(pearsonr = round(cor(get(y), get(x)), 3),
                        mae = round(Metrics::mae(get(y), get(x)), 2),
                        smape = round(Metrics::smape(get(y), get(x)), 2),
                        pbias = round(Metrics::percent_bias(get(y), get(x)), 2),
                        rsq = round(summary(mod)$r.squared, 3),
                        rsq_nooutliers = round(summary(mod_nooutliers)$r.squared, 3),
                        n_total = .N,
                        noutliers = nrow(outrows)
  )
  ,]
  
  return(outstats)
}

#------ getcombistats -------------------------------------
#Run getqstats function on all combination of pairs in cols
getcombistats <- function(in_tab, cols, log) {
  compalist <- t(combn(cols, 2))
  combistats <- mapply(function(x, y) {
    getqstats(in_tab[!is.na(get(x)) & !is.na(get(y))],
              x, y, rstudthresh= 3, log=log)},
    x = compalist[,1], y = compalist[,2]
  ) %>%
    t %>%
    as.data.table %>%
    .[, comp := paste(compalist[,1], compalist[,2], sep='-')]
  return(combistats)
}

#------ grab_downstreamnet------------------------------------------
grab_downstreamnet <- function(in_net, in_startID, DAratiolim) {
  print(in_startID)
  in_reach <- in_net[HYRIV_ID == in_startID,]
  #Create placeholder vector with max number of records
  downIDs <- rep(NA, in_net[MAIN_RIV == in_reach$MAIN_RIV, .N])
  DAstart <- in_reach$UPLAND_SKM
  DAratio <- 1
  i = 1
  
  while (DAratio > DAratiolim) {
    downIDs[i] <- in_reach$HYRIV_ID
    in_reach <- in_net[HYRIV_ID == in_reach$NEXT_DOWN,]
    DAratio <- max(0, in_reach[, DAstart/UPLAND_SKM]) #In case reach river mouth and in_reach is an empty table
    i = i + 1
  }
  
  outdt <- data.table(
    startID = in_startID,
    downIDs = downIDs[!is.na(downIDs)]
  )
  
  return(outdt)
}
#------ read_GRDCgauged_paths -----------------
#' Read file paths to streamflow data from GRDC gauging stations
#'
#' Based on selection of gauges, create a list of paths to streamflow data
#' associated with gauges.
#'
#' @param inp_GRDCgaugedir path to directory containing streamflow data GRDC standard files.
#' @param in_gaugep table containing column named \code{GRDC_NO} with the
#' gauge IDs that will be used to generate file path.
#'
#' @return vector of paths to GRDC-formatted streamflow time series tables, assuming
#' that files are called "GRDC_NO.txt", GRDC_NO being replaced with a 7-digit integer.
#'
#' @export
read_GRDCgauged_paths <- function(inp_GRDCgaugedir, in_gaugep, exclude_missing=T) { #, gaugeid = 'GRDC_NO' down the line
  #Get data paths of daily records for gauge stations
  fileNames <- file.path(inp_GRDCgaugedir,
                         paste(
                           in_gaugep[!is.na(in_gaugep$GRDC_NO),]$GRDC_NO,
                           ".txt", sep=""))
  #Check whether any GRDC record does not exist
  existl <- do.call(rbind, lapply(fileNames, file.exists))
  print(paste(sum(!existl),
              'GRDC records do not exist...'))
  
  if (exclude_missing) {
    fileNames <- fileNames[existl]
  }
  return(fileNames)
}

#------ readformatGRDC -----------------
#' Read and pre-format GRDC data
#'
#' Reads text file of daily discharge data for a single GRDC station.
#' Creates columns for year, month, and date of last non-zero flow day +
#' computes yearly number of days of missing data
#'
#' @param path (character) path to the text file of daily discharge data in
#'   standard GRDC format.
#'
#' @return \link[data.table]{data.table} of daily discharge data with additional columns
#'
#' @export
readformatGRDC<- function(path, verbose = F) {
  if (verbose) {
    print(path)
  }
  
  #extract GRDC unique ID by formatting path
  gaugeno <- strsplit(basename(path), '[.]')[[1]][1]
  
  #Read GRDC text data
  gaugetab <- cbind(fread(path, header = T, skip = 40, sep=";",
                          colClasses = c('character', 'character', 'numeric',
                                         'numeric', 'integer')),
                    GRDC_NO = gaugeno)%>%
    setnames('YYYY-MM-DD', 'dates') %>%
    setorder(GRDC_NO, dates)
  
  #Format data
  gaugetab[, `:=`(year = as.numeric(substr(dates, 1, 4)), #Create year column
                  month = as.numeric(substr(dates, 6, 7)), #create month column
                  integervalue = fifelse(Original == round(Original), 1, 0) #Flag integer discharge values]
  )]
  
  #For each record, compute date of last non-zero flow day
  gaugetab[, prevflowdate := gaugetab[zero_lomf(Original),'dates', with=F]] %>% #Get previous date with non-zero flow
    .[Original != 0, prevflowdate:=NA] #If non-zero flow, set prevflowdate to NA
  
  #Compute number of missing days per year, excluding NoData values
  gaugetab[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | is.na(Original)),
           `:=`(missingdays = diny(year)-.N,
                datadays = .N),
           by= 'year']
  
  return(gaugetab)
}

#------ compute_ef_smakhtin ----------------------------
compute_ef_smakhtin <- function(in_dt) {
  stopifnot(is.data.table(in_dt))
  
  #Detect if data are daily or monthly
  if (in_dt[, mean(length(unique(dates))), by=year] > 12) {
    tstep <- 'daily'
  } else {
    tstep <- 'monthly'
  }
  
  #If daily, compute monthly statistics (removing years with maxgap)
  if (tstep == 'daily') {
    in_dt <- in_dt[!(Original %in% c(-999, -99, -9999, 99, 999, 9999) | 
                       is.na(Original)), 
                   mean(Original),
                   by=.(month, year)] %>%
      .[, dates_m := as.Date(paste(year, month, '01', sep='-'))]
  }
    
  #
  
  
}

compute_efindices <- function(path, maxgap) {
  stopifnot(is.character(path))
  
  in_dt <- readformatGRDC(path)
  
  #Subset dataset based on maxgap and compute number of years kept
  dtsub <- in_dt[missingdays <= maxgap,]
  nyears <- dtsub[, length(unique(year))]
  
  if (nyears > 0) {
    compute_ef_smakhtin(in_dt = dtsub)
  }

}

#---------------------- Workflow function --------------------------------------
#------ readformat_eftab -------------
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
  
  #unique(eftab$efvol_ref)
  
  #Clean it up ---------------------------- do typo corrections directly on the spreadsheet
  eftab_format <- eftab %>%
    .[, (charcols) := lapply(.SD, function(x) str_trim(x, side = 'both')), 
      .SDcols = charcols] %>%
    .[, (catcolstonormalize) := lapply(.SD, function(x) str_to_sentence(x)), 
      .SDcols = catcolstonormalize] %>%
    .[!(Country %in% c(NA, 'SA/ Botswana/ Zimbabwe/ Mozambique', 
                       'Other African countries (from Retha)')),] %>% #Remove extraneous rows
    .[!(Latitude %in% c('Latitude', 'Coordinates')),] %>%
    .[eftype_ref == 'Intermdiate', eftype_ref := 'Intermediate'] %>%
    .[efname_ref %in% c('Desktop Reserve Model', 'RDRM'), efname_ref := 'DRM'] %>%
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
  
  #Standardize EF type field
  eftab_format[eftype_ref %in% c('Hydrological index'), eftype_ref := 'Hydrological'] %>%
    .[eftype_ref %in% c('Comprehensive', 'Comprehensive/holistic'),  eftype_ref := 'Comprehensive/Holistic'] %>%
    .[eftype_ref %in% c('Rapid', 'Rapid 1', 'Rapid 2', 'Rapid 3'), eftype_ref := 'Hydrological time series analysis'] %>%
    .[eftype_ref %in% c('Rapid 3- intermediate', 'Rapid/ intermediate', 'Rapid/intermediate'), eftype_ref := 'Intermediate'] %>% 
    .[efname_ref %in% c('PROBFLO'), eftype_ref := 'Comprehensive/Holistic']
  unique(eftab_format$eftype_ref)  
  
  #Standardize EMC field
  unique(eftab_format$ecpresent_ref)
  eftab_format[ecpresent_ref %in% c('A', 'A/B', 'B', 'B/C', 'C', 'C/D', 'D'),
               ecpresent_reftype := 'Standard'] 
  eftab_format[ecpresent_ref %in% c('EWW', 'MaxW', 'MinW', 'MinW1', 'MinW2', 'RI'),
               ecpresent_reftype := 'Australia'] 
  eftab_format[ecpresent_ref %in% c('EWS', 'ET', 'Other', 'SR', 'MinD'),
               ecpresent_reftype := 'Other'] 
  
  rectdf <- data.table(
    xmin = rep(-Inf, 4),
    xmax = rep(Inf, 4),
    ymin = c(0, 0.25, 0.35, 0.5, 0.75),
    ymax = c(0.25, 0.35, 0.5, 0.75, 1),
    label = c('E', 'D', 'C', 'B', 'A')
  )
  
  #Create unique ID for each study and site ##################################################Update with new data
  eftab_format[, IDef := .I]
  eftab_format[, IDefsite := as.numeric(factor(do.call(paste0,.SD))), .SDcols=c('Latitude', 'Longitude')] #Update with POINT_X and POINT_Y
  
  return(eftab_format)
}

#------ join_eftabp -------------
#Join points with Master Table - does not work. Not the same tables at all.
join_eftabp <- function(in_eftab, in_efp) {
  eftabp <- merge(in_eftab, in_efp, by=c('Latitude', 'Longitude'))
  
}

#------ joineftab_gefis  -------------
#Join and format GEFIS stats to master table
joineftab_gefis <- function(in_eftab, inp_gefistats, idcol) {
  #Read GEFIS stats
  gefistats <- fread(inp_gefistats) %>%
    .[, no := as.numeric(gsub('ws_', '', no))]
  jointab <- merge(in_eftab, gefistats, by=idcol)
  
  # jointab[, `:=`(
  #   MAR_EF_Probable_perc =  MAR_EF_Probable_mcm/MAR_Natural_Annual_Runoff_v2,
  #   mar_nat_wg22_ls_year = dis_nat_wg22_ls_year*31.56/1000,
  #   mar_ant_wg22_ls_year = dis_ant_wg22_ls_year*31.56/1000,
  #   MAR_EF_A_perc = MAR_A_v2/MAR_Natural_Annual_Runoff_v2,
  #   MAR_EF_B_perc = MAR_B_v2/MAR_Natural_Annual_Runoff_v2,
  #   MAR_EF_C_perc = MAR_C_v2/MAR_Natural_Annual_Runoff_v2,
  #   MAR_EF_D_perc = MAR_D_v2/MAR_Natural_Annual_Runoff_v2
  # )] %>% 
  #   .[,
  #     MAR_EF_ecpresentref_perc := fcase( #Compute the GEFIS EF perc based on the EMC of the country estimate
  #       ecpresent_ref == 'A', MAR_EF_A_perc,
  #       ecpresent_ref == 'A/b', (MAR_EF_A_perc + MAR_EF_B_perc)/2,
  #       ecpresent_ref == 'B', MAR_EF_B_perc,
  #       ecpresent_ref == 'B/C', (MAR_EF_B_perc + MAR_EF_C_perc)/2,
  #       ecpresent_ref == 'C', MAR_EF_C_perc,
  #       ecpresent_ref == 'C/D', (MAR_EF_C_perc + MAR_EF_D_perc)/2,
  #       ecpresent_ref == 'D', MAR_EF_D_perc
  #     )  
  #   ]
  
  #Categorize GEFIS average EMC score
  # eftab_format[, 
  #              ecpresent_gefisws := fcase(
  #                EMC_10Variable_2 < 0.25, 'E',
  #                0.25 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.35, 'D',
  #                0.35 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.5, 'C',
  #                0.5 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.75, 'B',
  #                0.75 <= EMC_10Variable_2, 'A'
  #              )]
  
  return(jointab)
}

#------ rformat_network ------------------
#' Format river network
#'
#' Format river network, joining data, selecting variables, computing derived variables, etc.
#
#' @param in_predvars data.table of predictor variable codes, names and attributes. See \link{selectformat_predvars}.
#' @param in_monthlydischarge optional data.table of long-term mean monthly 
#' discharge values (from WaterGAP v2.2) for global river reaches.
#' @param inp_riveratlasmeta path to RiverATLAS attributes metadata table.
#' @param inp_riveratlas path to RiverATLAS (v1.0) attribute table (for the entire global river network).
#' @param inp_riveratlas2 path to new attributes calculated for global river network following RiverATLAS methodology
#' 
#' 
#' @return data.table of identifiers and hydro-environmental attribute for all reaches in the
#' RiverATLAS global river
#' 
#' @export
rformat_network <- function(in_predvars, in_monthlydischarge=NULL,
                            inp_riveratlasmeta, inp_riveratlas) {
  cols_toread <-  unique(
    c("HYRIV_ID", "HYBAS_L12", "HYBAS_ID03", "LENGTH_KM",'INLAKEPERC',
      'PFAF_ID05',
      in_predvars[, varcode], 'pop_ct_csu', 'pop_ct_usu',
      'ele_mt_cav','ele_mt_uav', 'gwt_cm_cav', 'dor_pc_pva', 'ORD_STRA',
      #paste0('pre_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
      #paste0('cmi_ix_c', str_pad(1:12, width=2, side='left', pad=0)),
      paste0('pet_mm_c', str_pad(1:12, width=2, side='left', pad=0)),
      paste0('swc_pc_c', str_pad(1:12, width=2, side='left', pad=0)))
  )
  
  riveratlas <- fread_cols(file_name=inp_riveratlas,
                           cols_tokeep = cols_toread)
  
  if (!is.null(in_monthlydischarge)) {
    riveratlas <- merge(riveratlas, as.data.table(in_monthlydischarge),
                        by.x='HYRIV_ID', by.y='REACH_ID')
  }
  
  setorder(riveratlas, HYRIV_ID)
  
  riveratlas_format <- comp_derivedvar(riveratlas)
  
  return(riveratlas_format)
}
#------ efdb_hists  -------------
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

#------ country_summary  -------------
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

#------ layout_ggenvhist --------------------------
#' Layout plots of environmental histograms
#'
#' Run plotting functions across predictor variables and arrange plots
#' 
#' @param in_rivernetwork data.table of global river network's environmental attributes. Here, output from \link{rformat_network}.
#' @param in_sitedt selected gauging stations' environmental attributes. Here, output from \link{write_gaugepreds}.
#' @param in_predvars data.table of predictor variable codes, names and attributes. 
#' Output from \link{selectformat_predvars}.
#' 
#' @details function used to produce Extended Data Fig. 8  c-p in Messager et al. 2021.
#' 
#' @return ggplots with two density distributions of the environmental variable,
#' one for the gauging stations and one for the global river network. 
#' 
#' 
#' @export
layout_ggenvhist <- function(in_rivernetwork, in_sitedt, in_predvars) {
  varstoplot_hist <- c(
    "dis_m3_pyr", "sdis_ms_uyr", "UPLAND_SKM",
    "ari_ix_uav", "gwt_m_cav", "tmp_dc_uyr",
    "clz_cl_cmj", "lka_pc_use", "kar_pc_use", 
    "for_pc_use", "glc_pc_u16", "snw_pc_uyr", 
    "ppd_pk_uav", "urb_pc_use", "ire_pc_use"
  )
  
  if ("dis_m3_pyr" %in% varstoplot_hist) {
    setDT(in_rivernetwork)[, dis_m3_pyr := dis_m3_pyr + 1]
    setDT(in_sitedt)[, dis_m3_pyr := dis_m3_pyr + 1]
  }
  
  #Get legend
  pleg <- ggplot(in_sitedt, aes(x=dis_m3_pyr)) +
    geom_density(alpha=1/2) +
    scale_fill_manual(values=c('#2b8cbe', '#dd3497'),
                      name = 'Dataset',
                      labels=c('Global river network',
                               'Training gauges')) +
    theme(text=element_text(size=14))
  
  # tmp <- ggplot_gtable(ggplot_build(pleg))
  # leg <- tmp$grobs[[
  #   which(sapply(tmp$grobs, function(y) y$name) == "guide-box")
  # ]]
  
  #Get scales
  scalesenvhist <- formatscales(in_df=in_rivernetwork, varstoplot=varstoplot_hist)
  
  #Plot each facet
  penvhist_grobs <- lapply(varstoplot_hist, ggenvhist,
                           in_sitedt = in_sitedt,
                           in_rivdt = in_rivernetwork,
                           in_predvars = in_predvars,
                           scalesenvhist = scalesenvhist)
  #Add legend
  #penvhist_grobs[[length(penvhist_grobs) + 1]] <- leg
  
  #Plot
  grid.newpage()
  do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=5))
}

#------ ecoregions_summary -------------
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


#------ map_ef -------------
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
  #studyclus[, .N, by=.(cluster)]
  
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



# in_eftab_gefis <- tar_read(eftab_gefis)
# inp_grdcp <- file.path(resdir, 'GRDCstations_predbasic800.gpkg', 'GRDCstations_predbasic800')
# 
# inp_riveratlas <- tar_read(path_riveratlas)

#------ link_EFtoGRDC ------------------------------------------
link_EFtoGRDC <- function(in_eftab_gefis, inp_riveratlas, inp_grdcp) {
  in_net <- fread_cols(inp_riveratlas, 
                       cols_tokeep = c("HYRIV_ID", "NEXT_DOWN", 
                                       "MAIN_RIV", "UPLAND_SKM"))
  
  grdcdt <- as.data.table(
    st_read(dsn = dirname(inp_grdcp),
            layer = basename(inp_grdcp))
  ) %>%
    .[!is.na(GRDC_NO), .(GRDC_NO, firstYear_kept_o1800, lastYear_kept_o1800,
                         totalYears_kept_o1800, HYRIV_ID)] %>%
    merge(in_net, by='HYRIV_ID')
  

  #reduce the rivernetwork to only reaches in the same watershed as the starts
  netsub <- in_net[MAIN_RIV %in% unique(grdcdt$MAIN_RIV),]
  remove(in_net)
  
  #iterate through each gauge
  grdc_downnet <- grdcdt[, grab_downstreamnet(in_net = netsub, 
                                              in_startID = HYRIV_ID, 
                                              DAratiolim = 0.9),
                         by = HYRIV_ID] 
  
  eftab_grdc_join <- merge(in_eftab_gefis, grdc_downnet,
                           by.x='HYRIV_ID', by.y = 'downIDs')
  
  return(eftab_grdc_join)
}




compare_EFtoGRDC <- function(in_eftab_grdc_join) {
  
  "mar_ref"
  "dis_m3_pyr"
  ""
  
}


#---------------------- TO FINISH IMPLEMENTING ---------------------------------------
#Change column names
allplots_temporaryfunc <- function(tabs) { 
  
  #--------------------------Format fields for comparison ---------------------
  tar_load(eftab_gefis)
  names(eftab_gefis)
  
  namedt_iwmi <-  rbindlist(
    list(
      list('EFA_Type', 'eftype_ref'), 
      list('EFA_Method', 'efname_ref'), 
      list('Ecological', 'ecpresent_ref'),
      list('Mean_Annua', 'mar_ref'),
      list('Total_EFR_', 'efper_ref')
    )
  ) %>% setnames(c('old', 'new'))
  
  setnames(eftab_gefis, namedt_iwmi$old, namedt_iwmi$new)
  
  eftab_gefis[, `:=`(
    MAR_EF_Probable_perc =  MAR_EF_Probable_mcm/MAR_Natural_Annual_Runoff_v2,
    mar_nat_wg22_ls_year = dis_nat_wg22_ls_year*31.56/1000,
    mar_ant_wg22_ls_year = dis_ant_wg22_ls_year*31.56/1000,
    MAR_EF_A_perc = MAR_A_v2/MAR_Natural_Annual_Runoff_v2,
    MAR_EF_B_perc = MAR_B_v2/MAR_Natural_Annual_Runoff_v2,
    MAR_EF_C_perc = MAR_C_v2/MAR_Natural_Annual_Runoff_v2,
    MAR_EF_D_perc = MAR_D_v2/MAR_Natural_Annual_Runoff_v2
  )] %>% 
    .[,
      MAR_EF_ecpresentref_perc := fcase( #Compute the GEFIS EF perc based on the EMC of the country estimate
        ecpresent_ref == 'A', MAR_EF_A_perc,
        ecpresent_ref == 'A/b', (MAR_EF_A_perc + MAR_EF_B_perc)/2,
        ecpresent_ref == 'B', MAR_EF_B_perc,
        ecpresent_ref == 'B/C', (MAR_EF_B_perc + MAR_EF_C_perc)/2,
        ecpresent_ref == 'C', MAR_EF_C_perc,
        ecpresent_ref == 'C/D', (MAR_EF_C_perc + MAR_EF_D_perc)/2,
        ecpresent_ref == 'D', MAR_EF_D_perc
      )  
    ]
  
  eftab_gefis[, 
               ecpresent_gefisws := fcase(
                 EMC_10Variable_2 < 0.25, 'E',
                 0.25 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.35, 'D',
                 0.35 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.5, 'C',
                 0.5 <= EMC_10Variable_2 & EMC_10Variable_2 < 0.75, 'B',
                 0.75 <= EMC_10Variable_2, 'A'
               )]
  
  eftab_gefis[ecpresent_ref %in% c('A', 'A/B', 'B', 'B/C', 'C', 'C/D', 'D'),
              ecpresent_reftype := 'Standard'] 
  eftab_gefis[ecpresent_ref %in% c('EWW', 'MaxW', 'MinW', 'MinW1', 'MinW2', 'RI'),
              ecpresent_reftype := 'Australia'] 
  eftab_gefis[ecpresent_ref %in% c('EWS', 'ET', 'Other', 'SR', 'MinD'),
              ecpresent_reftype := 'Other'] 
  
  unique(eftab_gefis$mar_ref)
  eftab_gefis[, `:=`(mar_ref = as.numeric(mar_ref),
                     efper_ref = as.numeric(efper_ref)/100
  )]
  
  #Create new fields for comparison
  eftab_gefis[, `:=`(
    mar_spe_gefiswg = 200*(MAR_Natural_Annual_Runoff_v2 - mar_nat_wg22_ls_year)/(
      MAR_Natural_Annual_Runoff_v2 + mar_nat_wg22_ls_year),
    mar_spe_gefisref = 200*(MAR_Natural_Annual_Runoff_v2 - mar_ref)/(
      MAR_Natural_Annual_Runoff_v2 + mar_ref),
    ef_spe_gefisref_probable = 200*(MAR_EF_Probable_perc - efper_ref)/(
      MAR_EF_Probable_perc + efper_ref),
    ef_spe_gefisref_emcref = 200*(MAR_EF_ecpresentref_perc - efper_ref)/(
      MAR_EF_ecpresentref_perc + efper_ref)
  )]
  
  #--------------------------Compare reference MAR, WaterGAP MAR, and GEFIS MAR-----
  qstats_clz <- eftab_gefis[,
                        getcombistats(in_tab = .SD,
                                      cols = c('MAR_Natural_Annual_Runoff_v2', 'mar_nat_wg22_ls_year', 'mar_ref'),
                                      log = T
                        ), by=clz_cl_cmj
  ]
  
  qstats <- eftab_gefis[,
                        getcombistats(in_tab = .SD,
                                      cols = c('MAR_Natural_Annual_Runoff_v2', 'mar_nat_wg22_ls_year', 'mar_ref'),
                                      log = T
                        )
  ]
  fwrite(qstats, file.path(resdir, 'qstats_all.csv'))
  
  hydrop_wggefis <- ggplot(eftab_gefis, aes(x=mar_nat_wg22_ls_year, y=MAR_Natural_Annual_Runoff_v2)) +
    geom_point() + 
    geom_abline() +
    scale_x_log10(name = 'WaterGAP MAR (mcm)') + 
    scale_y_log10(name = 'GEFIS MAR (mcm)') +
    theme_bw()
  
  hydrop_wgref <- ggplot(eftab_gefis, 
         aes(x=mar_nat_wg22_ls_year, y=as.numeric(mar_ref))) +
    geom_point() + 
    geom_abline() +
    scale_x_log10(name = 'WaterGAP MAR (mcm)') + 
    scale_y_log10(name = 'Local EFA MAR (mcm)')+
    theme_bw()
  
  hydrop_gefisref <- ggplot(eftab_gefis, 
                            aes(x=MAR_Natural_Annual_Runoff_v2, y=as.numeric(mar_ref))) +
    geom_point() + 
    geom_abline() +
    scale_x_log10(name = 'GEFIS MAR (mcm)') + 
    scale_y_log10(name = 'Local EFA MAR (mcm)')+
    theme_bw()
  

  spemissing_gefisref <- ggplot(eftab_gefis, aes(x=100*(1-MAR_boolean), y=mar_spe_gefisref)) +
    geom_point() +
    #geom_smooth(method='lm') +
    scale_x_continuous(name='Missing hydro. information (% area)') +
    scale_y_continuous(name='Percentage error GEFIS vs. local EFA')
  
  ((hydrop_wggefis + hydrop_gefisref)/(hydrop_wgref | spemissing_gefisref))
  
  
  
  #Show missing data
  ggplot(eftab_gefis, aes(x=as.factor(Country), y=1-MAR_boolean)) +
    geom_boxplot()
  ggplot(eftab_gefis, aes(x=str_wrap(as.factor(Country), 15), y=100*(1-EMC_boolean))) +
    geom_boxplot() +
    scale_y_continuous('Percentage watershed area with missing EF information') +
    scale_x_discrete('Country')
  
  #Compare reference EMC and GEFIS EMC
  unique(eftab_gefis$Ecological)
  eftab_gefis$ecpresent_ref
  eftab_gefis$EMC_10Variable_2
  
  #Prepare dt for plot background showing the difference classes
  rectdf <- data.table(
    xmin = rep(-Inf, 4),
    xmax = rep(Inf, 4),
    ymin = c(0, 0.25, 0.35, 0.5, 0.75),
    ymax = c(0.25, 0.35, 0.5, 0.75, 1),
    label = c('E', 'D', 'C', 'B', 'A')
  )
  
  #Test differences among ckasses
  compute_tukeyletters <- function(in_dt, ectype) {
    dtsub <- in_dt[ecpresent_reftype == ectype]
    aov_out<- aov(EMC_10Variable_2 ~ ecpresent_ref, data=dtsub)    
    tukey_out <- TukeyHSD(aov_out)
    cld <- multcompView::multcompLetters4(aov_out, tukey_out)
    as.data.table(as.data.frame.list(cld$ecpresent_ref),keep.rownames = T)
  }
  
  tukeyletters <- lapply(unique(eftab_gefis$ecpresent_reftype), function(ectype) {
    compute_tukeyletters(in_dt = eftab_gefis, ectype = ectype)
  }) %>%
    rbindlist(fill=T) %>%
    merge(unique(eftab_gefis[, .(ecpresent_ref, ecpresent_reftype)]),
          by.x = 'rn', by.y = 'ecpresent_ref')
  
  ggplot(eftab_gefis) + 
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax, 
                               ymin=ymin, ymax=ymax, fill = label),
              color='black') +
    geom_boxplot(aes(x=ecpresent_ref, y=EMC_10Variable_2)) +
    geom_text(data=tukeyletters, aes(x=rn, y=0.7, label=Letters)) +
    # geom_text(data=rectdf, aes(x='A', y=ymin, label=label),
    #           vjust = -1.5, hjust=2) +
    scale_fill_brewer(name = str_wrap('Average GEFIS ecological class', 20),
                      palette='RdYlBu', direction = -1) +
    scale_y_continuous('Average GEFIS threat index', limits=c(0,1), 
                       expand = c(0,0), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete('Present ecological class from e-flow assessment') +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14)) +
    facet_wrap(~ecpresent_reftype, scales = 'free_x')
  
  #Compute classification statistics + confusion matrix
  eftab_gefis[ecpresent_reftype == 'Standard', `:=`(
              ecpresent_ref_simplemin = fcase(
                ecpresent_ref == 'A/B', 'A',
                ecpresent_ref == 'B/C', 'B',
                ecpresent_ref == 'C/D', 'C',
                ecpresent_ref %in% c('A', 'B', 'C', 'D'), ecpresent_ref
              ),
              ecpresent_ref_simplemax = fcase(
                ecpresent_ref == 'A/B', 'B',
                ecpresent_ref == 'B/C', 'C',
                ecpresent_ref == 'C/D', 'D',
                ecpresent_ref %in% c('A', 'B', 'C', 'D'), ecpresent_ref
              ),
              ecpresent_ref_num = fcase(
                ecpresent_ref == 'A', 0.875,
                ecpresent_ref == 'A/B', 0.75,
                ecpresent_ref == 'B', 0.675,
                ecpresent_ref == 'B/C', 0.65,
                ecpresent_ref == 'C', 0.425,
                ecpresent_ref == 'C/D', 0.5,
                ecpresent_ref == 'D', 0.3
              )
  )] %>%
    .[, ecpresent_refgefisdiff := EMC_10Variable_2 - ecpresent_ref_num]
              
  eftab_gefis[ecpresent_reftype == 'Standard', 
              caret::confusionMatrix(as.factor(ecpresent_gefisws), as.factor(ecpresent_ref_simplemin))
  ]
  eftab_gefis[ecpresent_reftype == 'Standard', 
              caret::confusionMatrix(as.factor(ecpresent_gefisws), as.factor(ecpresent_ref_simplemax))
  ]
  
  #Check whether difference is due to missing data (counter-intuitive pattern!)
  ggplot(eftab_gefis, aes(x=EMC_boolean, y=abs(ecpresent_refgefisdiff))) + 
    geom_point() + 
    geom_quantile(quantiles=c(0.1, 0.5, 0.9))
  
  
  #Compare HydroATLAS predictors/variables to ref classes
  tar_load(riveratlas_varsdt)
  eftab_ecstressorsmelt <- melt(
    eftab_gefis[, .(dor_pc_pva, glc_pc_c16, glc_pc_u16, crp_pc_cse, crp_pc_use,
                    pst_pc_cse, pst_pc_use, ire_pc_cse, ire_pc_use, ppd_pk_cav,
                    ppd_pk_uav, urb_pc_cse, urb_pc_use, rdd_mk_cav, rdd_mk_uav,
                    ecpresent_ref, ecpresent_reftype, EMC_10Variable_2, no
    )],
    id.vars = c('no', 'ecpresent_ref', 'ecpresent_reftype', 'EMC_10Variable_2')
  ) %>%
    merge(riveratlas_varsdt, by.x = 'variable', by.y = 'varcode', all.x=T)
    
  
  ggplot(eftab_ecstressorsmelt, 
         aes(x=ecpresent_ref, y=value, fill=ecpresent_reftype)) +
    geom_boxplot() +
    facet_wrap(~varname, scales='free_y') + 
    theme_bw()
  
  ggplot(eftab_ecstressorsmelt, 
         aes(x=value, y=EMC_10Variable_2)) +
    geom_point() + 
    geom_smooth() +
    facet_wrap(~varname, scales='free') +
    theme_bw()
  
  #Compare reference EF to GEFIS EF - % and vol
  ggplot(eftab_gefis[!(efper_ref %in% c('28-Jun', 'Jun-48', '2-Jan', '18-24')),],
         aes(x=MAR_EF_Probable_perc, y=as.numeric(efper_ref))) +
    geom_point(aes(color=Country)) + 
    geom_abline() +
    #geom_quantile() +
    scale_x_continuous(name = 'GEFIS e-flow (% of MAR)', limits=c(0, 1)) +
    scale_y_continuous(name = 'Country e-flow (% of MAR)', limits=c(0, 1)) +
    coord_fixed(expand=c(0,0)) +
    facet_wrap(~eftype_ref) +
    theme_bw()
  
  ggplot(eftab_gefis[!(efper_ref %in% c('28-Jun', 'Jun-48', '2-Jan', '18-24')),],
         aes(x=MAR_EF_ecpresentref_perc, y=as.numeric(efper_ref))) +
    geom_point(aes(color=Country)) + 
    geom_abline() +
    geom_quantile() +
    scale_x_continuous(name = 'GEFIS e-flow (% of MAR)', limits=c(0, 1)) +
    scale_y_continuous(name = 'Country e-flow (% of MAR)', limits=c(0, 1)) +
    coord_fixed(expand=c(0,0)) +
    facet_wrap(~eftype_ref) +
    theme_bw()
  
  #Compute statistics for different EF types and countries: make table with errors, bias, etc.
  efstats_country <- eftab_gefis[!(is.na(efper_ref) & !is.na(MAR_EF_Probable_perc)) &
                                   EMC_boolean > 0.5,
                                 getqstats(dt = .SD,
                                           x = 'MAR_EF_Probable_perc',
                                           y = 'efper_ref',
                                           log = F
                                 ), by=.(Country)
  ]
  
  efstats_eftype <- eftab_gefis[!(is.na(efper_ref) & !is.na(MAR_EF_Probable_perc)) &
                                  EMC_boolean > 0.5,
                                getqstats(dt = .SD,
                                          x = 'MAR_EF_Probable_perc',
                                          y = 'efper_ref',
                                          log = F
                                ), by=.(eftype_ref)
  ]
  
  efstats_eftype_emcref <- eftab_gefis[!(is.na(efper_ref)) & !(is.na(MAR_EF_ecpresentref_perc)) &
                                         (EMC_boolean > 0.5) & (eftype_ref %in% c('CH','HA')),
                                       getqstats(dt = .SD,
                                                 x = 'MAR_EF_ecpresentref_perc',
                                                 y = 'efper_ref',
                                                 log = F
                                       ), by=.(eftype_ref)
  ]
  
  #Assess error based on percentage missing
  ggplot(eftab_gefis, aes(x=mar_spe_gefisref, y=ef_spe_gefisref_probable)) +
    geom_point() +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0)
  
  ggplot(eftab_gefis, aes(x=mar_spe_gefisref, y=ef_spe_gefisref_refemc)) +
    geom_point() +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0)
  
  
  
  #Comparisons
  #Hydrology: (compare based on what ref is)
  # GEFIS~ref, GEFIS~WaterGAP, GEFIS~gauging stations
  # WaterGAP~ref, WaterGAP~gauging stations
  # ref~gauging stations
  
  #EMC class:
  #GEFIS~ref
  #HydroATLAS predictor variables~ref
  
  #E-flow as %
  # GEFIS~ref with GEFIS' original class, GEFIS~ref with ref class,
  # Gauging station Smakthin with ref class vs. ref; Gauging station Smakthin vs. GEFIS
  
  
  #Write about mismatches between routed threat and MAR
  #Compare point-base dpour EMC to ref EMC
}


