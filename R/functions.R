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
    'gwt_cm_cav',
    'inu_pc_umn',
    'inu_pc_umx',
    'inu_pc_cmn',
    'lka_pc_cse',
    'lka_pc_use',
    'dor_pc_pva', #anthropogenic - degree of regulation
    'gwt_m_cav',
    'ele_pc_rel',
    'ele_mt_cav',
    'ele_mt_uav',
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
    'pac_pc_cse', 
    'pac_pc_use',
    'swc_pc_uyr',
    'swc_pc_cyr',
    'swc_pc_cmn',
    'lit_cl_cmj',
    'kar_pc_use',
    'kar_pc_cse',
    'pop_ct_csu',
    'pop_ct_usu',
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
  
  in_dt2[, `:=`(
    glc_pc_u16 = as.numeric(glc_pc_u16),
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
    ari_ix_uav = scale_x_sqrt(expand=c(0,0), limits=c(0,20), breaks=c(0, 0.5, 1, 5)),
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
    lka_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    crp_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    pet_mm_uyr = scale_x_continuous(expand=c(0,0)),
    sdis_ms_uyr = scale_x_continuous(expand=c(0,0)),
    snw_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    #run_mm_cyr = scale_x_continuous(expand=c(0,0)),
    swc_pc_uyr = scale_x_continuous(labels=percent_format(scale=1),
                                    expand=c(0,0)),
    pac_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              labels=percent_format(scale=1),
                              expand=c(0,0)),
    ppd_pk_uav = scale_x_log10(breaks=c(0.001, 0.1, 1, 10^2,
                                        10^(0:log10(max(in_df$ppd_pk_uav)))),
                               labels=c(0, 0.1, 1, 10^2,
                                        10^(0:log10(max(in_df$ppd_pk_uav)))),
                               expand=c(0,0)),
    urb_pc_use = scale_x_sqrt(breaks=c(0, 5, 20, 50, 100),
                              limits=c(0,50),
                              labels=percent_format(scale=1),
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
    ire_pc_use = scale_x_sqrt(labels=percent_format(scale=1),
                              breaks=c(0, 5, 20, 50, 100),
                              limits=c(0,50),
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
  
  coordcart <- lapply(varstoplot, function(var) {
    coord_cartesian(xlim=as.data.table(in_df)[, c(min(get(var)+0.001, na.rm=T),
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
    xlim=c(0, 5))
  coordcart[['lka_pc_use']] <-  coord_cartesian(
    xlim=c(0, 50))
  coordcart[['urb_pc_use']] <-  coord_cartesian(
    xlim=c(0, 50))
  coordcart[['ire_pc_use']] <-  coord_cartesian(
    xlim=c(0, 50))
  
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
    
  # } else if (vartoplot == "glc_pc_u16") {
  #   rivclz <- in_rivdt[, sum(LENGTH_KM)/in_rivdt[,sum(LENGTH_KM)],
  #                      by=glc_pc_u16]
  #   gclz <- in_sitedt[,.N/in_sitedt[,.N],by=glc_pc_u16]
  #   bindclz <- rbind(rivclz, gclz, idcol='source')%>%
  #     setnames(c( 'source', vartoplot, 'density'))
  #   
  #   penvhist <- ggplot(bindclz, aes_string(x=vartoplot, y='density')) +
  #     geom_bar(aes(fill=as.factor(source)), stat='identity',
  #              position = 'identity', alpha=1/2, width=.6) +
  #     scale_fill_manual(values=c('#2b8cbe', '#dd3497'))
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
          plot.background = element_blank(),
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
  
  outstats <- dt[, list(#pearsonr = round(cor(get(y), get(x)), 3),
                        mae = round(Metrics::mae(get(y), get(x)), 2),
                        smape = round(Metrics::smape(get(y), get(x)), 2),
                        pbias = round(Metrics::percent_bias(get(y), get(x)), 2),
                        rsq = round(summary(mod)$r.squared, 3),
                        rsq_nooutliers = round(summary(mod_nooutliers)$r.squared, 3),
                        sig_coef = if (nrow(summary(mod)$coefficients) > 1) {
                          round(summary(mod)$coefficients[2,4], 3)} else {NaN},
                        n_total = .N,
                        noutliers = nrow(outrows)
  )
  ,]
  
  return(outstats)
}

#------ bin_dt -----------------
#' Bin data.table
#'
#' Bins a data.table over a numeric column.
#'
#' @param in_dt \link[data.table]{data.table} to bin.
#' @param binvar (character) column that will be used to define bins.
#' @param binfunc (character) binning approach. One of 'manual', 'equal_length', 'equal_freq'.
#' @param binarg (numeric) binning argument, depends on binning approach (\code{binfunc}).
#' @param bintrans (character or numeric) transformation of \code{binvar}, default is NULL.
#' @param ndigits (integer) number of decimals to keep for displaying formatted bin limits
#' @param na.rm (logical) whether to include NAs.
#' @param valuevar (character) if na.rm = FALSE, column to use to detect NAs and remove records.
#'
#' @return input \link[data.table]{data.table} with four new columns:
#' \itemize{
#'   \item bin - bin number (1 being the lowest value bin)
#'   \item bin_lmin - bin lower limit
#'   \item bin_lmax - bin higher limit
#'   \item bin_lformat - formatted character of bin limits
#'   (format: \code{round(bin_lmin, ndigits) - round(bin_lmax, ndigits))
#' }
#'
#' @details inspired from [rbin package](https://github.com/rsquaredacademy/rbin).
#' Differences include that it concentrates all binning approaches within a single
#' function and works on a data.table.
#'
#' binfunc: \cr
#' \itemize{
#'   \item 'manual' - bin continuous data manually. \code{binarg} sets the inner bin limits,
#'  such that the final table will have \code{length(binarg) + 1} bins. The lower end of the
#'  first bin is automatically set to be the minimum value in \code{binvar} and the upper end of
#'  the last bin is set to be the maximum value in \code{binvar}
#'
#'   \item 'equal_length' - Bin continuous data such that each bin has the same \code{binvar} interval length.
#'   If \code{bintrans} is not \code{NULL}, then interval length is computed on transformed scale.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#'
#'   \item 'equal_freq' - Bin continuous data such that each bin has the same number of records.
#'   \code{binarg} (automatically rounded to the nearest integer) sets the number of bins.
#' }
#'
#' bintrans: can either be 'log' (for natural log) or a numeric exponent to transform
#' according to x^bintrans.
#'
#' @seealso for examples, see applications in \code{\link{bin_misclass}},
#' \code{\link{eval_watergap}}, \code{\link{tabulate_globalsummary}},
#' \code{\link{formathistab}}, \code{\link{compare_us}}
#'
#' @export
bin_dt <- function(in_dt, binvar, binfunc, binarg,
                   bintrans=NULL, ndigits=2,
                   na.rm=FALSE, valuevar=NULL) {
  #Inspired from rbin, adapted to dt and simplified
  in_dt <- copy(in_dt)
  
  el_freq <- function(byd, bins) {
    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    append(min(byd, na.rm = TRUE), min(byd, na.rm = TRUE) + (bin_length * seq_len(bins)))[1:bins]
    
  }
  
  eu_freq <- function(byd, bins) {
    bin_length <- (max(byd, na.rm = TRUE) - min(byd, na.rm = TRUE)) / bins
    ufreq      <- min(byd, na.rm = TRUE) + (bin_length * seq_len(bins))
    n          <- length(ufreq)
    ufreq[n]   <- max(byd, na.rm = TRUE) + 1
    return(ufreq)
    
  }
  
  binvar_orig <- copy(binvar)
  
  #Remove NAs
  if (na.rm) {
    in_dt <- in_dt[!is.na(get(eval(valuevar))),]
  }
  
  #Transform data if trans
  if (!is.null(bintrans)) {
    transvar <- paste0(binvar, '_bintrans')
    if (bintrans == 'log') {
      nneg <- in_dt[get(eval(binvar)) <= 0, .N]
      warning(paste0('There are ', nneg, ' records with', binvar, ' <= 0...',
                     'removing them for log transformation'))
      in_dt[, eval(transvar) :=
              log(get(eval(binvar)))]
      
    } else if (is.numeric(bintrans)) {
      in_dt[, eval(transvar) :=
              get(eval(binvar))^eval(bintrans)]
      
    }
    binvar = transvar
  }
  
  byd <- in_dt[, get(eval(binvar))]
  
  if (binfunc == 'manual') {
    l_freq    <- append(min(byd), binarg)
    u_freq    <- c(binarg, (max(byd, na.rm = TRUE) + 1))
    bins      <- length(binarg) + 1
  }
  
  if (binfunc == 'equal_length') {
    bins = round(binarg)
    l_freq    <- el_freq(byd, bins)
    u_freq    <- eu_freq(byd, bins)
  }
  
  if (binfunc == 'equal_freq') {
    bins = round(binarg)
    bin_prop     <- 1 / bins
    bin_length   <- in_dt[, round(.N/bins)]
    first_bins   <- (bins - 1) * bin_length
    residual     <- in_dt[, .N - first_bins]
    bin_rep      <- c(rep(seq_len((bins - 1)), each = bin_length),
                      rep(residual, residual))
    l_freq        <- c(1, (bin_length * seq_len((bins - 1)) + 1))
    u_freq       <- c(bin_length * seq_len((bins - 1)), in_dt[,.N])
    setorderv(in_dt, cols= eval(binvar))
    in_dt[, binid := .I]
    binvar = 'binid'
  }
  
  for (i in seq_len(bins)) {
    in_dt[get(eval(binvar)) >= l_freq[i] & get(eval(binvar)) < u_freq[i],
          bin := i]
    in_dt[bin == i, `:=`(bin_lmin = min(get(eval(binvar_orig)), na.rm=T),
                         bin_lmax = max(get(eval(binvar_orig)), na.rm=T))] %>%
      .[bin == i, bin_lformat := paste(round(bin_lmin, ndigits),
                                       round(bin_lmax, ndigits),
                                       sep='-')]
    
    if (i == bins) {
      in_dt[get(eval(binvar)) == u_freq[i],  bin := i]
    }
  }
  
  if (binfunc == 'equal_freq') {in_dt[, binid := NULL]}
  
  return(in_dt)
}
#------ label_manualbins ------------
#' Label manual bins
#'
#' Utility function: label bin limits for \code{\link{formathistab}}.
#'
#' @param binarg (character vector) Arguments for bin_dt manual
#' @param minval (numeric) value to set for lower limit of first bin.
#'
#' @return vector of labels
#'
#' @export
label_manualbins <- function(binarg, minval, digits=1) {
  minlabel <- paste(round(minval, digits),
                    binarg[1], sep=" - ")
  otherlabels <- mapply(function(x, y) {paste(x, y-1, sep=" - ")},
                        binarg[1:(length(binarg)-1)], binarg[2:length(binarg)])
  return(c(minlabel, otherlabels))
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
  #print(in_startID)
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
#------ format_eftab -------------
#Clean up and format master table
format_eftab <- function(in_efp) {
  
  efpcopy <- copy(as.data.table(in_efp))
  
  #Change column names
  namedt <-  rbindlist(
    list(
      list("E-flow Location Name/No.", 'locname'),
      list("E-flow Method/Model Type", 'eftype_ref'), 
      list("E-flow Method/Model Name", 'efname_ref'), 
      list("System to Determine Ecological Condition", 'ecname_ref'), 
      list("Ecological Condition(s) in Present Day", 'ecpresent_ref'),
      list("Ecological Condition(s) Set as E-flow Objective for Future Management", 'ecfuture_ref'),
      list("Type of Hydrological Regime Used to Calculate E-flow", 'hydrotype_ref'),
      list("Natural/Naturalised Mean Annual Runoff at E-flow Location", 'mar_ref'),
      list("E-flow Requirement as Volume", 'efvol_ref'),
      list("E-flow as % of Natural/Naturalised Mean Annual Runoff", 'efper_ref')
    )
  ) %>% setnames(c('old', 'new'))
  
  setnames(efpcopy, namedt$old, namedt$new)
  
  efcopy <- clean_names(efpcopy) #Clean remaining names
  
  #Format cols
  charcols <- names(efpcopy)[efpcopy[,sapply(.SD, is.character)]]
  
  catcolstonormalize <- c('eftype_ref',
                          'hydrotype_ref')
  
  #Clean it up 
  efp_format <- efpcopy %>%
    .[, (charcols) := lapply(.SD, function(x) str_trim(x, side = 'both')), 
      .SDcols = charcols] %>%
    .[, (catcolstonormalize) := lapply(.SD, function(x) str_to_sentence(x)), 
      .SDcols = catcolstonormalize] %>%
    .[, ecpresent_ref := gsub("\\s*[(]*[nN]ot sure whether present day[)]*", "", 
                              ecpresent_ref, perl=T)] %>%
    .[!(EFUID %in% c(41)),] %>%  #Exclude EFA for Senegal river at Manantali. This is the amount only for one flooding event in the year
    .[!((Point_db == 'IWMI_3' & UID_Mathis %in% seq(603, 651)))] #Remove studies from Poland and Greece that for now cannot be used in analysis
  
  
  efp_format[, (charcols) := lapply(
    .SD,function(x){ifelse(x %in% c('', "Na", "NA", "N/A", 'Not specified', 'None'),
                           NA, x)}), 
    .SDcols = charcols]
  
  #UID_Mathis_Original is not unique as it is associated with each regional dataset
  #EFUID not unique anymore either because had to duplicate upper Ganga assessment for dorught and normal years

  #Convert mar_ref, efvol_ref, and efper_ref to numeric
  numcols <- c('mar_ref', 'efvol_ref', 'efper_ref')
  efp_format[, (numcols) := lapply(.SD, as.numeric), .SDcols = numcols]
  
  #For Brazil (double check that all records are in this case):
  #Convert Mean long term flow (m3/s) to Natural/Naturalised Mean Annual Runoff at E-flow Location (106m3)
  #Format efvol_ref
  efp_format[Country == 'Brasil', mar_ref := efvol_ref*2]

    
  unique(efp_format$ef_unit)
  #Make sure that all MAR values are in the same unit
  #31.5576 is the number of seconds in a year / 10^6
  efp_format[mar_unit == "m3 s-1", 
      `:=` (mar_ref = mar_ref*31.5576,
          mar_unit = "10^6m3 y-1")]
  efp_format[ef_unit == "m3 s-1", 
             `:=` (efvol_ref = efvol_ref*31.5576,
                   ef_unit = "10^6m3 y-1")]
  efp_format[mar_unit == "ML/day", 
             `:=` (mar_ref = mar_ref*365.25/1000,
                   mar_unit = "10^6m3 y-1")]
  efp_format[ef_unit == "ML/y", 
             `:=` (efvol_ref = efvol_ref/1000,
                   ef_unit = "10^6m3 y-1")]
  
  
  #There is no case when different MARs exist for the same sites
  #(despite the different scenarios)
  duplis <- efp_format[duplicated(paste0(Point_db, UID_Mathis)) | 
                         duplicated(paste0(Point_db, UID_Mathis), fromLast = T)
                         ,]

  ##############################################################################################Double-check why so many NAs in Victoria 

  efp_format[is.na(efper_ref) & !is.na(efvol_ref) & !is.na(mar_ref), 
             efper_ref := 100*efvol_ref/mar_ref]
  efp_format[is.na(efvol_ref) & !is.na(efper_ref) & !is.na(mar_ref), 
             efvol_ref := efvol_ref*mar_ref/100]

  #Standardize EMC field
  unique(efp_format$ecpresent_ref)
  efp_format[, ecpresent_ref_formatted := ecpresent_ref]
  efp_format[ecpresent_ref_formatted == 'I and II class', ecpresent_ref_formatted := 'A/B']
  efp_format[ecpresent_ref_formatted %in% c('Very good', 'High'), ecpresent_ref_formatted := 'A']
  efp_format[ecpresent_ref_formatted == 'Good', ecpresent_ref_formatted := 'B']
  efp_format[ecpresent_ref_formatted == 'Moderate', ecpresent_ref_formatted := 'C']
  efp_format[ecpresent_ref_formatted == 'Poor', ecpresent_ref_formatted := 'D']
  efp_format[ecpresent_ref_formatted == 'Bad', ecpresent_ref_formatted := 'E']
  
  #Format Australian system
  aus_scale_dt = data.table(
    aus_cat = c('excellent', 'good', 'moderate', 'poor', 'very poor', 'insufficient data'),
    cat_min = c(71, 51, 31, 11, 0, -99),
    cat_max = c(100, 70, 50, 30, 10, -99),
    equivalent_cat = c('A', 'B', 'C', 'D', 'E', NA)
  )
  
  efp_format[, 
             ecpresent_ref_formatted :=  fcase(
               as.numeric(ecpresent_ref_formatted) > 70, 'A',
               as.numeric(ecpresent_ref_formatted) > 50, 'B',
               as.numeric(ecpresent_ref_formatted) > 30, 'C',
               as.numeric(ecpresent_ref_formatted) > 10, 'D',
               as.numeric(ecpresent_ref_formatted) > 0, 'E',
               is.na(as.numeric(ecpresent_ref_formatted)), ecpresent_ref_formatted
             )
  ]
  
  #efp_format[grepl('Â', ecpresent_ref_formatted), ecpresent_ref_formatted := NA] #str_trim(gsub('Â', '', ecpresent_ref_formatted), side = 'both')
  efp_format[ecpresent_ref_formatted %in% c('A', 'A/B', 'B', 'B/C', 'C', 'C/D', 
                                            'D', 'D/E', 'E'),
               ecpresent_reftype := 'Standard'] 
  efp_format[ecpresent_ref_formatted %in% c('EWW', 'MaxW', 'MinW', 'MinW1', 
                                            'MinW2', 'RI', 'EWS', 'ET', 
                                            'Other', 'SR', 'MinD'),
               ecpresent_reftype := 'Other'] 
  
  rectdf <- data.table(
    xmin = rep(-Inf, 4),
    xmax = rep(Inf, 4),
    ymin = c(0, 0.25, 0.35, 0.5, 0.75),
    ymax = c(0.25, 0.35, 0.5, 0.75, 1),
    label = c('E', 'D', 'C', 'B', 'A')
  )
  
  #Compute numeric value for ordinal class
  efp_format[ecpresent_reftype == 'Standard', `:=`(
    ecpresent_ref_simplemin = fcase(
      ecpresent_ref_formatted == 'A/B', 'A',
      ecpresent_ref_formatted == 'B/C', 'B',
      ecpresent_ref_formatted == 'C/D', 'C',
      ecpresent_ref_formatted == 'D/E', 'D',
      ecpresent_ref_formatted %in% c('A', 'B', 'C', 'D', 'E'), ecpresent_ref_formatted
    ),
    ecpresent_ref_simplemax = fcase(
      ecpresent_ref_formatted == 'A/B', 'B',
      ecpresent_ref_formatted == 'B/C', 'C',
      ecpresent_ref_formatted == 'C/D', 'D',
      ecpresent_ref_formatted == 'D/E', 'E',
      ecpresent_ref_formatted %in% c('A', 'B', 'C', 'D', 'E'), ecpresent_ref_formatted
    ),
    ecpresent_ref_num = fcase(
      ecpresent_ref_formatted == 'A', 0.125,
      ecpresent_ref_formatted == 'A/B', 0.25,
      ecpresent_ref_formatted == 'B', 0.375,
      ecpresent_ref_formatted == 'B/C', 0.5,
      ecpresent_ref_formatted == 'C', 0.575,
      ecpresent_ref_formatted == 'C/D', 0.65,
      ecpresent_ref_formatted == 'D', 0.7,
      ecpresent_ref_formatted == 'D/E', 0.75,
      ecpresent_ref_formatted == 'E', 0.80
    )
  )]
  return(efp_format)
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
      'ORD_STRA', 'PFAF_ID05', in_predvars[, varcode]
    )
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
    "dis_m3_pyr", "UPLAND_SKM", "clz_cl_cmj", 
    "ari_ix_uav", "tmp_dc_uyr", "lka_pc_use",
    "for_pc_use", "crp_pc_use", "pac_pc_use",
    "ppd_pk_uav", "urb_pc_use", "ire_pc_use"
  )
  
  if ("dis_m3_pyr" %in% varstoplot_hist) {
    setDT(in_rivernetwork)[, dis_m3_pyr := dis_m3_pyr + 1]
    setDT(in_sitedt)[, dis_m3_pyr := dis_m3_pyr + 1]
  }
  
  #Remove duplicated sites
  in_sitedt <- unique(in_sitedt, by=c('POINT_X', 'POINT_Y'))
  
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
  do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=4, 
                               vp=viewport(width=0.97, height=0.97))
          )
}

#------ country_summary  -------------
#Histogram of number of sites x country x EFA_type
country_summary <- function(in_tab) {
  in_tab[, countryN := .N, by=Country]
  
  in_tab[is.na(eftype_ref), eftype_ref := 'None specified']
  
  country_hist <- ggplot(in_tab, aes(x=reorder(Country, countryN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=..count.., label=..count..), stat='count', geom="text", hjust=-.5) +
    coord_flip(clip='off') +
    scale_x_discrete('Country') + 
    scale_y_continuous('Number of EFAs') +
    scale_fill_brewer('Type of EFA', palette='Set2', na.value="grey") +
    theme_classic() +
    theme(#axis.title.y = element_text(vjust=-2),
          plot.margin = unit( c(0.1, 0.6, 0.1, 0.1), 'cm'),
          legend.position = c(0.7, 0.2))
  
  return(country_hist)
}

#------ ecoregions_summary -------------
#Histogram of number of sites x country x EFA_type
ecoregions_summary <- function(in_tab) {
  in_tab[, fecN := .N, by=fec_cl_cmj]
  
  ecoregions_hist <- ggplot(in_tab, aes(x=reorder(fec_cl_cmj, fecN))) + 
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

#------ map_ef -------------
#-- Plot studies
map_ef <- function(in_efp) {
  in_efp <- unique(in_efp, by=c('POINT_X', 'POINT_Y'))
                   
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
  
  clustcentro <- in_efp[, `:=`(cluster = cutree(studyclus, h=1000000))] %>%#,
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
    scale_size_continuous(name=str_wrap('Number of EFAs within 1000 km', 20),
                          range=c(1, 5), breaks=c(1, 5, 10, 50, 100)) +
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

#------ link_EFtoGRDC ------------------------------------------
link_EFtoGRDC <- function(in_eftab_gefis, inp_riveratlas, inp_grdcp) {
  in_eftab_gefis <- unique(in_eftab_gefis, by=c('POINT_X', 'POINT_Y'))
  
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

#------ compare_EFtoGRDC ------------------------------------------
compare_EFtoGRDC <- function(in_eftab_grdc_join) {
  
  "mar_ref"
  "dis_m3_pyr"
  ""
  
}


#------ join_efp_to_efmod ------------------------------------------------------
join_efp_to_efmod <- function(in_eftab, in_path_efp_mod) {
  efp_mod <- fread(in_path_efp_mod)
  
  #Format NAs
  efp_mod[, (names(efp_mod)) := lapply(
    .SD,function(x){ifelse(x %in% c(-9999, '-9999', ''),
                           NA, x)}), 
    .SDcols = names(efp_mod)]
  
  #Remove extraneous columns
  colstoremove <- c("V1", "OBJECTID", "merge_efextract_df", "eftype_2_x",
                    "Point_shift_mathis", "Comment_mathis", "POINT_X", "POINT_Y",
                    "path", "climate_scenario", "human_scenario", 'variable',
                    'eftype_2')
  efp_mod[, (colstoremove):=NULL]

  #Average monthly values
  efp_mod_format <- efp_mod[, list(value = mean(value)), 
                            by=c("UID_Mathis", "Point_db","ghm", "gcm", "var", "res",
                                 "run", "eftype", "eftype_format", "emc")]
  
  binlabels <- label_manualbins(binarg=c(10, 100, 1000, 10000, 10000000),
                                minval=min(in_eftab$up_area_skm_15s)) %>%
    data.table(bin_label = .,
               bin=seq(1,6))
  
  eftab_modjoin <- merge(in_eftab,
                         efp_mod_format,
                         by=c("UID_Mathis", "Point_db"),
                         allow.cartesian=TRUE) %>%  #Only keep maf estimates
    bin_dt(binvar = "up_area_skm_15s",
           binfunc='manual',
           binarg = c(10, 100, 1000, 10000, 10000000)) %>%
    merge(binlabels, by='bin') %>%
    .[eftype_format != 'mmf',] 
  
  #Correctly scale values for isimp2b
  eftab_modjoin[res=="0.5 arc-deg" & var=='qtot', value := value/1000]
  
  ##################################################################################Problem with Smahktin computations for isimp2b discharge
  
  return(eftab_modjoin)
}

#------ compare_hydrology ------------------------------------------
compare_hydrology <- function(in_efp_efmod_join) {
  eftab_maf <- in_efp_efmod_join[
    (Comment_mathis != 'Not well represented in HydroRIVERS') | 
      (is.na(Comment_mathis)),]  %>% #Remove those that cannot be compared
    .[!is.na(mar_ref),] %>% #remove sites without MAR value
    unique(by=c('POINT_X', 'POINT_Y', 'run', 'var', 'emc', 'eftype_format')) %>% #Keep unique sites
    .[eftype_format == 'maf']  %>%
    .[,`:=`(mar_ref = mar_ref/31.5576, 
            mar_unit = "m3 s-1")] 
  #--------------------------Compare reference and model MAR for discharge -----
  ggplot(eftab_maf[var=='dis',], aes(x=mar_ref, y=value)) +
    geom_point() +
    scale_x_log10(name= expression('Observed mean annual discharge'~(m^3~y^-1))) +
    scale_y_log10(name= expression('Predicted mean annual discharge'~(m^3~y^-1))) +
    geom_abline () +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic()
    
  ggplot(eftab_maf[var=='qtot',], aes(x=mar_ref, y=value)) +
    geom_point() +
    scale_x_log10(name= expression('Observed mean annual discharge'~(m^3~y^-1))) +
    scale_y_log10(name= expression('Predicted mean annual discharge'~(m^3~y^-1))) +
    geom_abline () +
    geom_smooth() +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic()
  
  
  qstats_all <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      x = 'mar_ref',
      y = 'value',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm')
  ]
  
  qstats_da <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      x = 'mar_ref',
      y = 'value',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'bin_label')
  ] %>%
    setnames('bin_label', "Drainage area")
  
  qstats_clz <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      x = 'mar_ref',
      y = 'value',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'clz_cl_cmj')
  ] 
  
  qstats_country <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      x = 'mar_ref',
      y = 'value',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'Country')
  ] 

  
  fwrite(qstats_all, file.path(resdir, 
                               paste0('qstats_all',
                                      format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_da, file.path(resdir, 
                           paste0('qstats_da',
                                  format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_clz, file.path(resdir, 
                              paste0('qstats_clz',
                                     format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_country, file.path(resdir, 
                               paste0('qstats_country',
                                      format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  
  return(list(plot = composite_plot, 
              table_all = qstats,
              table_clz = qstats_clz))
}

#------ compare_EMC ------------------------------------------
compare_EMC <- function(in_eftab, in_riveratlas_varsdt) {
  # unique(in_eftab$Ecological)
  # in_eftab$ecpresent_ref
  # in_eftab$EMC_10Variable_2
  
  #Prepare dt for plot background showing the difference classes
  rectdf <- data.table(
    xmin = rep(-Inf, 4),
    xmax = rep(Inf, 4),
    ymin = c(0, 0.25, 0.5, 0.65, 0.75),
    ymax = c(0.25, 0.5, 0.65, 0.75, 1),
    label = c('A', 'B', 'C', 'D', 'E')
  )
  
  #Test differences among ckasses
  compute_tukeyletters <- function(in_dt, ectype, var) {
    dtsub <- in_dt[ecpresent_reftype == ectype]
    aov_out<- aov(get(var) ~ ecpresent_ref_formatted, data=dtsub)    
    tukey_out <- TukeyHSD(aov_out)
    cld <- multcompView::multcompLetters4(aov_out, tukey_out)
    as.data.table(as.data.frame.list(cld$ecpresent_ref_formatted), keep.rownames = T)
  }
  
  #Check correspondence in WS
  tukeyletters_ws <- lapply(in_eftab[!is.na(ecpresent_reftype), 
                                  unique(ecpresent_reftype)], function(ectype) {
                                    compute_tukeyletters(in_dt = in_eftab, 
                                                         ectype = ectype,
                                                         var = 'EMC_10Variable_2_mean')
                                  }) %>%
    rbindlist(fill=T) %>%
    merge(unique(in_eftab[, .(ecpresent_ref_formatted, ecpresent_reftype)]),
          by.x = 'rn', by.y = 'ecpresent_ref_formatted') %>%
    merge(in_eftab[, .N, by=ecpresent_ref_formatted],
          by.x = 'rn', by.y = 'ecpresent_ref_formatted') %>%
    .[, plabel := fifelse(ecpresent_reftype == 'Standard', 
                          paste0(Letters, ' (', N, ')'),
                          paste0('(', N, ')')
                          )] %>%
    rbind(data.table(rn = c('None', NA),
                     plabel = c(in_eftab[ecpresent_ref_formatted == 'None', .N],
                                in_eftab[is.na(ecpresent_ref_formatted), .N])
                     ), 
                     fill=T)
  


  in_eftab[, ecpresent_reftype := factor(ecpresent_reftype,
                                          levels = c("Standard", "Other", "NA"))
  ]
  
  ECws_boxplot <- ggplot(in_eftab) + 
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax, 
                               ymin=ymin, ymax=ymax, fill = label),
              color=NA) +
    geom_boxplot(aes(x=ecpresent_ref_formatted, y=EMC_10Variable_2_mean)) +
    geom_text(data=tukeyletters_ws, aes(x=rn, y=0.7, label=plabel), hjust = 0) +
    # geom_text(data=rectdf, aes(x='A', y=ymin, label=label),
    #           vjust = -1.5, hjust=2) +
    scale_fill_brewer(name = str_wrap('Average GEFIS ecological class in catchment', 20),
                      palette='RdYlBu', direction = -1) +
    scale_y_continuous('Average Incident Biodiversity Threat (IBT) index in catchment', limits=c(0,1), 
                       expand = c(0,0), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete('Present-day ecological class from local EFA',
                     limits=rev) +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12)) +
    ggforce::facet_col(factor(ecpresent_reftype,levels = c("Standard", "Other", "NA"))~.,
                       scales = 'free_y', space='free')
  

  #Check correspondence at pourpoint
  tukeyletters_p <- lapply(in_eftab[!is.na(ecpresent_reftype), 
                                     unique(ecpresent_reftype)], function(ectype) {
                                       compute_tukeyletters(in_dt = in_eftab, 
                                                            ectype = ectype,
                                                            var = 'EMC_10Variable_2')
                                     }) %>%
    rbindlist(fill=T) %>%
    merge(unique(in_eftab[, .(ecpresent_ref_formatted, ecpresent_reftype)]),
          by.x = 'rn', by.y = 'ecpresent_ref_formatted') %>%
    merge(in_eftab[, .N, by=ecpresent_ref_formatted],
          by.x = 'rn', by.y = 'ecpresent_ref_formatted')  %>%
    .[, plabel := fifelse(ecpresent_reftype == 'Standard', 
                          paste0(Letters, ' (', N, ')'),
                          paste0('(', N, ')')
    )]
  
  ECp_boxplot <- ggplot(in_eftab) + 
    geom_rect(data=rectdf, aes(xmin=xmin, xmax=xmax, 
                               ymin=ymin, ymax=ymax, fill = label),
              color=NA) +
    geom_boxplot(aes(x=ecpresent_ref_formatted, y=EMC_10Variable_2)) +
    geom_text(data=tukeyletters_p, aes(x=rn, y=0.7, label=plabel), hjust = 0) +
    # geom_text(data=rectdf, aes(x='A', y=ymin, label=label),
    #           vjust = -1.5, hjust=2) +
    scale_fill_brewer(name = str_wrap('Average GEFIS ecological class in catchment', 20),
                      palette='RdYlBu', direction = -1) +
    scale_y_continuous('Average Incident Biodiversity Threat (IBT) in catchment', limits=c(0,1), 
                       expand = c(0,0), breaks=c(0, 0.25, 0.5, 0.75, 1)) +
    scale_x_discrete('Present-day ecological class from local EFA') +
    theme_classic() +
    coord_flip() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=12)) +
    ggforce::facet_col(factor(ecpresent_reftype,levels = c("Standard", "Other", "NA"))~.,
                       scales = 'free_y', space='free')
  
  #Compute classification statistics + confusion matrix
  confumat_min <- in_eftab[
    ecpresent_reftype == 'Standard', 
    caret::confusionMatrix(as.factor(ecpresent_gefisws), 
                           as.factor(ecpresent_ref_simplemin))
  ]
  confumat_max <- in_eftab[
    ecpresent_reftype == 'Standard', 
    caret::confusionMatrix(as.factor(ecpresent_gefisws), 
                           as.factor(ecpresent_ref_simplemax))
  ]
  
  #Check whether difference is due to missing data (counter-intuitive pattern!)
  # ggplot(in_eftab, aes(x=EMC_boolean_ratio, y=abs(ecpresent_refgefisdiff))) + 
  #   geom_point() +
  #   theme_bw()
  # # + 
  #   geom_quantile(quantiles=c(0.1, 0.5, 0.9))
  
  
  #Compare HydroATLAS predictors/variables to ref classes
  #Check protected area %
  eftab_ecstressorsmelt <- melt(
    in_eftab[, .(dor_pc_pva, crp_pc_use, pst_pc_use, ire_pc_use,
                 ppd_pk_uav, urb_pc_use,
                 ecpresent_ref_formatted, ecpresent_reftype, 
                 EMC_10Variable_2_mean, EFUID
    )],
    id.vars = c('EFUID', 'ecpresent_ref_formatted', 'ecpresent_reftype', 
                'EMC_10Variable_2_mean')
  ) %>%
    merge(in_riveratlas_varsdt, by.x = 'variable', by.y = 'varcode', all.x=T) %>%
    .[, varname := factor(varname, 
                          levels = c(
                            "Cropland Extent watershed Spatial extent (%)",
                            "Pasture Extent watershed Spatial extent (%)" ,
                            "Irrigated Area Extent (Equipped) watershed Spatial extent (%)",
                            "Degree of Regulation pour point Value",                 
                            "Population Density watershed Average",                        
                            "Urban Extent watershed Spatial extent (%)")
    )]
  
  drivers_boxplot <- ggplot(
    eftab_ecstressorsmelt[ecpresent_reftype == 'Standard',], 
    aes(x=ecpresent_ref_formatted, y=value)) +
    geom_boxplot() +
    facet_wrap(~varname, scales='free_y') +
    labs(x='Present-day ecological management class from local EFA',
         y='Value across upstream drianage area (catchment)') +
    theme_classic()
  
  # ggplot(eftab_ecstressorsmelt, 
  #        aes(x=value, y=EMC_10Variable_2_mean)) +
  #   geom_point() + 
  #   geom_smooth() +
  #   facet_wrap(~varname, scales='free') +
  #   theme_bw()
  
  return(list(
    ECws_boxplot = ECws_boxplot, 
    ECp_boxplot = ECp_boxplot,
    drivers_boxplot = drivers_boxplot,
    confumat_min = confumat_min,
    confumat_max = confumat_max
  ))
}

#------ compare_EFestimate ------------------------------------------
compare_EFestimate <- function(in_efp_efmod_join) {
  eftab <- in_efp_efmod_join[
    (Comment_mathis != 'Not well represented in HydroRIVERS') | 
      (is.na(Comment_mathis)),]  %>% #Remove those that cannot be compared
    merge(dcast(.[eftype_format == 'maf',],  #Create new column for MAF rather than as a record
                EFUID+gcm+ghm+var~eftype_format, 
                value.var = 'value'),
          by=c('EFUID', 'gcm', 'ghm', 'var')
    )%>%
    .[!(eftype_format %in% c("maf", "mmf")),] %>% #Remove maf records
    .[, efper_mod := fifelse(maf >0, 100*value/maf, 0)] %>% #COmpute percentage EF
    .[value>0, ef_refmod_diff := (value-efper_mod)/value] %>%
    .[ef_unit == "10^6m3 y-1", `:=`(
      ef_unit = "m3 s-1",
      efvol_ref = efvol_ref/31.5576
    )] %>%
    .[efper_mod > 100, efper_mod := 100] %>% #7 cases 
    .[efper_ref <= 100,] %>%
    merge(.[, unique(.SD),  .SDcols=c('Point_db', 'UID_Mathis', 'Country')][
      , list(Ncountry=.N), by='Country'],
      by='Country')
  
  #Analysis by volume
  #########Modify by volume
  eftab_worst <- eftab[order(ef_refmod_diff),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      eftab[order(ef_refmod_diff),], 
      by=c('POINT_X', 'POINT_Y', 'gcm', 'ghm', 'var', 'eftype_format'))),
  ] 
  
  eftab_best <- eftab[order(-ef_refmod_diff),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      eftab[order(-ef_refmod_diff),], 
      by=c('POINT_X', 'POINT_Y', 'gcm', 'ghm', 'var', 'eftype_format'))),
  ] 
  
  
  efvol_stats_country <- eftab_best[!is.na(efvol_ref) & 
                                      var=='qtot'& 
                                      Ncountry > 10
                                    ,
                                    getqstats(
                                      dt = .SD,
                                      x = 'efvol_ref',
                                      y = 'value',
                                      rstudthresh= 3,
                                      log = T
                                    ), 
                                    by=c('var', 'gcm', 'ghm', 
                                         'Country', 'eftype_format')
  ] 
  
  
  EFcompare_best_vol_France_p <- ggplot(eftab_best[var=='qtot' &
                                                     Country == 'France',], 
                                        aes(x=(efvol_ref), y=value)) +
    geom_abline(alpha=1/2) +
    geom_point() + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10() + 
    scale_y_log10() +
    geom_smooth() +
    facet_grid(c("run", "eftype_format"), 
               labeller = "label_both",
               scales='free')
  plot(EFcompare_best_vol_France_p)  
  
  
  #Analysis by percentage
  

  
  
  ####################################################################################Investigate
  check <- unique(eftab[efper_ref > 100,], by="EFUID")[,
  c("EFUID", "UID_Mathis", "Point_db", "Country", "Basin", "River",
  "locname", "mar_unit", "mar_ref", "ef_unit", "efvol_ref", "efper_ref"),
  with=F]
  
  eftab[is.na(eftype_ref), eftype_ref := 'None specified']
  eftab[, Country := factor(Country)]
  
    # unique(by=c('POINT_X', 'POINT_Y', 'run', 'var', 'emc', 'eftype_format')) %>% #Keep unique sites
    # .[eftype_format == 'maf'] %>%  #Only keep maf estimates
    # .[,`:=`(mar_ref = mar_ref/31.5576, #convert mar ref to m3 s-1
    #         mar_unit = "m3 s-1")] %>%
    # bin_dt(binvar = "up_area_skm_15s",
    #        binfunc='manual',
    #        binarg = c(10, 100, 1000, 10000, 10000000)) %>%
    # merge(binlabels, by='bin')
    
    # #Create new fields for comparison
    # efp_format[, `:=`(
    #   mar_spe_gefiswg = 200*(mar_gefis - mar_nat_wg22_ls_year)/(
    #     mar_gefis + mar_nat_wg22_ls_year),
    #   mar_spe_gefisref = 200*(mar_gefis - mar_ref)/(
    #     mar_gefis + mar_ref),
    #   ef_spe_gefisref_probable = 200*(MAR_EF_Probable_perc - efper_ref)/(
    #     MAR_EF_Probable_perc + efper_ref),
    #   ef_spe_gefisref_emcref = 200*(MAR_EF_ecpresentref_perc - efper_ref)/(
    #     MAR_EF_ecpresentref_perc + efper_ref)
    # )]
    # 
    # #Create unique ID for each study and site
    # efp_format[, SiteUID := as.numeric(factor(do.call(paste0,.SD))),
    #            .SDcols=c('POINT_X', 'POINT_Y')] #Update with POINT_X and POINT_Y
  
  #When multiple scenarios and one is clearly present-day, only keep that one 
  #in_eftabsub_pre <- in_eftab[!(EFUID %in% c(2, 3, 4, 10, 35, 36, 37)),]
  

  
  
  #Compare reference EF to GEFIS EF - % and vol
  EFcompare_best_p <- ggplot(eftab_best[var=='qtot' & eftype_ref =='Holistic',], 
         aes(x=(efper_ref), y=efper_mod)) +
    geom_abline(alpha=1/2) +
    geom_point(aes(group=Country, color=Country)) + #aes(group=Country,color=Country, shape=Country)
    geom_smooth() +
    facet_grid(c("run", "eftype_format"), 
               labeller = "label_both",
               scales='free')
  
  plot(EFcompare_best_p)
  
  EFcompare_best_France_p <- ggplot(eftab_best[var=='qtot' &
                                          Country == 'France',], 
                             aes(x=(efper_ref), y=efper_mod)) +
    geom_abline(alpha=1/2) +
    geom_point() + #aes(group=Country,color=Country, shape=Country)
    geom_smooth() +
    facet_grid(c("run", "eftype_format"), 
               labeller = "label_both",
               scales='free')
  plot(EFcompare_best_France_p)  
  

  
    
    scale_shape_manual(values=rep(c(15, 16, 17, 18), ceiling(nlevels(in_eftabsub_best$Country)/4))) +
    #geom_quantile() +
    scale_x_continuous(name = 'GEFIS e-flow (% of MAR)', limits=c(0, 100)) +
    scale_y_continuous(name = 'Locally-determined e-flow (% of MAR)', limits=c(0, 100)) +
    coord_fixed(expand=c(0,0), clip='off') +
    facet_wrap(~eftype_ref, nrow=3) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(),
      #strip.background=element_rect(colour="white", fill='lightgray')
    )
  
  #in_eftabsub[!is.na(MAR_EF_ecpresentref_perc),]
  EFcompare_bestec_p <- ggplot(in_eftabsub_best[((eftype_ref != 'Habitat simulation') | is.na(eftype_ref)),],
         aes(x=MAR_EF_ecpresentref_perc, y=as.numeric(efper_ref), color=Country)) +
    geom_abline(alpha=1/2) +
    geom_point(aes(group=Country,color=Country, shape=Country)) +
    scale_shape_manual(values=rep(c(15, 16, 17, 18), ceiling(nlevels(in_eftabsub_best$Country)/4))) +
    # geom_quantile(data =in_eftabsub[!is.na(MAR_EF_ecpresentref_perc) &
    #                                   (Country == 'South Africa'),],
    #               color='black', alpha=1/2) +
    geom_smooth(method='lm') +
    scale_x_continuous(name = 'GEFIS e-flow (% of MAR)', limits=c(0, 100)) +
    scale_y_continuous(name = 'Locally-determined e-flow (% of MAR)', limits=c(0, 100)) +
    coord_fixed(expand=c(0,0)) +
    facet_wrap(~eftype_ref) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(),
      #strip.background=element_rect(colour="white", fill='lightgray')
    )
  
  #Compute statistics for different EF types and countries: make table with errors, bias, etc.
  efstats_country <- in_eftabsub_best[!(is.na(efper_ref) & !is.na(MAR_EF_Probable_perc)) &
                                   EMC_boolean_ratio > 0.3,
                                 getqstats(dt = .SD,
                                           x = 'MAR_EF_Probable_perc',
                                           y = 'efper_ref',
                                           log = F
                                 ), by=.(Country)
  ] %>%
    .[n_total < 20, `:=`(
      rsq = NA,
      rsq_nooutliers = NA,
      sig_coef = NA
    )]
  
  efstats_eftype_best <- in_eftabsub_best[!(is.na(efper_ref) & !is.na(MAR_EF_Probable_perc)) &
                                  EMC_boolean_ratio > 0.3,
                                getqstats(dt = .SD,
                                          x = 'MAR_EF_Probable_perc',
                                          y = 'efper_ref',
                                          log = F
                                ), by=.(eftype_ref)
  ] %>%
    .[n_total < 20, `:=`(
      rsq = NA,
      rsq_nooutliers = NA,
      sig_coef = NA
    )]
  
  efstats_eftype_worst <- in_eftabsub_worst[!(is.na(efper_ref) & !is.na(MAR_EF_Probable_perc)) &
                                            EMC_boolean_ratio > 0.3,
                                          getqstats(dt = .SD,
                                                    x = 'MAR_EF_Probable_perc',
                                                    y = 'efper_ref',
                                                    log = F
                                          ), by=.(eftype_ref)
  ] %>%
    .[n_total < 20, `:=`(
      rsq = NA,
      rsq_nooutliers = NA,
      sig_coef = NA
    )]
  
  
  efstats_eftype_emcref <- in_eftabsub_pre[!(is.na(efper_ref)) & !(is.na(MAR_EF_ecpresentref_perc)) &
                                         (EMC_boolean_ratio > 0.3),
                                       getqstats(dt = .SD,
                                                 x = 'MAR_EF_ecpresentref_perc',
                                                 y = 'efper_ref',
                                                 log = F
                                       ), by=.(eftype_ref)
  ] %>%
    .[n_total < 20, `:=`(
      rsq = NA,
      rsq_nooutliers = NA,
      sig_coef = NA
    )]
  
  #Analyze the impact of missing information from masking
  efmaskper_p <- ggplot(in_eftab[is.na(efper_ref)],
                        aes(x=100*(1-EMC_boolean_ratio), 
                            y=ef_spe_gefisref_probable)) +
    geom_point(aes(color=Country, shape=Country)) +
    scale_shape_manual(values=rep(c(15, 16, 17, 18), ceiling(nlevels(in_eftabsub_best$Country)/4))) +
    # geom_quantile(data =in_eftabsub[!is.na(MAR_EF_ecpresentref_perc) &
    #                                   (Country == 'South Africa'),],
    #               color='black', alpha=1/2) +
    geom_smooth(method='lm') +
    #scale_x_continuous(name = 'GEFIS e-flow (% of MAR)', limits=c(0, 100)) +
    #scale_y_continuous(name = 'Locally-determined e-flow (% of MAR)', limits=c(0, 100)) +
    coord_fixed(expand=c(0,0)) +
    facet_wrap(~eftype_ref) +
    theme_classic() +
    theme(
      panel.grid.major = element_line(),
      #strip.background=element_rect(colour="white", fill='lightgray')
    )
  
  #Examine error based on river size, MAR SPE
  #aes(group=Country,color=Country, shape=Country))
  
  return(list(
    EFcompare_originalec_p  = EFcompare_originalec_p ,
    EFcompare_bestec_p = EFcompare_bestec_p ,
    table_country =  efstats_country,
    table_eftype_best =  efstats_eftype_best,
    table_eftype_worst =  efstats_eftype_worst,
    table_emcref = efstats_eftype_emcref
  ))
}

#------ check_missing ------------------------------------------
check_masking <- function(in_eftab) {
  in_eftab[, plabel_mask := paste0(Country, 
                                   ' (', 
                                   length(unique(paste0(POINT_X, POINT_Y))), 
                                   ')'), 
           by=Country]
  
  #Show missing data
  ggplot(in_eftab, aes(x=as.factor(paste0(Country), y=1-MAR_boolean_ratio))) +
    geom_boxplot()
    
  efmaskper_p <- ggplot(in_eftab,
         aes(x=str_wrap(as.factor(plabel_mask), 20), 
             y=100*(1-EMC_boolean_ratio),
             fill=Country, color=Country)) +
    #tar_makegeom_hline(yintercept=30) +
    geom_jitter(alpha=1/3) +
    geom_boxplot(alpha=1/2) +
    scale_y_continuous(str_wrap('Percentage of the catchment area that is masked', 40), 
                       limits=c(0, 100)) +
    scale_x_discrete('Country', limits=rev) +
    coord_flip(expand=c(0,0), clip='off') +
    theme_bw() +
    theme(legend.position='none',
          plot.margin = unit( c(0.1, 0.5, 0.1, 0.1), 'cm'))
  
  
  #Assess error based on percentage missing
  ggplot(in_eftab, aes(x=mar_spe_gefisref, y=ef_spe_gefisref_probable)) +
    geom_point() +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0)
  
  ggplot(in_eftab, aes(x=mar_spe_gefisref, y=ef_spe_gefisref_refemc)) +
    geom_point() +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0)
  
  return(efmaskper_p)
}


######################### MAKE COMPARISONS WITH BEST AND WORST MATCH ######################
#Desktop methods tend to be more precautionnary and holistic assessments tend to make lower recommendations
#Generate a range of performance indices based on best and worst match when lacking A-D.
