#---------------------- Utility functions --------------------------------------
# function for number of observations 
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

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
      stat_count(aes(y=after_stat(count), label=after_stat(count)), stat='count', geom="text", hjust=-.5) +
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
      #stat_count(aes(y=after_stat(count), label=after_stat(count)), stat='count', geom="text", hjust=-.5) +
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
    ari_ix_cav = ari_ix_cav/100,
    ari_ix_uav = ari_ix_uav/100,
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
# dt <- eftab_maf[var=='qtot' & gcm=="gfdl-esm2m" & ghm=='h08',]
# predicted <- 'mar_mod'
# actual <- 'mar_ref'
# log=TRUE
# rstudthresh= 3
# mae_decomp = T

getqstats <- function(dt, predicted, actual, rstudthresh= 3,
                      mae_decomp=FALSE, log=FALSE) {
  if (log) {
    in_form = paste0('log10(', actual, '+0.1)~log10(', predicted, '+0.1)')
  } else {
    in_form = paste0(actual,'~',predicted)
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
    mae = round(Metrics::mae(get(actual), get(predicted)), 2),
    smape = 100*round(Metrics::smape(get(actual), get(predicted)), 2),
    pbias = 100*round(-Metrics::percent_bias(get(actual), get(predicted)), 2), #Take the negative of percent_bias so that (predicted-actual) / abs(actual)
    rsq = round(summary(mod)$r.squared, 3),
    rsq_nooutliers = round(summary(mod_nooutliers)$r.squared, 3),
    sig_coef = if (nrow(summary(mod)$coefficients) > 1) {
      round(summary(mod)$coefficients[2,4], 3)} else {NaN},
    n_total = .N,
    noutliers = nrow(outrows)
  )
  ,]
  
  #Compute MAE decomposition metrics from Robeson and Willmott 2023 
  #Decomposition of the mean absolute error (MAE) into systematic and unsystematic
  #components.https://doi.org/ 10.1371/journal.pone.0279774
  if (mae_decomp) {
    dt_copy <- copy(dt)
    
    mae_log10 <- dt_copy[, Metrics::mae(log10(get(actual)+0.1), 
                                   log10(get(predicted)+0.1)
    )]
    
    #Compute bias component of MAElog
    MBE <- dt_copy[, mean(log10(get(predicted)+0.1)) - mean(log10(get(actual)+0.1))] #Mean bias error
    dt_copy[, p_prime := log10(get(predicted)+0.1) - MBE] #unbiased predicted values
    bias_error_robesonwillmott = abs(MBE)
    
    #Compute proportionality error for each obs
    dt_copy[, p_hat := predict(mod)]
    dt_copy[, p_hatprime := p_hat - MBE] #unbiased predicted regression values 
    dt_copy[, proportionality_error_i := abs(p_hatprime-log10(get(actual)+0.1))] #observation-specific weight of the relative importance of proportionality error
    
    #Compute unsystematic error for each obs
    dt_copy[, unsystematic_error_i := abs(p_prime - p_hatprime)] #absolute difference between unbiased predictions and unbiased regression values
    
    mae_log10_bias <- dt_copy[
      , 
      sum(
        abs(log10(get(predicted)+0.1)-log10(get(actual)+0.1)) 
        * bias_error_robesonwillmott
        /(bias_error_robesonwillmott 
          + proportionality_error_i
          + unsystematic_error_i)
      )/.N
    ]
    
    mae_log10_proportionality_error <- dt_copy[
      , 
      sum(
        abs(log10(get(predicted)+0.1)-log10(get(actual)+0.1)) 
        * proportionality_error_i
        /(bias_error_robesonwillmott 
          + proportionality_error_i
          + unsystematic_error_i)
      )/.N
    ]
    
    mae_log10_unsystematic_error <- dt_copy[
      , 
      sum(
        abs(log10(get(predicted)+0.1)-log10(get(actual)+0.1)) 
        * unsystematic_error_i
        /(bias_error_robesonwillmott 
          + proportionality_error_i
          + unsystematic_error_i)
      )/.N
    ]        
    
    outstats <- append(
      outstats, 
      list(
        mae_log10=mae_log10,
        mae_log10_bias=mae_log10_bias, 
        mae_log10_proportionality_error=mae_log10_proportionality_error,
        mae_log10_unsystematic_error = mae_log10_unsystematic_error
      ))
  }
  
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
  
  #Simplify column names
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
    .[, (charcols) := lapply(.SD, function(x) str_trim(x, side = 'both')),  #Remove trailing spaces
      .SDcols = charcols] %>%
    .[, (catcolstonormalize) := lapply(.SD, function(x) str_to_sentence(x)),  #Make upper and lower cases consistebt
      .SDcols = catcolstonormalize] %>%
    .[, ecpresent_ref := gsub("\\s*[(]*[nN]ot sure whether present day[)]*", "", 
                              ecpresent_ref, perl=T)] %>%
    .[!(EFUID %in% c(41)),] %>%  #Exclude EFA for Senegal river at Manantali. This is the amount only for one flooding event in the year
    .[!(EFUID %in% c(40)),] %>%  #Exclude EFA for Colorado River Delta. It doesn't make sense
    .[!(Point_db == 'IWMI_1' & Country == 'Australia'),] %>%  #Remove Australian data from the original IWMI database to keep homogeneous data for the country (and because didn't QC them enough)
    .[!(Point_db == 'IWMI_1' & UID_Mathis %in% c(60, 87)),] %>%  #Deleted in upstream steps. But didn't re-run everything (not well represented in network. Toosmall and imprecise)
    .[!((Point_db == 'IWMI_3' & UID_Mathis %in% seq(603, 651))),] %>% #Remove studies from Poland and Greece that for now cannot be used in analysis
    .[!((Point_db == "Rhône_eaufrance" & UID_Mathis == 460)),] %>% #Deleted in upstream steps. But didn't re-run everything
    .[!((Point_db == 'Mexico_WWF' & UID_Mathis == 2622)),] %>% #Anomalous. The prescribed e-flow is 0.4L/s for a river with MAR of 1 m3/s
    .[!((Point_db == 'Mexico_WWF' & UID_Mathis == 1252)),] %>% #Anomalous. The Mexican model estimates a discharge of 229 m3/s for a 350 km2 basin
    .[!((Point_db == 'Mexico_WWF' & UID_Mathis == 3077)),] %>% #Anomalous. The Mexican model estimates a discharge of 670 m3/s for a 1000 km2 basin
    .[!((Point_db == 'Mexico_WWF' & UID_Mathis == 1237)),]  %>% #Deleted in upstream steps. But didn't re-run everything
    .[Point_db == 'Victoria_VEWH' & UID_Mathis == 61, mar_ref := NA] #MAR is simulated flow but not natural flow
  
  #Convert all forms of NA to R-compatible NA
  efp_format[, (charcols) := lapply(
    .SD,function(x){ifelse(x %in% c('', "Na", "NA", "N/A", 'Not specified', 'None'),
                           NA, x)}), 
    .SDcols = charcols]
  
  #UID_Mathis_Original is not unique as it is associated with each regional dataset
  
  #Convert mar_ref, efvol_ref, and efper_ref to numeric
  numcols <- c('mar_ref', 'efvol_ref', 'efper_ref')
  efp_format[, (numcols) := lapply(.SD, as.numeric), .SDcols = numcols]
  
  
  #--- Format efvol_ref ---------------------------
  #For Brazil :
  #Convert Mean long term flow (m3/s) to Natural/Naturalised Mean Annual Runoff at E-flow Location (106m3)
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
  
  efp_format[!is.na(mar_unit ) & !is.na(ef_unit) & !(mar_unit == ef_unit),] #There are no cases of differing units
  efp_format[is.na(efper_ref) & !is.na(efvol_ref) & !is.na(mar_ref), 
             efper_ref := 100*efvol_ref/mar_ref]
  efp_format[is.na(efvol_ref) & !is.na(efper_ref) & !is.na(mar_ref), 
             `:=`(efvol_ref = efper_ref*mar_ref/100,
                  ef_unit = mar_unit)]
  
  #--- Format EMC ref-----------------------------------------------------------
  
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
  
  #efp_format[grepl('Â', ecpresent_ref_formatted), ecpresent_ref_formatted := NA] 
  #str_trim(gsub('Â', '', ecpresent_ref_formatted), side = 'both')
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
layout_ggenvhist <- function(in_rivernetwork, in_sitedt, 
                             in_predvars, out_path) {
  varstoplot_hist <- c(
    "UPLAND_SKM", "clz_cl_cmj", "ari_ix_uav", 
    "tmp_dc_uyr", "for_pc_use", "crp_pc_use", 
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
  scalesenvhist <- formatscales(in_df=in_rivernetwork, 
                                varstoplot=varstoplot_hist)
  
  #Plot each facet
  penvhist_grobs <- lapply(varstoplot_hist, ggenvhist,
                           in_sitedt = in_sitedt,
                           in_rivdt = in_rivernetwork,
                           in_predvars = in_predvars,
                           scalesenvhist = scalesenvhist)
  #Add legend
  #penvhist_grobs[[length(penvhist_grobs) + 1]] <- leg
  
  #Plot
  out_plot <- do.call("grid.arrange", list(grobs=penvhist_grobs, nrow=3, 
                                           vp=viewport(width=0.97, height=0.97))
  )
  
  pdf(out_path,width=8, height=6)
  print(out_plot)
  dev.off()
  
  return(out_plot)
}

#------ country_summary  -------------
#Histogram of number of sites x country x EFA_type
country_summary <- function(in_tab, out_path) {
  in_tab_nodupli <- in_tab[
    !duplicated(in_tab, by=c("Point_db", "UID_Mathis", "eftype_ref")) &
      !is.na(efvol_ref) | !is.na(efper_ref),
  ]
  
  in_tab_nodupli[, countryN := .N, by=Country]
  
  in_tab_nodupli[is.na(eftype_ref), eftype_ref := 'None specified']
  
  country_hist <- ggplot(in_tab_nodupli, aes(x=reorder(Country, countryN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=after_stat(count), label=after_stat(count)), stat='count', geom="text", hjust=-.5) +
    coord_flip(clip='off') +
    scale_x_discrete('Country') + 
    scale_y_continuous('Number of EFAs') +
    scale_fill_brewer('Type of EFA', palette='Set2', na.value="grey") +
    theme_classic() +
    theme(#axis.title.y = element_text(vjust=-2),
      plot.margin = unit( c(0.1, 0.6, 0.1, 0.1), 'cm'),
      legend.position = c(0.7, 0.2))
  
  ggsave(out_path, plot=country_hist, width=4, height=4)
  
  return(country_hist)
}

#------ ecoregions_summary -------------
#Histogram of number of sites x country x EFA_type
ecoregions_summary <- function(in_tab) {
  in_tab[, fecN := .N, by=fec_cl_cmj]
  
  ecoregions_hist <- ggplot(in_tab, aes(x=reorder(fec_cl_cmj, fecN))) + 
    geom_histogram(aes(fill=eftype_ref), stat='count', geom='text') +
    stat_count(aes(y=after_stat(count), label=after_stat(count)), stat='count', geom="text", hjust=-.5) +
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
map_ef <- function(in_efp, out_path) {
  in_efp <- unique(in_efp[!is.na(efvol_ref) | !is.na(efper_ref)], 
                   by=c('POINT_X', 'POINT_Y'))
  
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
    scale_size_continuous(name=str_wrap('Number of e-flow assessment sites within 1000 km', 25),
                          range=c(1, 6), breaks=c(1, 10, 100, 436)) +
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
          legend.background = element_blank()
    )
  
  ggsave(out_path, dpi=600, width=6, height=4, units='in')
  print(studiesmap)
  dev.off()
  
  return(studiesmap)
  
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
  
  #Bin records by upstream drainage area
  binlabels <- label_manualbins(binarg=c(100, 1000, 10000, 30000, 10000000),
                                minval=min(in_eftab$up_area_skm_15s)) %>%
    data.table(bin_label = .,
               bin=seq(1,5))
  
  #Merge points with records
  eftab_modjoin_preprocessed <- merge(in_eftab,
                                      efp_mod_format,
                                      by=c("UID_Mathis", "Point_db"),
                                      allow.cartesian=TRUE) %>%  
    bin_dt(binvar = "up_area_skm_15s",
           binfunc='manual',
           binarg = c(100, 1000, 10000, 30000, 10000000)) %>%
    merge(binlabels, by='bin') %>%
    .[(Comment_mathis != 'Not well represented in HydroRIVERS') |  #Remove those whose location is unsure
        (is.na(Comment_mathis)),] %>%
    .[!is.na(mar_ref), #Convert untis for reference mean annual runoff
      `:=`(mar_ref = mar_ref/31.5576, 
           mar_unit = "m3 s-1")] %>% 
    merge(.[, unique(.SD),  .SDcols=c('Point_db', 'UID_Mathis', 'Country')][ #Compute number of unique sites per country
      , list(Ncountry=.N), by='Country'],
      by='Country') %>%
    .[, Country_labelo10 := fifelse(Ncountry>=10, Country, 'Other')] %>% #Label all countries with less than 10 sites as "Other" for subsequent displays
    .[res=="0.5 arc-deg" & var=='qtot', value := value/1000] %>%  #Correctly scale flow values for isimp2b
    .[is.na(eftype_ref), eftype_ref := 'None specified'] %>%
    merge(dcast(.[eftype_format == 'maf',],  #Create new column for MAF rather than as a record
                EFUID+gcm+ghm+var~eftype_format,
                value.var = 'value'),
          by=c('EFUID', 'gcm', 'ghm', 'var')
    ) %>%
    .[!(eftype_format %in% c("maf", "mmf")),] %>% #Remove maf records
    .[ef_unit == "10^6m3 y-1", `:=`( #Standardize units for reference e-flows
      ef_unit = "m3 s-1",
      efvol_ref = efvol_ref/31.5576
    )]  %>%
    .[eftype == "smakthinef", eftype:="smakhtinef"] %>% #Correct typo in Smakhtin
    setnames(c('value', 'maf'), 
             c('efvol_mod', 'mar_mod')) %>%
    .[, SHAPE := NULL]
  
  #The highest e-flows (Smakhtin_a_ for the Amazon yielded NAs because the values 
  #were outside of the possible range for the raster.
  #Compute them based on the ratio in EF value of Smakhtin_b
  #between the actual pourpoint and the downstream-most location with data for the missing layer
  eftab_modjoin_preprocessed[EFUID==427 & eftype_format=='smakthin_a' &
                               run == 'h08_hadgem2-es_picontrol_1860soc' & 
                               var=='qtot',
                             efvol_mod := (144014.784/132200.856)*170483.088
  ]
  
  eftab_modjoin_preprocessed[EFUID==427 & eftype_format=='smakthin_a' &
                               run == 'h08_miroc5_picontrol_1860soc'& 
                               var=='qtot',
                             efvol_mod := (141730.848/129882.280)*169503.664
  ]
  
  #----- Match Smakhtin's estimates to reference vols based on ref EMC ---------
  eftab_mod_smakthin <- dcast(
    eftab_modjoin_preprocessed[ecpresent_ref_simplemin %in% c("A", "B", "C", "D", "E") &
                                 (eftype == "smakhtinef") & (var == 'qtot')],
    EFUID+UID_Mathis+Point_db+ecpresent_ref_formatted+ecpresent_ref_simplemin+
      ecpresent_ref_simplemax+run+var~eftype_format,
    value.var='efvol_mod')
  
  eftab_mod_smakthin[, smakhtin_emcmin := fcase(
    ecpresent_ref_simplemin == 'A', smakthin_a,
    ecpresent_ref_simplemin == 'B', smakthin_b,
    ecpresent_ref_simplemin == 'C', smakthin_c,
    ecpresent_ref_simplemin == 'D', smakthin_d,
    ecpresent_ref_simplemin == 'E', smakthin_d
  )]
  
  eftab_mod_smakthin[, smakhtin_emcmax := fcase(
    ecpresent_ref_simplemax == 'A', smakthin_a,
    ecpresent_ref_simplemax == 'B', smakthin_b,
    ecpresent_ref_simplemax == 'C', smakthin_c,
    ecpresent_ref_simplemax == 'D', smakthin_d,
    ecpresent_ref_simplemax == 'E', smakthin_d
  )]
  
  eftab_mod_smakthin[, smakhtin_emcavg := (smakhtin_emcmin+smakhtin_emcmax)/2]
  
  
  eftab_smakhtin <- merge(eftab_modjoin_preprocessed, 
                          eftab_mod_smakthin[, c("EFUID","UID_Mathis", "Point_db", 
                                                 "ecpresent_ref_formatted", "run", "var",
                                                 "smakhtin_emcavg"), 
                                             with=F],
                          by=c("EFUID","UID_Mathis", "Point_db", 
                               "ecpresent_ref_formatted", "run", "var"),
                          all.x=F) %>%
    .[, c( "eftype", "eftype_format", "emc", "efvol_mod"):=NULL] %>%
    melt(id.vars=names(.)[(names(.) != "smakhtin_emcavg")]) %>%
    setnames(c("variable", "value"), c("eftype_format", "efvol_mod")) 
  
  
  eftab_modjoin <- rbind(eftab_smakhtin, eftab_modjoin_preprocessed, 
                         use.names=T, fill=T) %>%
    .[, efper_mod := fifelse(mar_mod >0, 100*efvol_mod/mar_mod, 0)] #Compute percentage EF
  
  #----- Order factors ----------------------------------------------------------
  
  #Order eftype_format
  eftab_modjoin[, eftype_format := factor(
    eftype_format,
    levels=c('tennant', 'tessmann', 'q90q50', 'vmf',
             'smakthin_a',  'smakthin_b',  'smakthin_c',
             'smakthin_d',  'smakhtin_emcavg'),
    labels = c('Tennant', 'Tessmann', 'Q90_Q50', 'VMF',
               'Smakhtin - class A', 'Smakhtin - class B',
               'Smakhtin - class C', 'Smakhtin - class D',
               'Smakhtin - matched class')
  )]
  
  #Order Country_labelo10
  eftab_modjoin[, Country_labelo10 := factor(
    Country_labelo10,
    levels=c('Australia', 'Brazil', 'Canada', 'France',
             'Mexico', 'South Africa', 'Other')
  )]
  
  return(eftab_modjoin)
}

#------ compare_hydrology ------------------------------------------
# in_efp_efmod_join <- tar_read(efp_efmod_join)
# in_clz_labels <- tar_read(clz_labels)
# outdir = file.path(resdir, 'tables', 'hydrological_comparison')

compare_hydrology <- function(in_efp_efmod_join, in_clz_labels,
                              outdir = file.path(resdir, 'tables', 'hydrological_comparison')) {
  
  #Format data
  in_efp_efmod_join[is.na(mar_ref) &
                      (!is.na(efvol_ref) | !is.na(efper_ref)), 
                    nrow(unique(.SD)), 
                    .SDcols=c("UID_Mathis", "Point_db")] #164 sites without MAR value (114 without MAR but with EF)
  
  eftab_maf <- in_efp_efmod_join[!is.na(mar_ref),] %>% #remove sites without MAR value
    unique(by=c('POINT_X', 'POINT_Y', 'run', 'var', 'emc')) %>% #Keep unique sites
    .[gcm != 'cru-era',] #Remove high-res pcr-globwb
  
  eftab_maf[
    , n_sites := .SD[!duplicated(.SD, by=c('Point_db', 'UID_Mathis')) 
                     & !is.na(mar_ref), 
                     .N]
    , by='Country']
  
  #Compute percentage error
  eftab_maf[, `:=` (APE = 100*abs((mar_mod-mar_ref)/mar_ref),
                    PE = 100*((mar_mod-mar_ref)/mar_ref)
  )]
  
  eftab_maf[, DAo3000km2 := fifelse(up_area_skm_15s >= 3000, 'o3000', 'u3000')]
  
  
  #prepare data to add plus-minus one order of magnitude ribbons on plots
  magnitude_lines_10_dis <-seq(log10(eftab_maf[var=='dis', min(mar_ref, na.rm=T)/2]),
                               log10(eftab_maf[var=='dis', 2*max(mar_ref, na.rm=T)]),
                               0.1) %>%
    data.table(
      x = 10^.,
      y_low = 10*(10^.),
      y_high = (10^.)/10
    )
  
  magnitude_lines_10_qtot <-seq(log10(eftab_maf[var=='qtot', min(mar_ref, na.rm=T)/2]),
                                log10(eftab_maf[var=='qtot', 2*max(mar_ref, na.rm=T)]),
                                0.1) %>%
    data.table(
      x = 10^.,
      y_low = 10*(10^.),
      y_high = (10^.)/10
    )
  
  #--------------------------Compare reference and model MAR for raw discharge -----
  plot_obspred_discharge <- ggplot(eftab_maf[var=='dis',]) +
    geom_ribbon(data=magnitude_lines_10_dis, aes(x, ymin=y_low, ymax=y_high), 
                alpha=1/3, fill='#6CDAE7') +
    geom_point(aes(x=mar_ref, y=mar_mod), alpha=1/3) +
    scale_x_log10(name= expression('Mean annual flow estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_log10(name= expression('Mean annual flow estimated by global model'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    geom_abline () +
    coord_fixed() +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  plot_APEarea_discharge <- ggplot(
    eftab_maf[var=='dis',], aes(x=up_area_skm_15s, y=APE)) +
    geom_point() +
    scale_x_log10(name= 'Upstream drainage area'~(km^2),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name= "Absolute error (%)", 
                  breaks=c(10,100,1000,1000000),
                  labels = trans_format("log10", math_format(10^.x))) +
    #geom_hline(yintercept=100) +
    geom_vline(xintercept=c(3000), linetype='dashed') +
    geom_smooth(size=1, se=T, span=1, method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic()+ 
    theme(panel.grid.major.y = element_line())
  
  #--------------------------Compare reference and model MAR for downscaled discharge -----
  plot_obspred_downscaledqtot <-   ggplot(eftab_maf[var=='qtot',]) +
    geom_ribbon(data=magnitude_lines_10_qtot, aes(x, ymin=y_low, ymax=y_high), 
                alpha=1/3, fill='#6CDAE7') +
    geom_point(aes(x=mar_ref, y=mar_mod), alpha=1/3) +
    scale_x_log10(name= expression('Mean annual flow estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), expand=c(0,0),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name= expression('Mean annual flow estimated by global model'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), expand=c(0,0),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_abline () +
    geom_smooth(aes(x=mar_ref, y=mar_mod), se=FALSE, method='lm') +
    coord_fixed() +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major = element_line())
  
  #APE~Drainage area
  plot_APEarea_downscaledqtot <- ggplot(
    eftab_maf[var=='qtot',], 
    aes(x=up_area_skm_15s, y=APE)) +
    geom_point(alpha=1/5) + #aes(color=Country_labelo10)
    scale_x_log10(name= 'Upstream drainage area'~(km^2),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name= "Absolute percentage error (%)", 
                  breaks=c(1,10,100,1000,10000),
                  labels = trans_format("log10", math_format(10^.x))) +
    #scale_color_discrete(name = 'Country') +
    geom_vline(xintercept=c(3000), linetype='dashed') +
    #geom_quantile(quantiles=c(0.1, 0.5, 0.9), size=1) +
    geom_smooth(size=1, se=T, span=1, method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  
  #APE~Aridity Index
  plot_APEaridity_downscaledqtot <- ggplot(
    eftab_maf[var=='qtot' & gcm != 'cru-era',], 
    aes(x=ari_ix_uav, y=APE)) +
    geom_point(alpha=1/5, aes(color=Country_labelo10)) +
    scale_x_continuous(name= 'Global aridity index') +
    scale_y_log10(name= "Absolute percentage error (%)",
                  breaks=c(1,10,100,1000,10000),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_discrete(name = 'Country') +
    #geom_vline(xintercept=c(100, 2500), linetype='dashed') +
    geom_smooth(size=1, se=FALSE, span=1, color='black', method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  #APE~urb_pc_use
  plot_APEurban_downscaledqtot <- ggplot(
    eftab_maf[var=='qtot' & gcm != 'cru-era',], 
    aes(x=urb_pc_use, y=APE)) +
    geom_point(alpha=1/5, aes(color=Country_labelo10)) +
    scale_x_continuous(name= 'Upstream urban cover (%)', breaks=c(0,25,50,75,100)) +
    scale_y_log10(name= "Absolute error (%)", breaks=c(1,10,100,1000,10000),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_discrete(name = 'Country') +
    #geom_vline(xintercept=c(100, 2500), linetype='dashed') +
    #geom_smooth(size=1, se=TRUE, span=1) +
    # geom_smooth(size=1, se=FALSE, span=1, color='black', method='gam',
    #             formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    geom_quantile(quantiles=c(0.1, 0.5, 0.9), color='black') +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  #APE~crp_pc_use
  plot_APEcrop_downscaledqtot <- ggplot(
    eftab_maf[var=='qtot' & gcm != 'cru-era',], 
    aes(x=crp_pc_use, y=APE)) +
    geom_point(alpha=1/5, aes(color=Country_labelo10)) +
    scale_x_continuous(name= 'Upstream crop cover (%)', breaks=c(0,25,50,75,100)) +
    scale_y_log10(name= "Absolute error (%)", breaks=c(1,10,100,1000,10000),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_discrete(name = 'Country') +
    #geom_vline(xintercept=c(100, 2500), linetype='dashed') +
    #geom_smooth(width=1, se=TRUE, span=1) +
    geom_smooth(size=1, se=FALSE, span=1, color='black', method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 4)) +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  #APE~dor_pc_pva
  plot_APEreservoir_downscaledqtot <- ggplot(
    eftab_maf[var=='qtot' & gcm != 'cru-era',], 
    aes(x=dor_pc_pva, y=APE)) +
    geom_point(alpha=1/5, aes(color=Country_labelo10)) +
    scale_x_sqrt(name= 'Degree of regulation', breaks=c(0,10,100,1000)) +
    scale_y_log10(name= "Absolute error (%)", breaks=c(1,10,100,1000,10000),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_discrete(name = 'Country') +
    #geom_vline(xintercept=c(100, 2500), linetype='dashed') +
    #geom_smooth(width=1, se=TRUE, span=1) +
    geom_smooth(size=1, se=FALSE, span=1, color='black', method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 4)) +
    facet_grid(c("gcm", "ghm"), 
               labeller = "label_both") + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  #--------------------------Compare reference and model MAR for original and downscaled discharge for largest drainage basins -----
  plot_large_rivers_APE_dis_vs_qtot <- ggplot(
    dcast(
      unique(eftab_maf[gcm == 'gfdl-esm2m' & ghm == 'watergap2-2c' & 
                         up_area_skm_15s >= 30000,], by=c('EFUID', "var")),
      EFUID + Point_db + UID_Mathis + up_area_skm_15s + Country_labelo10 ~ var,
      value.var = 'APE'
    ),
    aes(x=dis, y=qtot)) +
    geom_point(aes(color=Country_labelo10), alpha=1) +
    scale_x_log10(name= expression('Absolute error at 30 arc-min (%)'),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(name= expression('Absolute error at 0.25 arc-min (%)'),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_abline () +
    #geom_smooth(aes(color=Country_labelo10), se=F, method='lm') +
    coord_fixed() +
    #facet_wrap(~var) + 
    theme_classic() + 
    theme(panel.grid.major.y = element_line())
  
  #-------------------------- Compute statistics -------------------------------
  qstats_all <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm')
  ] %>%
    .[, N := eftab_maf[, nrow(unique(.SD[,c('Point_db', 'UID_Mathis')]))]]
  
  #Compute average difference between worst- and best-performing GCM models per GHM
  avg_range_gcms <- qstats_all[, lapply(.SD, function(x) {max(x)-min(x)}),
                               .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
                               by=c('var', 'ghm')][
                                 , lapply(.SD, mean),
                                 .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
                                 by='var']
  
  #Compute average difference between worst- and best-performing GhM models per GCM
  avg_range_ghms <- qstats_all[, lapply(.SD, function(x) {max(x)-min(x)}),
                               .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
                               by=c('var', 'gcm')][
                                 , lapply(.SD, mean),
                                 .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
                                 by='var']
  
  range_all <- qstats_all[, lapply(.SD, function(x) {max(x)-min(x)}),
                          .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
                          by=c('var')]
  
  #COmp
  qstats_pixelthresh <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'DAo3000km2')
  ] %>%
    .[gcm == 'gfdl-esm2m' & ghm == 'watergap2-2c']  %>%
    merge(eftab_maf[, eftab_maf[, nrow(unique(.SD[,c('Point_db', 'UID_Mathis')])),
                                by='DAo3000km2']
    ],
    by='DAo3000km2') %>%
    setnames('V1', 'N')
  # %>%
  #   .[ ,
  #      lapply(.SD, mean), 
  #      .SDcols = c('smape', 'pbias', 'rsq', 'rsq_nooutliers'),
  #      by=c('var', 'DAo3000km2')] %>%
  
  
  qstats_o30000 <- eftab_maf[up_area_skm_15s >= 30000
                             ,
                             getqstats(
                               dt = .SD,
                               actual = 'mar_ref',
                               predicted =  'mar_mod',
                               rstudthresh= 3,
                               log = T
                             ), 
                             by=c('var', 'gcm', 'ghm')
  ] %>%
    .[, N := eftab_maf[up_area_skm_15s >= 30000,
                       nrow(unique(.SD[,c('Point_db', 'UID_Mathis')]))]]%>%
    .[gcm == 'gfdl-esm2m' & ghm == 'watergap2-2c'] 
  
  
  qstats_da <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'bin_label')
  ] %>%
    merge(eftab_maf[!duplicated(eftab_maf[,c("Point_db", "UID_Mathis")]) & !is.na(mar_ref), 
                    .N, by='bin_label'],
          by='bin_label') %>%
    setnames('bin_label', "Drainage area")
  
  qstats_clz <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'clz_cl_cmj')
  ]  %>%
    merge(in_clz_labels, by.x='clz_cl_cmj', by.y='GEnZ_ID') %>%
    merge(eftab_maf[!duplicated(eftab_maf[,c("Point_db", "UID_Mathis")]) & !is.na(mar_ref), 
                    .N, by='clz_cl_cmj'],
          by='clz_cl_cmj') %>%
    .[, clz_label := paste0(GEnZ_Name, ", n=", N)] 
  
  
  qstats_country <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T,
      mae_decomp = TRUE
    ), 
    by=c('var', 'gcm', 'ghm', 'Country')
  ] %>%
    merge(eftab_maf[!duplicated(eftab_maf[,c("Point_db", "UID_Mathis")]) & !is.na(mar_ref), 
                    .N, by='Country'],
          by='Country') %>%
    .[, Country_label := paste0(Country, ", n=", N)]
  
  

  qstats_u10 <- eftab_maf[n_sites < 10
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T,
      mae_decomp = TRUE
    ), 
    by=c('var', 'gcm', 'ghm')
  ]
  
  #--------------- Plot statistics ---------------------------------------------
  qstats_dis_DA_plotformat <- melt(qstats_da[var=='dis', 
                                             c('gcm', 'ghm',  'Drainage area',
                                               'smape', 'pbias', 'rsq_nooutliers'), 
                                             with=F],
                                   id.vars=c('gcm', 'ghm', 'Drainage area'))
  
  qstats_dis_DA_plot <- ggplot(qstats_dis_DA_plotformat,
                               aes(x=ghm, y=value, color=gcm)) + 
    geom_point() +
    facet_grid(variable~`Drainage area`, scales='free') +
    theme_classic() +      theme(panel.grid.major = element_line())
  
  
  qstats_qtot_all_plotformat <- melt(
    qstats_all[var=='qtot', 
               c('gcm', 'ghm', 'var',
                 'smape', 'pbias', 'rsq_nooutliers'), with=F],
    id.vars=c('gcm', 'ghm', 'var')) %>%
    .[, `:=`(max_value = max(value),
             min_value = min(value)),
      by=c('ghm', 'variable')]
  levels(qstats_qtot_all_plotformat$variable) <- c('sMAPE',
                                                   expression(`%`~"Bias"),
                                                   expression(R^2))
  
  qstats_qtot_all_plot <- ggplot(qstats_qtot_all_plotformat,
                                 aes(x=ghm, y=value)) + 
    geom_segment(aes(xend=ghm, y=min_value, yend=max_value)) +
    geom_point(aes(color=gcm), size=3) +
    facet_grid(variable~var, scales='free',
               labeller = label_parsed) +
    scale_y_continuous(name='Value') +
    theme_classic() +     
    theme(panel.grid.major = element_line(),
          panel.border = element_rect(fill=NA),
          strip.text.x = element_blank())
  
  qstats_qtot_DA_plotformat <- melt(qstats_da[var=='qtot', 
                                              c('gcm', 'ghm',  'Drainage area',
                                                'smape', 'pbias', 'rsq_nooutliers'), 
                                              with=F],
                                    id.vars=c('gcm', 'ghm', 'Drainage area'))
  
  qstats_qtot_DA_plot <- ggplot(qstats_qtot_DA_plotformat,
                                aes(x=`Drainage area`, y=value, color=gcm)) + 
    geom_point() +
    facet_grid(variable~ghm, scales='free') +
    theme_classic() + 
    theme(panel.grid.major.y = element_line(),
          axis.text.x = element_text(angle=30, hjust=1))
  
  qstats_qtot_country_plotformat <- melt(
    qstats_country[var=='qtot' & N>10 & gcm != 'cru-era',
                   c('gcm', 'ghm', 'var', 'Country_label',
                     'smape', 'pbias', 'rsq_nooutliers'), with=F],
    id.vars=c('gcm', 'ghm', 'var', 'Country_label')) %>%
    .[, `:=`(max_value = max(value),
             min_value = min(value)),
      by=c('ghm', 'variable', 'Country_label')]
  
  qstats_qtot_country_plot <- ggplot(qstats_qtot_country_plotformat,
                                     aes(x=ghm, y=value, color=gcm)) + 
    geom_segment(aes(xend=ghm, y=min_value, yend=max_value)) +
    geom_point() +
    facet_grid(variable~Country_label, scales='free') +
    theme_classic() + 
    theme(panel.grid.major.y = element_line(),
          axis.text.x = element_text(angle=30, hjust=1))
  
  
  
  #---- Select for each country the most performant combination of gcm and ghm -----
  #Create a table for display
  qstats_country_melt <- melt(qstats_country[var=='qtot' & N >= 10,],
                              id.vars=c('Country', 'var', 'gcm', 'ghm'),
                              measure.vars=c('smape', 'pbias', 'rsq')
  )

  qstats_other <- melt(qstats_country[var=='qtot' & N >= 10,],
                       id.vars=c('Country', 'var', 'gcm', 'ghm'),
                       measure.vars=c('smape', 'pbias', 'rsq')
  )

  best_model_comb <- rbind(
    qstats_country_melt[variable %in% c('smape', 'pbias'),] %>%
      .[order(value), .SD[1,], by=c('Country', 'variable')],
    qstats_country_melt[variable=='rsq',] %>%
      .[order(-value), .SD[1,], by=c('Country', 'variable')]
  ) %>%
    .[, gcm_ghm := paste(gcm, ghm, sep=' + ')] %>%
    dcast(Country~variable, value.var='gcm_ghm')

  ######################### OLD METHOD ###########################
  # #For each country with n>10, 
  # #select the best combination of ghm and gcm
  # #accoording to each of three criteria (pbias, smape and R2)
  # #Compute weights for later averaging by country
  # multicriterion_hydroweights_o10 <- rbindlist(list(
  #   qstats_country[var=='qtot' & n_total >= 10, .SD[order(pbias),][1,], by=c('Country')],
  #   qstats_country[var=='qtot' & n_total >= 10, .SD[order(smape),][1,], by=c('Country')],
  #   qstats_country[var=='qtot' & n_total >= 10, .SD[order(-rsq_nooutliers),][1,], by=c('Country')])
  # ) %>%
  #   .[, list(model_weight = .N), by=c('Country', 'ghm', 'gcm')]
  # 
  # 
  # #Compute overall weight for countries with n<10
  # multicriterion_hydroweights_u10 <- qstats_country[n_total<10, 
  #                                                   unique(Country)] %>%
  #   cbind(
  #     multicriterion_hydroweights_o10[, sum(model_weight), by=c('ghm', 'gcm')][
  #       rep(1:.N, each=length(.)),]
  #   ) %>%
  #   setnames(c('.', 'V1'), c('Country', 'model_weight'))
  # 
  # #Produce weighted average model estimates
  # ensemble_mod <- merge(
  #   in_efp_efmod_join[var=='qtot', 
  #                     c('Country', 'EFUID', 'UID_Mathis', 'Point_db',
  #                       'ghm', 'gcm',
  #                       'mar_mod', 'efvol_mod', 'efper_mod', 'eftype_format')],
  #   rbind(multicriterion_hydroweights_o10,
  #         multicriterion_hydroweights_u10,
  #         use.names=T),
  #   by=c("Country", "ghm", "gcm"),
  #   all.x=F
  # ) %>%
  #   .[, list(mar_mod = weighted.mean(mar_mod, model_weight),
  #            efvol_mod = weighted.mean(efvol_mod, model_weight),
  #            efper_mod = weighted.mean(efper_mod, model_weight)),
  #     by = c("EFUID", "UID_Mathis", "Point_db", "eftype_format")]
  ######################### UPDATE WITH Sanderson et al. 2017 method ###########################

  #Skill weighting based on MAE decomposition 
  #--- individually for countries with more than 10 records
  qstats_country[
    var=='qtot' & n_total >= 10, 
    `:=`
    (w_bias = exp(-(mae_log10_bias/(0.4*min(mae_log10_bias)))^2),
      w_proportionality_error= exp(-(mae_log10_proportionality_error/(0.4*min(mae_log10_proportionality_error)))^2),
      w_unsystematic_error = exp(-(mae_log10_unsystematic_error/(0.4*min(mae_log10_unsystematic_error)))^2)
    ),
    by=Country]
  
  qstats_country[
    var=='qtot' & n_total >= 10,
    model_weight := (
      (w_bias*1/sum(w_bias)) +
        (w_proportionality_error*1/sum(w_proportionality_error)) +
        (w_unsystematic_error*1/sum(w_unsystematic_error))
    )/3,
    by=Country]
  
  #Altogether for countries with less than 10 records
  qstats_u10[
    var=='qtot', 
    `:=`
    (w_bias = exp(-(mae_log10_bias/(0.4*min(mae_log10_bias)))^2),
      w_proportionality_error= exp(-(mae_log10_proportionality_error/(0.4*min(mae_log10_proportionality_error)))^2),
      w_unsystematic_error = exp(-(mae_log10_unsystematic_error/(0.4*min(mae_log10_unsystematic_error)))^2)
    )]
  
  qstats_u10[
    var=='qtot',
    model_weight := (
      (w_bias*1/sum(w_bias)) +
        (w_proportionality_error*1/sum(w_proportionality_error)) +
        (w_unsystematic_error*1/sum(w_unsystematic_error))
    )/3]
  
  
  qstats_u10_format <- qstats_country[var=='qtot' & n_total<10, unique(Country)] %>%
    cbind(
      qstats_u10[var=='qtot', model_weight, by=c('ghm', 'gcm')][
        rep(1:.N, each=length(.)),]
    ) %>% setnames('.', 'Country') %>%
    .[order(Country),]
  
  multicriterion_hydroweights <- rbind(
    qstats_country[var=='qtot' & n_total >= 10, .(Country, ghm, gcm, model_weight)],
    qstats_u10_format)

  #Produce weighted average model estimates
  ensemble_mod <- merge(
    in_efp_efmod_join[var=='qtot', 
                      c('Country', 'EFUID', 'UID_Mathis', 'Point_db',
                        'ghm', 'gcm',
                        'mar_mod', 'efvol_mod', 'efper_mod', 'eftype_format')],
    multicriterion_hydroweights,
    by=c("Country", "ghm", "gcm"),
    all.x=F
  ) %>%
    .[, list(mar_mod = weighted.mean(mar_mod, model_weight),
             efvol_mod = weighted.mean(efvol_mod, model_weight),
             efper_mod = weighted.mean(efper_mod, model_weight)),
      by = c("EFUID", "UID_Mathis", "Point_db", "eftype_format")]
  
  multicriterion_hydroweights[, gcm_ghm := paste(gcm,ghm,sep='+')]
  
  ########################################################################################
  
  keep_cols <- names(in_efp_efmod_join)[
    !(names(in_efp_efmod_join) %in% 
        c('ghm', 'gcm', 'run', 'mar_mod', 'efvol_mod', 
          'efper_mod'))]
  
  eftab_ensemblemod <- merge(
    in_efp_efmod_join[var=='qtot',][
      !duplicated(in_efp_efmod_join[var=='qtot', keep_cols, with=F]),
    ][, keep_cols, with=F]
    ,
    ensemble_mod,
    by=c("EFUID", "UID_Mathis", "Point_db", "eftype_format")
  )
  
  
  #---- Compute statistics -------------
  eftab_ensemblemod_maf <- eftab_ensemblemod[!is.na(mar_ref),] %>% #remove sites without MAR value
    unique(by=c('POINT_X', 'POINT_Y')) %>% #Keep unique sites
    .[, `:=` (APE = 100*abs((mar_mod-mar_ref)/mar_ref),
              PE = 100*((mar_mod-mar_ref)/mar_ref)
    )]
  #Compute statistics
  qstats_ensemble_all <- getqstats(
    dt = eftab_ensemblemod_maf,
    actual = 'mar_ref',
    predicted =  'mar_mod',
    rstudthresh= 3,
    log = T
  )
  
  qstats_ensemble_da <- eftab_ensemblemod_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('bin_label')
  ] %>%
    setnames('bin_label', "Drainage area")
  
  qstats_ensemble_clz <- eftab_ensemblemod_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('clz_cl_cmj')
  ]  %>%
    merge(in_clz_labels, by.x='clz_cl_cmj', by.y='GEnZ_ID') %>%
    .[, clz_label := paste0(GEnZ_Name, ", n=", n_total)] 
  
  
  qstats_ensemble_country <- eftab_ensemblemod_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('Country_labelo10')
  ]
  
  
  #---- Plot resulting model performance-----
  plot_obspred_ensemble <- ggplot(
    eftab_ensemblemod_maf, aes(x=mar_ref, y=mar_mod)) +
    geom_point(aes(color=Country_labelo10), alpha=1/2) +
    scale_x_log10(name= expression('Mean annual flow estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_log10(name= expression('Mean annual flow estimated by global model'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    geom_abline () +
    geom_smooth(se=FALSE,span=1, color='black') +
    theme_classic() +
    facet_wrap(~Country_labelo10) +
    theme(legend.position = 'none')
  
  #PE~reference MAR
  plot_PEmar_downscaledqtot <- ggplot(
    eftab_ensemblemod_maf, 
    aes(x=mar_ref, y=PE, 
        color=Country_labelo10)) +
    geom_hline(yintercept=0) +
    geom_point(alpha=1/2) +
    scale_x_log10(name= expression('Mean annual flow estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000), 
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_continuous(name= "Percentage error (%)") +
    geom_smooth(size=1, se=F, span=1, color='black', alpha=1/2)  +
    facet_wrap(~Country_labelo10, scales='free_y') +
    theme_classic() + 
    theme(panel.grid.major.y = element_line(),
          legend.position = 'none')
  
  
  #PE~Drainage area
  plot_PEarea_downscaledqtot <- ggplot(
    eftab_ensemblemod_maf, 
    aes(x=up_area_skm_15s, y=PE, 
        color=Country_labelo10)) +
    geom_hline(yintercept=0) +
    geom_point() +
    scale_x_log10(name= 'Upstream drainage area'~(km^2)) +
    scale_y_continuous(name= "Percentage error (%)") +
    #geom_smooth(size=1, se=T, span=1, color='black', alpha=1/2)  +
    geom_smooth(width=1, se=F,method='lm', color='black', alpha=1/2)  +
    facet_wrap(~Country_labelo10, scales='free_y') +
    theme_classic() + 
    theme(panel.grid.major.y = element_line(),
          legend.position = 'none')
  
  ####NOT A MATTER OF DRAINAGE AREA
  
  
  #---- Compare model performance by country -----
  qstats_format_ensemble_vs_best <- eftab_maf[
    ,
    getqstats(
      dt = .SD,
      actual = 'mar_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('var', 'gcm', 'ghm', 'Country_labelo10')
  ] %>% 
    .[gcm == 'gfdl-esm2m' & ghm == 'watergap2-2c' & var=='qtot',] %>%
    rbind(qstats_ensemble_country, use.names=T, fill = T) %>%
    .[, approach_label := fifelse(is.na(gcm), 'ensemble', paste(gcm, ghm, sep=' + '))] %>%
    melt(id.vars=c('Country_labelo10', 'approach_label'), 
         measure.vars=c('smape', 'pbias', 'rsq')
    )
  
  qstats_format_ensemble_vs_best_segment <- dcast(
    qstats_format_ensemble_vs_best, 
    Country_labelo10+variable~approach_label,
    value.var = 'value'
  ) %>%
    .[variable %in% c('smape', 'pbias'), 
      outcome := fifelse((abs(ensemble) - abs(`gfdl-esm2m + watergap2-2c`)) < 0,
                         'improvement', 'deterioration')] %>%
    .[variable == 'rsq', 
      outcome := fifelse((ensemble-`gfdl-esm2m + watergap2-2c`) > 0,
                         'improvement', 'deterioration')] 
  
  plot_ensemble_vs_best <- ggplot(qstats_format_ensemble_vs_best,
                                  aes(x=Country_labelo10)) + 
    geom_segment(data=qstats_format_ensemble_vs_best_segment,
                 aes(xend= Country_labelo10, 
                     y=ensemble, yend=`gfdl-esm2m + watergap2-2c`,
                     color=outcome), size=2) + #x=x, ,
    geom_point(aes(y=value, shape=approach_label), 
               size=4, alpha=1/2) + #, height=1, width=0) +
    scale_shape_discrete(name = 'Model') + 
    scale_color_discrete(name = 'Outcome from using ensemble model') +
    facet_wrap(~variable, scales='free') + 
    theme_classic() +
    theme(axis.title.x=element_blank())
  
  ####NOT A MATTER OF DRAINAGE AREA
  
  #--------------------------- Get mean and SD rank of model combinations across country ------
  qstats_all_ensemble_country_rbind <- rbind(
    eftab_maf[
      ,
      getqstats(
        dt = .SD,
        actual = 'mar_ref',
        predicted =  'mar_mod',
        rstudthresh= 3,
        log = T
      ), 
      by=c('var', 'gcm', 'ghm', 'Country_labelo10')
    ],
    qstats_ensemble_country[, `:=`(gcm='ensemble', ghm='ensemble', var='qtot')]
  )

  qstats_rank_all_vs_ensemble <- qstats_all_ensemble_country_rbind[
    gcm != 'cru-era',
    list(
      ghm,
      gcm,
      ranking_smape = frank(smape),
      ranking_pbias = frank(pbias),
      ranking_rsq = frank(-rsq)
    ),
    by=c('Country_labelo10')
  ][, list(mean_rank = mean(c(ranking_smape,
                              ranking_pbias,
                              ranking_rsq)),
           sd_rank = sd(c(ranking_smape,
                          ranking_pbias,
                          ranking_rsq))
  ),
  by=c('ghm', 'gcm')] 
  
  
  #--------------- Write out statistics ----------------------------------------
  if (!dir.exists(outdir)) {
    dir.create(outdir)  
  }
  
  fwrite(qstats_all, file.path(outdir, 
                               paste0('qstats_all',
                                      format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_da, file.path(outdir, 
                              paste0('qstats_da',
                                     format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_clz, file.path(outdir, 
                               paste0('qstats_clz',
                                      format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_country, file.path(outdir, 
                                   paste0('qstats_country',
                                          format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_ensemble_all, file.path(outdir, 
                                        paste0('qstats_ensemble_all',
                                               format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_ensemble_da, file.path(outdir, 
                                       paste0('qstats_ensemble_da',
                                              format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_ensemble_country, file.path(outdir, 
                                            paste0('qstats_ensemble_country',
                                                   format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(dcast(multicriterion_hydroweights, Country~gcm_ghm, value.var = 'model_weight'),
         file.path(outdir, 
                   paste0('multicriterion_hydroweights',
                          format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(best_model_comb, file.path(outdir, 
                                    paste0('best_model_comb',
                                           format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  fwrite(qstats_rank_all_vs_ensemble,
         file.path(outdir, 
                   paste0('qstats_rank_all_vs_ensemble',
                          format(Sys.Date(), '%Y%m%d'), '.csv')))
  
  #--------------- Return results ----------------------------------------------
  return(list(plot_obspred_discharge = plot_obspred_discharge, 
              plot_APEarea_discharge = plot_APEarea_discharge,
              plot_obspred_downscaledqtot = plot_obspred_downscaledqtot, 
              plot_APEarea_downscaledqtot = plot_APEarea_downscaledqtot,
              plot_APEaridity_downscaledqtot = plot_APEaridity_downscaledqtot,
              plot_APEurban_downscaledqtot = plot_APEurban_downscaledqtot,
              plot_APEcrop_downscaledqtot = plot_APEcrop_downscaledqtot,
              plot_APEreservoir_downscaledqtot = plot_APEreservoir_downscaledqtot,
              plot_obspred_ensemble = plot_obspred_ensemble,
              plot_PEarea_downscaledqtot = plot_PEarea_downscaledqtot,
              qstats_dis_DA_plot = qstats_dis_DA_plot,
              qstats_qtot_all_plot = qstats_qtot_all_plot,
              qstats_qtot_DA_plot = qstats_qtot_DA_plot,
              qstats_qtot_country_plot = qstats_qtot_country_plot,
              qstats_ensemble_vs_best = plot_ensemble_vs_best,
              table_all = qstats_all,
              table_all_avg = qstats_all[, lapply(.SD, mean), 
                                         .SDcols = c('smape', 'pbias', 'rsq',
                                                     'rsq_nooutliers'),
                                         by='var'],
              table_pixel_thresh = qstats_pixelthresh,
              table_clz = qstats_clz,
              table_da = qstats_da,
              table_country = qstats_country,
              table_range_gcms = avg_range_gcms,
              table_range_ghms = avg_range_ghms,
              table_range_all = range_all,
              table_ensemble_all = qstats_ensemble_all,
              table_ensemble_da = qstats_ensemble_da,
              table_ensemble_country = qstats_ensemble_country,
              eftab_ensemblemod = eftab_ensemblemod
  ))
}

#------ compare_EFestimate ------------------------------------------
# in_efp_efmod_join = tar_read(efp_efmod_join)
# in_eftab_ensemble = tar_read(hydrology_comparison)$eftab_ensemblemod
# outdir = file.path(resdir, 'tables', 'EFestimate_comparison')

compare_EFestimate <- function(in_efp_efmod_join,
                               in_eftab_ensemble,
                               outdir = file.path(resdir, 'tables', 'EFestimate_comparison')) {
  
  #Format data for ensemble estimates and all hydro models ---------------------
  format_eftab_forcomparison <- function(in_dt) {
    
    eftab <- in_dt %>% 
      .[efper_mod > 100, efper_mod := 100] %>% #9 cases 
      .[efper_ref <= 100,] %>% #27 sites are excluded.
      .[(efvol_ref > 0) | is.na(efvol_ref),]
    
    #Compute difference between model estimates and reference values
    eftab[!is.na(efper_ref) & efper_ref>0, 
          efper_pmae := abs(efper_mod-efper_ref)/efper_ref] 
    eftab[!is.na(efper_ref) & efper_ref>0, 
          efper_diff := efper_mod-efper_ref] 
    eftab[!is.na(efper_ref) & efper_ref>0, 
          efper_adiff := abs(efper_mod-efper_ref)] 
    eftab[!is.na(efvol_ref) & efvol_ref>0, 
          efvol_pmae := abs(efvol_mod-efvol_ref)/efvol_ref] 
    eftab[!is.na(efvol_ref) & efvol_ref>0, 
          efvol_diff := efvol_mod-efvol_ref]
    eftab[!is.na(efvol_ref) & efvol_ref>0, 
          efvol_adiff := abs(efvol_mod-efvol_ref)]
    eftab[!is.na(mar_ref) & mar_ref>0, 
          mar_diff := mar_mod-mar_ref] 
    eftab[!is.na(mar_ref) & mar_ref>0,
          mar_adiff := abs(mar_mod-mar_ref)] 
    eftab[!is.na(mar_ref) & mar_ref>0, 
          mar_aperdiff := 2*abs(mar_mod-mar_ref)/(mar_mod+mar_ref)]
    
    return(eftab)
  }
  
  
  formatted_eftab_all <- format_eftab_forcomparison(in_efp_efmod_join) %>%
    .[!(hydrotype_ref %in% c("Wet year", "Dry year", "Drought year")),]
  
  #Select best and worst scenarios based on PMAE
  eftab_best_all <- formatted_eftab_all[order(efvol_pmae),] %>%
    .[!(duplicated( #Remove duplicate sites (multiple scenarios, keeping the one with least EC difference when standard system)
      formatted_eftab_all[order(efvol_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'gcm', 'ghm', 'var', 'eftype_format'))),
    ] 
  
  eftab_bestper_all <- formatted_eftab_all[order(efper_pmae),][
    !(duplicated( #Remove duplicate sites (multiple scenarios, keeping the one with least EC difference when standard system)
      formatted_eftab_all[order(efper_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'gcm', 'ghm', 'var', 'eftype_format'))),
  ] 
  
  formatted_eftab_ensemble <- format_eftab_forcomparison(in_eftab_ensemble)
  
  #Select best and worst scenarios based on PMAE for the ensemble
  eftab_best_ensemble <- formatted_eftab_ensemble[order(efvol_pmae),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      formatted_eftab_ensemble[order(efvol_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'eftype_format'))),
  ] 
  
  eftab_bestper_ensemble <- formatted_eftab_ensemble[order(efper_pmae),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      formatted_eftab_ensemble[order(efper_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'eftype_format'))),
  ] 
  
  eftab_worst_ensemble <- formatted_eftab_ensemble[order(-efvol_pmae),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      formatted_eftab_ensemble[order(-efvol_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'eftype_format'))),
  ] 
  
  eftab_worstper_ensemble <- formatted_eftab_ensemble[order(-efper_pmae),][
    !(duplicated( #Remove duplicate sites (multiple scenarions, keeping the one with least EC difference when standard system)
      formatted_eftab_ensemble[order(-efper_pmae),], 
      by=c('POINT_X', 'POINT_Y', 'eftype_format'))),
  ] 
  
  
  #--------------- Check whether there is correspondence between EMC and % EF -----------
  #Test differences among classes
  compute_tukeyletters <- function(in_dt, var_x, var_y) {
    aov_out<- aov(as.formula(paste(var_y, '~', var_x)),
                  data=in_dt)    
    tukey_out <- TukeyHSD(aov_out)
    cld <- multcompView::multcompLetters4(aov_out, tukey_out)
    as.data.table(as.data.frame.list(cld[[1]]), keep.rownames = T)
  }
  
  ecref_tukeyletters <- lapply(
    c('Australia', 'France', 'South Africa'),
    function(c) {
      formatted_eftab_all[
        ecpresent_ref_formatted %in% c('A', 'A/B', 'B', 'B/C', 
                                       'C', 'C/D', 'D', 'D/E', 'E') &
          #Ncountry > 10 & 
          Country == c &
          !duplicated(EFUID), 
        compute_tukeyletters(in_dt = .SD, 
                             var_x = 'ecpresent_ref_formatted',
                             var_y='efper_ref')[, c(1,2)]
      ] %>%
        cbind(c)
    }
  ) %>% rbindlist %>%
    setnames(c('c', 'rn'), 
             c('Country', 'ecpresent_ref_formatted'))
  
  plot_efper_ecpresent_ref <- ggplot(
    formatted_eftab_all[
      ecpresent_ref_formatted %in% c('A', 'A/B', 'B', 'B/C', 
                                     'C', 'C/D', 'D', 'D/E', 'E') &
        Ncountry > 10 &
        !duplicated(EFUID),],
    aes(x=ecpresent_ref_formatted, y=efper_ref)) + 
    geom_boxplot() +
    scale_x_discrete(name = 'Present-day ecological condition class') +
    scale_y_continuous(
      name = 'Reference e-flow (% of mean annual flow)',
      limits=c(0,105), expand=c(0,0)) +
    stat_summary(fun.data = give.n, geom = "text", fun = median) +
    geom_text(data=ecref_tukeyletters, aes(y=100, label=Letters)) +
    facet_wrap(~Country) + 
    theme_classic()
  
  ####################### Analysis by volume ##########################################################
  #--------------------- Stats -----------------------------------------
  efvol_stats_all <- eftab_best_all[!is.na(efvol_ref) & 
                                      var=='qtot'& 
                                      gcm != 'cru-era' &
                                      #Ncountry > 10 &
                                      efvol_ref >0
                                    ,
                                    getqstats(
                                      dt = .SD,
                                      actual = 'efvol_ref',
                                      predicted =  'efvol_mod',
                                      rstudthresh= 3,
                                      log = T
                                    ), 
                                    by=c('var', 'gcm', 'ghm', 
                                         'eftype_format')
  ] 
  
  efvol_stats_country <- eftab_best_all[!is.na(efvol_ref) & 
                                          var=='qtot'& 
                                          gcm != 'cru-era' &
                                          Ncountry > 10
                                        ,
                                        getqstats(
                                          dt = .SD,
                                          actual = 'efvol_ref',
                                          predicted =  'efvol_mod',
                                          rstudthresh= 3,
                                          log = T
                                        ), 
                                        by=c('var', 'gcm', 'ghm', 
                                             'Country', 'eftype_format')
  ] 
  
  efvol_stats_eftype <- eftab_best_all[!is.na(efvol_ref) & 
                                         var=='qtot'& 
                                         gcm != 'cru-era' &
                                         Ncountry > 10
                                       ,
                                       getqstats(
                                         dt = .SD,
                                         actual = 'efvol_ref',
                                         predicted =  'efvol_mod',
                                         rstudthresh= 3,
                                         log = T
                                       ), 
                                       by=c('var', 'gcm', 'ghm', 
                                            'eftype_ref',
                                            'eftype_format')
  ] 
  
  #Compute average difference in e-flow volume estimate among all model combinations
  avg_range_all <- eftab_best_all[
    !is.na(efvol_ref) &
      var=='qtot' &
      gcm != 'cru-era' &
      efvol_mod > 0,
    list(
      efvol_mod_range = 100*(max(efvol_mod)-min(efvol_mod))/min(efvol_mod)
    ),
    by=c("EFUID")
  ]
                                   
  #Compute median percentage difference between highest and lowest e-flowvolume estimate across GHM-GCM combinations
  #for a given e-flow estimation method
  median_range_hydromods <- eftab_best_all[
    !is.na(efvol_ref) &
      var=='qtot' &
      gcm != 'cru-era' &
      efvol_mod > 0,
    list(
      efvol_mod_range = 100*(max(efvol_mod)-min(efvol_mod))/min(efvol_mod)
    ),
    by=c("EFUID", "eftype_format")
  ] %>%
    .[, median(efvol_mod_range),]
  
  #Compute median percentage difference between highest and lowest e-flow volume estimate across 
  #e-flow estimation methods for a given GHM-GCM combination
  median_range_efmods <- eftab_best_all[
    !is.na(efvol_ref) &
      var=='qtot' &
      gcm != 'cru-era' &
      efvol_mod > 0,
    list(
      efvol_mod_range = 100*(max(efvol_mod)-min(efvol_mod))/min(efvol_mod)
    ),
    by=c("EFUID", "run")
  ] %>%
    .[, median(efvol_mod_range),]
  
  #----- Only for ensemble estimates
  efvol_stats_ensemble_all <- eftab_best_ensemble[!is.na(efvol_ref) & 
                                                    var=='qtot'& 
                                                    #Ncountry > 10 &
                                                    efvol_ref >0
                                                  ,
                                                  getqstats(
                                                    dt = .SD,
                                                    actual = 'efvol_ref',
                                                    predicted =  'efvol_mod',
                                                    rstudthresh= 3,
                                                    log = T
                                                  ), 
                                                  by=c('eftype_format')
  ] 
  
  
  efvol_stats_ensemble_country_best <- eftab_best_ensemble[!is.na(efvol_ref) & 
                                                        var=='qtot'& 
                                                        Ncountry > 10
                                                      ,
                                                      getqstats(
                                                        dt = .SD,
                                                        actual = 'efvol_ref',
                                                        predicted =  'efvol_mod',
                                                        rstudthresh= 3,
                                                        log = T
                                                      ), 
                                                      by=c('Country', 'eftype_format')
  ] 
  
  efvol_stats_ensemble_country_best[efvol_stats_ensemble_country_best[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], by=Country]$V1,
    label := eftype_format]
  
  efvol_stats_ensemble_country_worst <- eftab_worst_ensemble[!is.na(efvol_ref) & 
                                                             var=='qtot'& 
                                                             Ncountry > 10
                                                           ,
                                                           getqstats(
                                                             dt = .SD,
                                                             actual = 'efvol_ref',
                                                             predicted =  'efvol_mod',
                                                             rstudthresh= 3,
                                                             log = T
                                                           ), 
                                                           by=c('Country', 'eftype_format')
  ] 
  
  efvol_stats_ensemble_country_worst[efvol_stats_ensemble_country_worst[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], by=Country]$V1,
    label := eftype_format]
  
  efvol_stats_ensemble_eftype <- eftab_best_ensemble[!is.na(efvol_ref) & 
                                                       var=='qtot'& 
                                                       Ncountry > 10 &
                                                       eftype_ref != 'Other'
                                                     ,
                                                     getqstats(
                                                       dt = .SD,
                                                       actual = 'efvol_ref',
                                                       predicted =  'efvol_mod',
                                                       rstudthresh= 3,
                                                       log = F
                                                     ), 
                                                     by=c('eftype_ref',
                                                          'eftype_format')
  ] 
  
  efvol_stats_ensemble_eftype[efvol_stats_ensemble_eftype[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], by=eftype_ref]$V1,
    label := eftype_format]
  
  efvol_stats_ensemble_eftypecountry <- eftab_best_ensemble[!is.na(efvol_ref) & 
                                                              var=='qtot'& 
                                                              Ncountry > 10 &
                                                              eftype_ref != 'Other'
                                                            ,
                                                            getqstats(
                                                              dt = .SD,
                                                              actual = 'efvol_ref',
                                                              predicted =  'efvol_mod',
                                                              rstudthresh= 3,
                                                              log = T
                                                            ), 
                                                            by=c('eftype_ref',
                                                                 'eftype_format',
                                                                 'Country')
  ] 
  
  efvol_stats_ensemble_eftypecountry[efvol_stats_ensemble_eftypecountry[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], 
    by=c('eftype_ref', 'Country')]$V1,
    label := eftype_format]
  
  #Compare the performance stats between using the global EF estimate as 
  #predictor vs. simply using the global MAR estimate as predictor
  efvol_stats_comparison_ensemble_comparison_marefvol<- eftab_best_ensemble[
    !is.na(efvol_ref) & 
      var=='qtot'& 
      #Ncountry > 10 &
      efvol_ref >0
    ,
    getqstats(
      dt = .SD,
      actual = 'efvol_ref',
      predicted =  'mar_mod',
      rstudthresh= 3,
      log = T
    ), 
    by=c('eftype_format')
  ]  %>%
    merge(efvol_stats_ensemble_all, by=c('eftype_format'),
          suffixes = c('_marmod', '_efmod'))
  
  #Assess whether ensemble performs better than others----------------------------
  efvol_stats_all_ensemble_country_rbind <- rbind(
    efvol_stats_country[, label:=NA],
    efvol_stats_ensemble_country_best[, `:=`(gcm='ensemble', ghm='ensemble', var='qtot')]
  )
  
  stats_all_vs_ensemble <- efvol_stats_all_ensemble_country_rbind[
    gcm != 'cru-era', 
    lapply(.SD, function(x) {mean(x)}),
    by=c('gcm', 'ghm', 'Country'), 
    .SDcols=c('smape', 'pbias', 'rsq', 'rsq_nooutliers')]
  
  rank_all_vs_ensemble <- efvol_stats_all_ensemble_country_rbind[
    gcm != 'cru-era',
    list(
      ghm,
      gcm,
      ranking_smape = frank(smape),
      ranking_pbias = frank(pbias),
      ranking_rsq = frank(-rsq)
    ),
    by=c('eftype_format', 'Country')
  ][, list(mean_rank_smape = mean(ranking_smape),
           sd_rank_smape = sd(ranking_smape),
           mean_rank_pbias = mean(ranking_pbias),
           sd_rank_pbias = sd(ranking_pbias),
           mean_rank_rsq = mean(ranking_rsq),
           sd_rank_rsq= sd(ranking_rsq)
           ),
    by=c('ghm', 'gcm')] 
  
  #Analyze relationship between e-flow diff and discharge diff ----------------------------
  ggplot(eftab_best_ensemble,
         aes(x=mar_adiff, 
             y=efvol_adiff)) + 
    geom_point(aes(color=Country_labelo10)) + 
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth(method='lm', color='black', se=F) +
    facet_grid(c("Country_labelo10", "eftype_format"),
               scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    #coord_fixed() +
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  
  ggplot(eftab_best_ensemble,
         aes(x=mar_adiff, 
             y=efvol_adiff)) + 
    geom_point()  +
    scale_x_log10() +
    scale_y_log10() +
    geom_smooth() +
    coord_fixed() +
    theme_classic()
  
  error_correlation_mar_evol_eftype_country <- eftab_best_ensemble[
    !is.na(efvol_diff) &
      !is.na(mar_diff) &
      Ncountry > 10
    ,
    getqstats(
      dt = .SD,
      actual = 'efvol_adiff',
      predicted =  'mar_adiff',
      rstudthresh= 3,
      log = T
    ), 
    by=c('eftype_format',
         'Country')
  ] 
  
  ggplot(error_correlation_mar_evol_eftype_country,
         aes(x=Country, y=rsq, color=eftype_format)) + 
    geom_point() + 
    theme_classic()
  
  error_correlation_mar_evol_all <- eftab_best_ensemble[!is.na(efvol_diff) &
                                                                  !is.na(mar_diff) &
                                                                  Ncountry > 10
                                                                ,
                                                                getqstats(
                                                                  dt = .SD,
                                                                  actual = 'efvol_adiff',
                                                                  predicted =  'mar_adiff',
                                                                  rstudthresh= 3,
                                                                  log = T
                                                                )
  ] 
  
  #----------- Compute % of over- and underestimation for each method and country -------
  efvol_error_over_and_underestimation <- eftab_best_ensemble[
    !is.na(efvol_ref) & 
      var=='qtot',
    list(
      efvol_over = .SD[efvol_diff >= 10, .N]/.N,
      efvol_under = .SD[efvol_diff <= -10, .N]/.N
    ),
    by=c('Country_labelo10', 'eftype_format')
  ]
  
  eftab_best_ensemble[
    !is.na(efvol_ref) & 
      var=='qtot',
    list(
      efvol_over = .SD[efvol_diff >= 0, .N]/.N,
      efvol_under = .SD[efvol_diff <= -0, .N]/.N
    ),
    by=c('eftype_format')
  ]
  
  eftab_best_ensemble[
    !is.na(efvol_ref) & 
      var=='qtot',
    .SD[efvol_diff > 0, .N]/.N]
  
  eftab_best_ensemble[
    !is.na(efvol_ref) & 
      var=='qtot',
    .SD[efvol_diff > -10, .N]/.N,
    by=c('eftype_format')
  ]
  
  
  #--------------------- Plots -----------------------------------------------------------------
  magnitude_lines_10 <-seq(log10(eftab_best_ensemble[var=='qtot', min(efvol_ref, na.rm=T)/2]),
                           log10(eftab_best_ensemble[var=='qtot', 2*max(efvol_ref, na.rm=T)]),
                           0.1) %>%
    data.table(
      x = 10^.,
      y_low = 10*(10^.),
      y_high = (10^.)/10
    )
  
  magnitude_lines_2 <-seq(log10(eftab_best_ensemble[var=='qtot', min(efvol_ref, na.rm=T)/2]),
                          log10(eftab_best_ensemble[var=='qtot', 2*max(efvol_ref, na.rm=T)]),
                          0.1) %>%
    data.table(
      x = 10^.,
      y_low = 2*(10^.),
      y_high = (10^.)/2
    )
  
  #All results----------------
  EFcompare_ensemble_best_vol_scatter <- ggplot(eftab_best_ensemble[var=='qtot',]) +
    geom_ribbon(data=magnitude_lines_10, aes(x=x, ymin=y_low, ymax=y_high), 
                alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2) +
    geom_point(aes(x=efvol_ref, y=efvol_mod), alpha=1/3) + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10(name= expression('Mean annual e-flows estimated by local e-flow assessment'~(m^3~s^-1)),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_log10(name= expression('Mean annual e-flows estimated by global models'~(m^3~s^-1)),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    geom_smooth(aes(x=efvol_ref, y=efvol_mod), span=1) +
    coord_fixed() +
    facet_wrap("eftype_format") +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  #plot(EFcompare_ensemble_best_vol_scatter)
  
  # EFcompare_all_best_vol_scatter <- ggplot(eftab_best_all[var=='qtot',],
  #                                               aes(x=(efvol_ref), y=efvol_mod)) +
  #   geom_abline(alpha=1/2) +
  #   geom_point() + #aes(group=Country,color=Country, shape=Country)
  #   scale_x_log10(name= expression('Reference e-flow volume'~(m^3~s^-1))) +
  #   scale_y_log10(name= expression('Modeled e-flow volume'~(m^3~s^-1))) +
  #   facet_grid(c("run", "eftype_format"), 
  #              labeller = "label_both",
  #              scales='free')
  # plot(EFcompare_all_best_vol_scatter)
  
  #By country--------------
  EFcompare_ensemble_best_vol_country_scatter <- ggplot(
    eftab_best_ensemble[var=='qtot',]) +
    # geom_ribbon(data=magnitude_lines_10, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='grey') + #'#6CDAE7'
    geom_abline(alpha=1/2) +
    geom_point(aes(x=(efvol_ref), y=efvol_mod, color=eftype_ref), 
               alpha=1/6) + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10(name= expression('Mean annual e-flows estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_log10(name= expression('Mean annual e-flows estimated by global models'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_color_discrete(name='Type of e-flow assessment') +
    geom_smooth(aes(x=efvol_ref, y=efvol_mod, color=eftype_ref),
                size=1, se=FALSE, 
                #color='black', 
                method='gam',
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    # geom_smooth(aes(x=efvol_ref, y=efvol_mod),
    #             span=1, se=F, color='black') +
    facet_grid(c("Country_labelo10", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    coord_fixed() +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  #plot(EFcompare_ensemble_best_vol_country_scatter)  
  
  #By eftype--------------
  EFcompare_ensemble_best_vol_eftype_scatter <- ggplot(
    eftab_best_ensemble[var=='qtot' &
                          !(eftype_ref %in% c('None specified', 'Other',
                                              'Hydraulic rating')),]) +
    # geom_ribbon(data=magnitude_lines_10, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='grey') + #'#6CDAE7'
    geom_abline(alpha=1/2) +
    geom_point(aes(x=(efvol_ref), y=efvol_mod, color=Country_labelo10),
               alpha=1/6) + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10(name= expression('Mean annual e-flows estimated by local e-flow assessment'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_y_log10(name= expression('Mean annual e-flows estimated by global models'~(m^3~s^-1)),
                  breaks=c(0.1, 10, 1000, 100000),
                  labels = trans_format("log10", math_format(10^.x)),
                  expand=c(0,0)) +
    scale_color_discrete(name='Country') +
    geom_smooth(aes(x=efvol_ref, y=efvol_mod, color=Country_labelo10),
                method='gam', se=F,
                formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)) +
    facet_grid(c("eftype_ref", "eftype_format"),
               scales = 'free_x',
               labeller = label_wrap_gen(width=15))+
    theme_classic() + 
    theme(panel.grid.major = element_line())
  #plot(EFcompare_ensemble_best_vol_eftype_scatter)  
  
  
  #All stats plots--------------
  
  statsplot_all_vol_smape <- ggplot(
    melt(efvol_stats_all, 
         id.vars=c('gcm', 'ghm', 'eftype_format'),
         measure.vars = c('smape', 'pbias', 'rsq')
    ),
    aes(x=eftype_format, y=value
        )) + 
    geom_point(aes(color=ghm, shape=gcm)) +
    geom_point(
      data = melt(efvol_stats_ensemble_all, 
           id.vars=c('eftype_format'),
           measure.vars = c('smape', 'pbias', 'rsq')
      ),
      color='black',
      size=3,
      alpha=1/2) +
    facet_wrap(~variable, scales='free') +
    theme_classic() +     
    theme(panel.grid.major = element_line())
  #statsplot_all_vol_smape
  
  
  efvol_stats_ensemble_country_bestworst <- merge(
    efvol_stats_ensemble_country_best,
    efvol_stats_ensemble_country_worst[, c('Country', 'eftype_format', 
                                           'smape', 'pbias')],
    by = c('Country', 'eftype_format'),
    suffixes = c("_best", "_worst"))
  
  statsplot_country_pbiassmape <- ggplot(
    efvol_stats_ensemble_country_bestworst, 
    aes(color=eftype_format)) + 
    geom_abline(slope=0, alpha=1/2, linetype='dashed') +
    # geom_point(aes(x=smape_worst, y=pbias_worst), 
    #            size=3, alpha=1/2) +
    geom_point(aes(x=smape_best, y=pbias_best), 
               size=3, alpha=1/3) +
    # geom_segment(aes(x=smape_best, xend=smape_worst,
    #                  y=pbias_best, yend=pbias_worst)) + 
    ggrepel::geom_text_repel(aes(x=smape_worst, y=pbias_worst, 
                                 label=str_wrap(label, width=10)),
                             #max.overlaps = Inf,
                             min.segment.length = Inf) +
    scale_x_continuous(name = 'Symmetric mean absolute percentage error (%)',
                       limits=c(0,125), breaks=seq(0,125,25), expand=c(0,0)) +
    scale_y_continuous(name = 'Bias (%)',
                       limits=c(-100,1000), breaks=seq(-100,1000, 100), 
                       expand=c(0,0)) +
    scale_color_discrete(name=str_wrap('Global e-flow estimation method',
                                       20)) +
    facet_wrap(~Country) +
    theme_classic() +      
    theme(panel.grid.major = element_line())
  #statsplot_country_pbias_hydrosub
  
  
  # statsplot_eftype_pbiassmape <- ggplot(
  #   efvol_stats_ensemble_eftype, 
  #   aes(x=smape, y=pbias, color=eftype_format)) + 
  #   geom_abline(slope=0, alpha=1/2) +
  #   geom_point(size=3, alpha=1/2) +
  #   ggrepel::geom_text_repel(aes(label=label)) +
  #   scale_x_continuous(limits=c(0,150), breaks=seq(0,150,25), expand=c(0,0)) +
  #   scale_y_continuous(limits=c(-100,1000), breaks=seq(-100,1000, 100), 
  #                      expand=c(0,0)) +
  #   facet_wrap(~eftype_ref) +
  #   theme_classic() +     
  #   theme(panel.grid.major = element_line())
  # #statsplot_eftype_pbiassmape
  
  statsplot_eftypecountry_pbiassmape <- ggplot(
    efvol_stats_ensemble_eftypecountry[n_total > 3,],
    aes(x=smape, y=pbias)) +
    geom_abline(slope=0, alpha=1/2) +
    #geom_line(aes(group=Country), alpha=1/2) +
    geom_point(aes(color=eftype_format, shape=Country), size=3) +
    # ggrepel::geom_text_repel(aes(label=str_wrap(label,15),
    #                              color=eftype_format,
    #                              max.overlaps = Inf)) +
    scale_shape_manual(values=c(8, 15,16,17,7, 18,9)) +
    scale_x_continuous(limits=c(50,150), breaks=seq(50,150,25), expand=c(0,0)) +
    scale_y_continuous(limits=c(-100,1000), breaks=seq(-100,1000, 100),
                       expand=c(0,0)) +
    facet_wrap(~eftype_ref) +
    theme_classic() + 
    theme(panel.grid.major = element_line())
  #statsplot_eftypecountry_pbiassmape
  
  
  ####################### Analysis by percentage ######################################################
  #--------------------- Stats -----------------------------------------
  #----- Only for ensemble estimates
  efper_stats_all <- eftab_best_all[
    !is.na(efper_ref) & 
      var=='qtot'& 
      gcm != 'cru-era' &
      #Ncountry > 10 &
      efper_ref >0
    ,
    getqstats(
      dt = .SD,
      actual = 'efper_ref',
      predicted =  'efper_mod',
      rstudthresh= 3,
      log = F
    ), 
    by=c('var', 'gcm', 'ghm', 
         'eftype_format')
  ] 
  
  efper_stats_ensemble_all <- eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot'& 
      #Ncountry > 10 &
      efper_ref >0
    ,
    getqstats(
      dt = .SD,
      actual = 'efper_ref',
      predicted =  'efper_mod',
      rstudthresh= 3,
      log = F
    ), 
    by=c('eftype_format')
  ] 
  

  efper_stats_ensemble_country <- eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot'& 
      Country_labelo10 != 'Brazil'
    ,
    getqstats(
      dt = .SD,
      actual = 'efper_ref',
      predicted =  'efper_mod',
      rstudthresh= 3,
      log = F
    ), 
    by=c('Country_labelo10', 'eftype_format')
  ] 
  
  efper_stats_ensemble_country[efper_stats_ensemble_country[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], by=Country_labelo10]$V1,
    label := eftype_format]
  
  efper_stats_ensemble_eftype <- eftab_best_ensemble[!is.na(efper_ref) & 
                                                       var=='qtot'& 
                                                       Ncountry > 10 &
                                                       eftype_ref != 'Other'
                                                     ,
                                                     getqstats(
                                                       dt = .SD,
                                                       actual = 'efper_ref',
                                                       predicted =  'efper_mod',
                                                       rstudthresh= 3,
                                                       log = F
                                                     ), 
                                                     by=c('eftype_ref',
                                                          'eftype_format')
  ] 
  
  efper_stats_ensemble_eftype[efper_stats_ensemble_eftype[
    , .I[(smape == min(smape)) | (abs(pbias) == min(abs(pbias)))], by=eftype_ref]$V1,
    label := eftype_format]
  
  #----------- Compute avg and sd of e-flow% by method -------
  efper_mod_summarystats_byeftype <- eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    list(efper_mod_mean = mean(efper_mod),
         efper_mod_sd = sd(efper_mod),
         efper_mod_min = min(efper_mod),
         efper_mod_max = max(efper_mod),
         .N),
    by=eftype_format
  ]
  
  eftab_best_ensemble[!is.na(efper_ref) & 
                        var=='qtot' &
                        !duplicated(EFUID),
                      list(efper_ref_mean = mean(efper_mod),
                           efper_ref_sd = sd(efper_mod),
                           efper_ref_min = min(efper_mod),
                           efper_ref_max = max(efper_mod),
                           N = .N),
                      by=eftype_ref
  ]
  
  efper_ref_summarystats_bycountry<- eftab_best_ensemble[!is.na(efper_ref) & 
                        var=='qtot'&
                        !duplicated(EFUID),
                      list(efper_ref_mean = mean(efper_ref),
                           efper_ref_sd = sd(efper_ref),
                           efper_ref_min = min(efper_ref),
                           efper_ref_max = max(efper_ref),
                           N = .N),
                      by=Country_labelo10
  ]
  
  eftab_best_ensemble[!is.na(efper_ref) & 
                        var=='qtot'&
                        !duplicated(EFUID),
                      .SD[efper_ref > 75, .N]/.N,
                      by=eftype_ref
  ]
  
  #----------- Compute % of over- and underestimation for each method and country -------
  efper_over_and_underestimation_ensemble_all <- eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    list(
      efper_over = .SD[efper_diff >= 10, .N]/.N,
      efper_under = .SD[efper_diff <= -10, .N]/.N
    ),
    by=c('eftype_format')
  ]
  
  efper_over_and_underestimation_ensemble_country <- eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    list(
      efper_over = .SD[efper_diff >= 10, .N]/.N,
      efper_under = .SD[efper_diff <= -10, .N]/.N
    ),
    by=c('Country_labelo10', 'eftype_format')
  ]
  
  eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    list(
      efper_over = .SD[efper_diff >= 10, .N]/.N,
      efper_under = .SD[efper_diff <= -10, .N]/.N
    ),
    by=c('eftype_format')
  ]
  
  eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    .SD[efper_diff > -10, .N]/.N]
  
  
  eftab_best_ensemble[
    !is.na(efper_ref) & 
      var=='qtot',
    .SD[efper_diff > -10, .N]/.N,
    by=c('eftype_format')
  ]
  
  #----------- Analyze relationship between e-flow % diff and discharge diff ----------------------------
  error_correlation_mar_efper_eftype_country_plot <- ggplot(eftab_best_ensemble,
         aes(x=mar_aperdiff, 
             y=efper_adiff)) + 
    scale_x_continuous('Absolute percentage difference between global 
                       and local estimates of mean annual flow') +
    scale_y_continuous('Absolute difference between global 
                       and local e-flow estimates as a percentage of MAF') +
    geom_point(aes(color=Country_labelo10)) + 
    geom_smooth(method='lm', color='black', se=F) +
    facet_grid(c("Country_labelo10", "eftype_format"),
               scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    #coord_fixed() +
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  
  error_correlation_mar_efper_eftype_country <- eftab_best_ensemble[
    !is.na(efper_diff) &
      !is.na(mar_aperdiff) &
      Ncountry > 10
    ,
    getqstats(
      dt = .SD,
      actual = 'efper_adiff',
      predicted =  'mar_aperdiff',
      rstudthresh= 3,
      log = T
    ), 
    by=c('eftype_format',
         'Country')
  ] 
  
  error_correlation_mar_efper_all <- eftab_best_ensemble[!is.na(efper_diff) &
                                                          !is.na(mar_aperdiff) &
                                                          Ncountry > 10
                                                        ,
                                                        getqstats(
                                                          dt = .SD,
                                                          actual = 'efper_adiff',
                                                          predicted =  'mar_aperdiff',
                                                          rstudthresh= 3,
                                                          log = T
                                                        )
  ] 
  
  
  #--------------------- Plots -----------------------------------------------------------------
  #All results----------------
  EFcompare_ensemble_best_per_scatter <- ggplot(
    eftab_best_ensemble[var=='qtot',]) +
    geom_abline(alpha=1/2) +
    geom_point(aes(x=efper_ref, y=efper_mod), alpha=1/3) + #aes(group=Country,color=Country, shape=Country)
    scale_x_continuous(name='E-flow as a % of mean annual flow estimated by local assessments',
                       breaks=c(0,20,40,60,80,100), limits=c(0,100)) +
    scale_y_continuous(name='E-flow as a % of mean annual flow estimated by global models',
                       breaks=c(0,20,40,60,80,100), limits=c(0,100)) +
    #scale_color_discrete(name='Country') +
    #geom_smooth(span=1) +
    #coord_cartesian(expand=c(0,0)) +
    facet_wrap("eftype_format") +
    theme_classic() +     
    theme(panel.grid.major = element_line())
  #plot(EFcompare_ensemble_best_per_scatter)
  
  # EFcompare_all_best_per_scatter <- ggplot(eftab_best_all[var=='qtot',],
  #                                               aes(x=(efper_ref), y=efper_mod)) +
  #   geom_abline(alpha=1/2) +
  #   geom_point() + #aes(group=Country,color=Country, shape=Country)
  #   scale_x_log10(name= expression('Reference e-flow volume'~(m^3~s^-1))) +
  #   scale_y_log10(name= expression('Modeled e-flow volume'~(m^3~s^-1))) +
  #   facet_grid(c("run", "eftype_format"), 
  #              labeller = "label_both",
  #              scales='free')
  # plot(EFcompare_all_best_per_scatter)
  
  #By country--------------
  
  EFcompare_ensemble_best_per_country_scatter <- ggplot(
    eftab_best_ensemble[var=='qtot',]) +
    # geom_ribbon(data=magnitude_lines, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2) +
    geom_point(aes(x=(efper_ref), y=efper_mod, color=Country_labelo10),
               alpha=1/3) + #aes(group=Country,color=Country, shape=Country)
    scale_x_continuous(name= expression(
      'E-flow as a % of mean annual flow estimated by local assessments',
      breaks=c(0,30,60,90)
    )) +
    scale_y_continuous(name= expression(
    'E-flow as a % of mean annual flow estimated by global models',
    breaks=c(0,30,60,90)
    )) +
    geom_smooth(aes(x=efper_ref, y=efper_mod),
                size=1, se=FALSE, 
                color='black', 
                method='lm'
                #method='gam',
                #formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)
                ) + 
    facet_grid(c("Country_labelo10", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    coord_fixed() +
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  #plot(EFcompare_ensemble_best_per_country_scatter)  
  
  
  #By eftype--------------
  EFcompare_ensemble_best_per_eftype_scatter <- ggplot(
    eftab_best_ensemble[var=='qtot' &
                          !(eftype_ref %in% c('None specified', 'Other',
                                              'Hydraulic rating')),]) +
    geom_abline(alpha=1/2) +
    geom_point(aes(x=efper_ref, y=efper_mod, 
                   color=Country_labelo10),
               alpha=1/6) + #aes(group=Country,color=Country, shape=Country)
    scale_x_continuous(name= expression(
      'E-flow as a % of mean annual flow estimated by local assessments',
      breaks=c(0,30,60,90)
    )) +
    scale_y_continuous(name= expression(
      'E-flow as a % of mean annual flow estimated by global models',
      breaks=c(0,30,60,90)
    )) +
    # geom_smooth(aes(x=efper_ref, y=efper_mod, color=Country_labelo10),
    #             se=F, size=1.5,
    #             method='lm') +
    geom_smooth(aes(x=efper_ref, y=efper_mod),
                size=1, se=FALSE, 
                color='black', 
                method='lm'
                #method='gam',
                #formula = y ~ s(x, bs = "cs", fx = TRUE, k = 3)
    ) + 
    facet_grid(c("eftype_ref", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    coord_fixed() +
    theme_classic() +
    theme(panel.grid.major = element_line())
          #legend.position = 'none')
  #plot(EFcompare_ensemble_best_per_eftype_scatter)  
  
  #Difference in percentage by MAR--------------
  EFMARdiff_ensemble_best_per_country <- ggplot(
    eftab_best_ensemble[var=='qtot',],
    aes(x=mar_ref, y=efper_diff)) +
    # geom_ribbon(data=magnitude_lines, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2, slope=0) +
    geom_point(aes(color=Country_labelo10)) + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10(name=expression('Mean annual flow estimated by local assessment'~m^3~s-1),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_y_continuous(name= 'Difference between global and local e-flow estimates (%)') +
    geom_smooth(method='lm', color='black', se=F) +
    facet_grid(c("Country_labelo10", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  #plot(EFMARdiff_ensemble_best_per_country )  
  
  #Difference in percentage by DA--------------
  EFDAdiff_ensemble_best_per_country <- ggplot(
    eftab_best_ensemble[var=='qtot',],
    aes(x=up_area_skm_15s, y=efper_diff)) +
    # geom_ribbon(data=magnitude_lines, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2, slope=0) +
    geom_point(aes(color=eftype_ref)) + #aes(group=Country,color=Country, shape=Country)
    scale_x_log10(name=expression('Upstream drainage area'~km^2),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_y_continuous(name= 'Difference between global and local e-flow estimates (%)') +
    geom_smooth(method='lm', color='black', se=F) +
    facet_grid(c("Country_labelo10", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  
  
  #Difference in percentage by aridity between countries --------------
  EFGAIdiff_ensemble_best_per_country <- ggplot(
    eftab_best_ensemble[var=='qtot',],
    aes(x=ari_ix_uav, y=efper_diff)) +
    # geom_ribbon(data=magnitude_lines, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2, slope=0) +
    geom_point(aes(color=Country_labelo10)) + #aes(group=Country,color=Country, shape=Country)
    scale_x_continuous(name = 'Global aridity index') + 
    scale_y_continuous(name= 'Difference between global and local e-flow estimates (%)') +
    geom_smooth(method='lm', color='black', se=F) +
    facet_grid(c("Country_labelo10", "eftype_format"),
               #scales = 'free_y',
               labeller = label_wrap_gen(width=15))+
    theme_classic() +
    theme(panel.grid.major = element_line(),
          legend.position = 'none')
  #plot(EFGAIdiff_ensemble_best_per_country)  
  
  #Difference in percentage by aridity between eftypes --------------
  EFGAIdiff_ensemble_best_per_eftype <- ggplot(
    eftab_best_ensemble[var=='qtot'  &
                          !(eftype_ref %in% c('None specified', 'Other',
                                              'Hydraulic rating')),],
    aes(x=ari_ix_uav, y=efper_diff)) +
    # geom_ribbon(data=magnitude_lines, aes(x=x, ymin=y_low, ymax=y_high), 
    #             alpha=1/3, fill='#6CDAE7') +
    geom_abline(alpha=1/2, slope=0) +
    geom_point(aes(color=Country_labelo10)) + #aes(group=Country,color=Country, shape=Country)
    scale_x_continuous() + 
    scale_y_continuous() +
    geom_smooth(span=1) +
    facet_grid(c("eftype_ref", "eftype_format"), 
               labeller = "label_both", scales = 'free') +
    theme_classic() +      theme(panel.grid.major = element_line())
  #plot(EFGAIdiff_ensemble_best_per_eftype)  
  
  
  #All stats plots--------------
  statsplot_all_per_smape <- ggplot(
    efper_stats_all, 
    aes(x=eftype_format, y=mae, 
        color=ghm, shape=gcm)) + 
    geom_point() +
    theme_classic() +      
    theme(panel.grid.major = element_line())
  #statsplot_all_per_smape
  
  
  # statsplot_country_pbiassmape <- ggplot(
  #   efper_stats_ensemble_country,
  #   aes(x=smape, y=pbias, color=eftype_format)) +
  #   geom_abline(slope=0, alpha=1/2) +
  #   geom_point(size=3, alpha=1/2) +
  #   ggrepel::geom_text_repel(aes(label=label)) +
  #   scale_x_continuous(limits=c(0,150), breaks=seq(0,150,25), expand=c(0,0)) +
  #   scale_y_continuous(limits=c(-100,1000), breaks=seq(-100,1000, 100),
  #                      expand=c(0,0)) +
  #   facet_wrap(~Country) +
  #       theme_classic() +      theme(panel.grid.major = element_line())
  # statsplot_country_pbias_hydrosub
  # 
  # 
  # statsplot_eftype_pbiassmape <- ggplot(
  #   efper_stats_ensemble_eftype,
  #   aes(x=smape, y=pbias, color=eftype_format)) +
  #   geom_abline(slope=0, alpha=1/2) +
  #   geom_point(size=3, alpha=1/2) +
  #   ggrepel::geom_text_repel(aes(label=label)) +
  #   scale_x_continuous(limits=c(0,150), breaks=seq(0,150,25), expand=c(0,0)) +
  #   scale_y_continuous(limits=c(-100,1000), breaks=seq(-100,1000, 100),
  #                      expand=c(0,0)) +
  #   facet_wrap(~eftype_ref) +
  #       theme_classic() +      theme(panel.grid.major = element_line())
  # statsplot_eftype_pbiassmape
  
  #------------- Write out statistics ------------------------------------------
  if (!dir.exists(outdir)) {
    dir.create(outdir)  
  }
  
  fwrite(efvol_stats_all, file.path(
    outdir, 
    paste0('efvol_stats_all',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_country, file.path(
    outdir, 
    paste0('efvol_stats_country',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_eftype, file.path(
    outdir, 
    paste0('efvol_stats_eftype',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_ensemble_all, file.path(
    outdir, 
    paste0('efvol_stats_ensemble_all',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_ensemble_country_best, file.path(
    outdir, 
    paste0('efvol_stats_ensemble_country_best',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_ensemble_country_worst, file.path(
    outdir, 
    paste0('efvol_stats_ensemble_country_worst',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_ensemble_eftype, file.path(
    outdir, 
    paste0('efvol_stats_ensemble_eftype',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_ensemble_eftypecountry, file.path(
    outdir, 
    paste0('efvol_stats_ensemble_eftypecountry',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efvol_stats_comparison_ensemble_comparison_marefvol, file.path(
    outdir, 
    paste0('efvol_stats_comparison_ensemble_comparison_marefvol',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_stats_all, file.path(
    outdir, 
    paste0('efper_stats_all',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_stats_ensemble_all, file.path(
    outdir, 
    paste0('efper_stats_ensemble_all',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_stats_ensemble_country, file.path(
    outdir, 
    paste0('efper_stats_ensemble_country',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_stats_ensemble_eftype, file.path(
    outdir, 
    paste0('efper_stats_ensemble_eftype',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_mod_summarystats_byeftype, file.path(
    outdir, 
    paste0('efper_mod_summarystats_byeftype',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_ref_summarystats_bycountry, file.path(
    outdir, 
    paste0('efper_ref_summarystats_bycountry',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_over_and_underestimation_ensemble_country , file.path(
    outdir, 
    paste0('efper_over_and_underestimation_ensemble_country',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(efper_over_and_underestimation_ensemble_all, file.path(
    outdir, 
    paste0('efper_over_and_underestimation_ensemble_all',
           format(Sys.Date(), '%Y%m%d'), '.csv')))
  fwrite(error_correlation_mar_efper_eftype_country, file.path(
    outdir, 
    paste0('error_correlation_mar_efper_eftype_country',
           format(Sys.Date(), '%Y%m%d'), '.csv')))

  
  #------------- Return output -------------------------------------------------
  return(list(
    efcomparison_general_stats = list(
      median_range_all_vol = avg_range_all[, median(efvol_mod_range)],
      range_all_vol_o1000 = avg_range_all[, .SD[efvol_mod_range>1000,.N]/.N] #Proportion of sites for which e-flow volume estimates ranged an order of magnitude.
    ),
    plot_efper_ecpresent_ref= plot_efper_ecpresent_ref,
    EFcompare_ensemble_best_vol_scatter = EFcompare_ensemble_best_vol_scatter,
    EFcompare_ensemble_best_vol_country_scatter = EFcompare_ensemble_best_vol_country_scatter,
    EFcompare_ensemble_best_vol_eftype_scatter = EFcompare_ensemble_best_vol_eftype_scatter,
    statsplot_all_vol_smape = statsplot_all_vol_smape, 
    # statsplot_country_pbiassmape,
    # statsplot_eftype_pbiassmape,
    # statsplot_eftypecountry_pbiassmape,
    EFcompare_ensemble_best_per_scatter = EFcompare_ensemble_best_per_scatter,
    EFcompare_ensemble_best_per_country_scatter = EFcompare_ensemble_best_per_country_scatter,
    EFcompare_ensemble_best_per_eftype_scatter = EFcompare_ensemble_best_per_eftype_scatter, 
    error_correlation_mar_efper_eftype_country_plot = error_correlation_mar_efper_eftype_country_plot,
    EFMARdiff_ensemble_best_per_country = EFMARdiff_ensemble_best_per_country,
    EFDAdiff_ensemble_best_per_country = EFDAdiff_ensemble_best_per_country, 
    EFGAIdiff_ensemble_best_per_country = EFGAIdiff_ensemble_best_per_country,
    EFGAIdiff_ensemble_best_per_eftyp = EFGAIdiff_ensemble_best_per_eftype,
    statsplot_all_per_smape = statsplot_all_per_smape 
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
    theme_classic() +      theme(panel.grid.major = element_line()) +
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
