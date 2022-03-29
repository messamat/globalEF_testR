
source('R/packages.R')
source('R/functions.R')

figdir <- file.path(resdir, 'figures')
if (!dir.exists(figdir)) {
  dir.create(figdir)
}

tar_load(eftab)
qplot(eftab$station_river_distance)

  
#Number of EF assessments
eftab[, .N]

#Number of unique sites:
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), .N]

#Number of countries: 
eftab[, length(unique(Country))]
  
#Characteris of unique sites
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), median(UPLAND_SKM)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), median(dis_m3_pyr)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(ppd_pk_uav)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(crp_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(urb_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(for_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), sd(for_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), sd(for_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(pac_pc_use)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]) & pac_pc_use > 50, .N]/eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]),.N]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]) & pac_pc_use < 10, .N]/eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]),.N]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), mean(dor_pc_pva)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]), median(dor_pc_pva)]
eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]) & dor_pc_pva > 10, .N]/eftab[!duplicated(eftab[, .(POINT_X, POINT_Y)]),.N]

#EF type
eftab[, .N, by=eftype_ref]
eftab[, length(unique(efname_ref))]

eftab[, length(unique(ecpresent_ref_formatted))]

eftab[, min(dis_m3_pyr)]
eftab[, min(mar_ref, na.rm=T)/31.56]
eftab[, max(dis_m3_pyr)]
eftab[, max(mar_ref, na.rm=T)/31.56]

#
eftab[ecpresent_ref_formatted == 'A/B', mean(dor_pc_pva)]
eftab[ecpresent_ref_formatted == 'D', mean(dor_pc_pva)]

tar_load(envhist)
pdf(file.path(figdir, paste0('envhist', format(Sys.Date(), '%Y%m%d'), '.pdf')),
    width = 10,
    height = 8)
grid.draw(envhist)
dev.off()


tar_load(countrytab)
png(file.path(figdir, paste0('countryhist', format(Sys.Date(), '%Y%m%d'), '.png')),
    width = 5, height = 5, unit='in', res=600)
grid.draw(countrytab)
dev.off()


tar_load(efmap)
png(file.path(figdir, paste0('efmap', format(Sys.Date(), '%Y%m%d'), '.png')),
    width = 6, height = 5, unit='in', res=600)
grid.draw(efmap)
dev.off()


tar_load(hydrology_comparison)
png(file.path(figdir, paste0('hydrology_comparison', format(Sys.Date(), '%Y%m%d'), '.png')),
    width = 6, height = 6, unit='in', res=600)
grid.draw(hydrology_comparison$plot)
dev.off()

gt(hydrology_comparison$table_all) %>%
  gtsave(file.path(figdir, paste0('hydrology_comparison_tableall',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))


gt(hydrology_comparison$table_clz) %>%
  gtsave(file.path(figdir, paste0('hydrology_comparison_tableclz',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))

tar_load(EMC_comparison)
png(file.path(figdir, paste0('EMC_comparison_ws', format(Sys.Date(), '%Y%m%d'), '.png')),
    width = 6, height = 8, unit='in', res=600)
grid.draw(EMC_comparison$ECws_boxplot)
dev.off()

png(file.path(figdir, paste0('EMC_comparison_p', format(Sys.Date(), '%Y%m%d'), '.png')),
    width = 6, height = 8, unit='in', res=600)
grid.draw(EMC_comparison$ECp_boxplot)
dev.off()

pdf(file.path(figdir, paste0('EMC_humanstressors', format(Sys.Date(), '%Y%m%d'), '.pdf')),
    width = 7,
    height = 6)
grid.draw(EMC_comparison$drivers_boxplot)
dev.off()

gt(EMC_comparison$confumat_max) %>%
  gtsave(file.path(figdir, paste0('EMC_comparison_tableall',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))


tar_load(EFestimate_comparison)
pdf(file.path(figdir, paste0('EFcompare_originalc', format(Sys.Date(), '%Y%m%d'), '.pdf')),
    width = 6,
    height = 6)
grid.draw(EFestimate_comparison$EFcompare_originalec_p)
dev.off()

pdf(file.path(figdir, paste0('EFcompare_bestc', format(Sys.Date(), '%Y%m%d'), '.pdf')),
    width = 7,
    height = 7)
grid.draw(EFestimate_comparison$EFcompare_bestec_p)
dev.off()


gt(EFestimate_comparison$table_eftype_best) %>%
  gtsave(file.path(figdir, paste0('EFestimate_comparisontable_eftype_best',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))
gt(EFestimate_comparison$table_eftype_worst) %>%
  gtsave(file.path(figdir, paste0('EFestimate_comparisontable_eftype_worst',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))
gt(EFestimate_comparison$table_emcref) %>%
  gtsave(file.path(figdir, paste0('EFestimate_comparisontable_emcref',
                                  format(Sys.Date(), '%Y%m%d'), '.html')))
