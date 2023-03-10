source("R/packages.R")
library(rprojroot)
library(magrittr)


rootdir <- rprojroot::find_root(has_dir('src'))
datdir <- file.path(rootdir, 'data')
flowdat <- fread(file.path(datdir,'X0434010_QmnJ.csv'))



quantile(c(rep(0,116), 0.8411439, 0.1940759, 0.9927341, 0.6467704),
         c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 
           0.3, 0.4, 0.5 ,0.6 ,0.7,
           0.8 ,0.9 ,0.95 ,0.99 ,0.999, 0.9999))

get_fdc <- function(in_dt, na_value, datecol, discol, min_n) {
  dt <- copy(in_dt)
  dt_format <- dt[, year_month := substr(get(datecol), 1, 7)] %>%
    .[, list(mean_dis = mean(.SD[get(discol) != na_value, get(discol)]),
             n = sum(get(discol) != na_value)), 
      by=year_month]
  
  out_fdc <- data.table(quant = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 
                                    0.3, 0.4, 0.5 ,0.6 ,0.7,
                                    0.8 ,0.9 ,0.95 ,0.99 ,0.999, 0.9999))
  out_fdc[, `:=`(discharge_natural = quantile(dt_format[n > min_n, mean_dis],
                                              quant),
                 exceedance_prob = 1-quant)]
  return(list(fdc=out_fdc, ts_format = dt_format))
}

fdc_gefis <- get_fdc(in_dt = flowdat, na_value = -999, 
                     datecol = 'Date (TU)', discol = 'Valeur (en l/s)', min_n = 28)

fdc_in <- fdc_gefis
n_shift <- 1

shift_ts <- function(fdc_in, n_shift) {
  newcol <- paste0('discharge_', LETTERS[n_shift])
  newcol_bin <- paste0('bin', newcol)
  fdc_in$fdc[, (newcol_bin) := shift(discharge_natural, n=n_shift)] %>%
    .[, I := .I]
  
  #When the flow is between 99.99% and 0.01% exceedance probability
  ts_binned <- fdc_in$ts_format[mean_dis >= min(fdc_in$fdc$discharge_natural) &
                                  mean_dis <= max(fdc_in$fdc$discharge_natural), 
                                bin := findInterval(mean_dis, fdc_in$fdc[,discharge_natural])]

  ts_binned[!is.na(bin), 
            (newcol) := 
              stats::approx(x = log(c(fdc_in$fdc[bin, discharge_natural],
                                      fdc_in$fdc[bin+n_shift, discharge_natural])),
                            y = c(fdc_in$fdc[bin, get(newcol_bin)],
                                  fdc_in$fdc[bin+n_shift, get(newcol_bin)]), 
                            xout = log(mean_dis), method = "linear")$y
              
  ]
  
  
  


ggplot(ts_binned, aes(x=as.Date(paste0(year_month, '-01')),
                                y=mean_dis
)) + 
  geom_line() +
  geom_line(aes(y=discharge_C), color='red') +
  labs(x='Month', y=expression(paste('Monthly mean discharge  ',(m^3~s^-1)))) +
  scale_y_log10() +
  coord_cartesian(expand=c(0,0)) +
  theme_classic()

monthly_plot <- ggplot(fdc_gefis$ts_format, aes(x=as.Date(paste0(year_month, '-01')),
                                        y=mean_dis
)) + 
  geom_line() +
  labs(x='Month', y=expression(paste('Monthly mean discharge  ',(m^3~s^-1)))) +
  coord_cartesian(expand=c(0,0)) +
  theme_classic()

ggplot(fdc_gefis$fdc, aes(x=as.character(100*exceedance_prob), y=discharge_natural)) + 
  geom_point() +
  geom_line() +
  labs(x='Time flow exceeded (%)', 
       y=expression(paste('Discharge  ',(m^3~s^-1)))) +
  scale_y_log10() +
  theme_classic()
}

# 
# quantile(ts_binned$mean_dis, 0.9999)
# max(ts_binned$mean_dis)
