source('R/packages.R')
rootdir = rprojroot::find_root(has_dir('src'))
outdir <- file.path(rootdir, 'data', 'GEFIS_test_data', 'Data by Country', 'Europe', 'France', 'Rhone')


studies_page <- "https://www.rhone-mediterranee.eaufrance.fr/gestion-de-leaugestion-quantitative-de-la-ressource-en-eauetudes-volumes-prelevables/etudes-0"
studies_html <- rvest::read_html(studies_page)
rapports_elems <- rvest::html_elements(studies_html, "ul") %>%
  html_elements('a')

rapports_elems_id <- rapports_elems %>%
  html_attr('data-nodeid') %>%
  is.na %>%
  !.

studies_page_host <- httr::parse_url(studies_page)
studies_page_host$path <- NULL

region_pages <- rapports_elems[rapports_elems_id] %>%
  html_attr('href') %>%
  xml2::url_absolute(., base=httr::build_url(studies_page_host))

process_rhone_eflow_page <- function(url, outdir, verbose = T) {
  if (verbose) print(url)
  study_page <- rvest::read_html(url)
  title <- rvest::html_element(study_page, '.title') %>%
    html_text() %>%
    gsub("\r?\n|\r", '', .) %>% #Remove line break at the end of title
    gsub("[/]", "-", .)  
    
  url_host <- httr::parse_url(url)
  url_host$path <- NULL
  
  pdfs_url <- rvest::html_elements(study_page, 'a') %>%
    rvest::html_attr('href') %>%
    grep("[.]pdf$", ., value=T) %>%
    xml2::url_absolute(., base=httr::build_url(url_host))
  
  
  outdir_full <- file.path(outdir, title)
  if (!(dir.exists(outdir_full))) {
    dir.create(outdir_full)
  }
  
  lapply(pdfs_url, function(in_url) {
    out_pdf <- file.path(outdir_full, 
                         tail(strsplit(in_url, split='/')[[1]], n=1))
    if (!file.exists(out_pdf)) {
      download.file(url = in_url, destfile = out_pdf, mode = 'wb')
    } else {
      print(paste(out_pdf, 'was already downloaded.'))
    }
  })
  
  return(
    list(title = title,
         pdfs = pdfs_url
    )
  )
}

lapply(region_pages, process_rhone_eflow_page, outdir=outdir)

############### Download shapefile of Masses d'eau ##############################
sdage_masses_url = 'https://www.rhone-mediterranee.eaufrance.fr/sites/sierm/files/content/2019-07/RefSdage2016_mdoriv_FRD.zip'
download.file(sdage_masses_url, file.path(outdir, basename(sdage_masses_url)))
unzip(file.path(outdir, basename(sdage_masses_url)), exdir = outdir)

############### Get coordinates - to be continued ##############################
etat_cours_deau_db <- "https://www.rhone-mediterranee.eaufrance.fr/surveillance-des-eaux/qualite-des-cours-deau/donnees-detat-des-cours-deau-superficiels?search_api_fulltext=&field_sq_region=All&field_sq_departement=All&field_sq_commune=All&field_sq_watershed=All&field_sq_subwatershed=All&field_sq_stream=All&items_per_page=90&view_mode=medium&page=1"
rvest::read_html(etat_cours_deau_db)
