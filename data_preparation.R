####----------------------Library loading----------------------------------####
library(preprocessMA)
library(compareRhythms)
library(genefilter)
library(SummarizedExperiment)



####----------------------High resolution liver data------------------------####

GEO <- "GSE11923"
geo <- get_geo(GEO)
raw_file_name <- get_raw_files(GEO)
pkg <- install_brainarray("mouse4302", annot = "ensg")
cel_files <- as.character(SummarizedExperiment::colData(geo)[["supplementary_file"]])
eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
data <- update_se(geo, eset)

expt_design <- do.call(rbind, strsplit(as.character(data$title), 
                                                  " ", fixed = TRUE))

data$time <- as.numeric(expt_design[, 3])

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")

####----------------------High fat diet liver data-------------------------####

GEO <- "GSE52333"
geo <- get_geo(GEO)
raw_file_name <- get_raw_files(GEO)
pkg <- install_brainarray("mogene10st", annot = "ensg")



