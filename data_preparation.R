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

cel_files <- as.character(colData(geo)[["supplementary_file"]])
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

cel_files <- as.character(colData(geo)[["supplementary_file"]])
#--fix misreferenced file in supplementary table--#
err_file_name <- cel_files[grepl("GSM1263252", cel_files)]
cel_files[grepl("GSM1263252", cel_files)] <- sub("pli[\\w\\-\\.]+", "CEL.gz", sub("bA", "bB", err_file_name), perl=TRUE)
eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
data <- update_se(geo, eset)

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")


####----------------------Aging in liver data------------------------------####

GEO <- "GSE93903"
geo <- get_geo(GEO)
raw_file_name <- get_raw_files(GEO)
pkg <- install_brainarray("htmg430pm", annot = "ensg")

cel_files <- as.character(colData(geo)[["supplementary_file"]])
eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
data <- update_se(geo, eset)

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")

####----------------------Lung Cancer effects on liver data------------------------------####

GEO <- "GSE73222"
geo <- get_geo(GEO)
raw_file_name <- get_raw_files(GEO)
pkg <- install_brainarray("mogene20st", annot = "ensg")

cel_files <- as.character(colData(geo)[["supplementary_file"]])
eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
data <- update_se(geo, eset)

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")

####----------------------Ketogenic diet on Liver and Gut data------------------------------####

GEO <- "GSE87425"
geo <- get_geo(GEO)
raw_file_name <- get_raw_files(GEO)
pkg <- install_brainarray("mogene10st", annot = "ensg")

cel_files <- as.character(colData(geo)[["supplementary_file"]])
eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
data <- update_se(geo, eset)

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")



