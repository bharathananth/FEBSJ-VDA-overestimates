####----------------------Library loading----------------------------------####
library(preprocessMA)
library(compareRhythms)
library(genefilter)
library(SummarizedExperiment)
library(tximport)
library(readr)
library(stringr)

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

####----------------------Lung inflammation data-------------------------####

# No brainarray yet for this dataset

# GEO <- "GSE59404"
# geo <- get_geo(GEO)
# raw_file_name <- get_raw_files(GEO)
# pkg <- install_brainarray("mogene10st", annot = "ensg")
# 
# cel_files <- as.character(colData(geo)[["supplementary_file"]])
# #--fix misreferenced file in supplementary table--#
# err_file_name <- cel_files[grepl("GSM1263252", cel_files)]
# cel_files[grepl("GSM1263252", cel_files)] <- sub("pli[\\w\\-\\.]+", "CEL.gz", sub("bA", "bB", err_file_name), perl=TRUE)
# eset <- read_celfiles(raw_file_name, cel_files,platform_design = pkg)
# data <- update_se(geo, eset)
# 
# saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")



####------------------------Make tx2gene dataframe for tximport----------------------------####
# txDb <- GenomicFeatures::makeTxDbFromGFF("~/Downloads/gencode.vM24.primary_assembly.annotation.gff3",
#                                          organism = "Mus musculus")
# k <- keys(txDb, keytype = "TXNAME")
# tx2gene <- select(txDb, k, "GENEID", "TXNAME")
# write.csv(tx2gene, file="data/tx2gene.gencode.vM24.csv", row.names = FALSE)

tx2gene <- read_csv("data/tx2gene.gencode.vM24.csv")

# annot_db <- biomaRt::useEnsembl("ensembl",
#                                 dataset = "mmusculus_gene_ensembl",
#                                 mirror = "useast")
# 
# attributes <- c('ensembl_gene_id', "entrezgene_id", "mgi_symbol",
#                 "chromosome_name", "strand",
#                 "start_position", "end_position")
# 
# row_data <- biomaRt::getBM(attributes = attributes,
#                            filters = 'ensembl_gene_id',
#                            values = str_match(tx2gene$GENEID, "^ENSMUSG\\d+")[,1],
#                            mart = annot_db)
# 
# row_data <- row_data[!BiocGenerics::duplicated(row_data$ensembl_gene_id), ]
# 
# rownames(row_data) <- row_data$ensembl_gene_id
# 
# write.csv(row_data, file="data/row_data.csv")

row_data <- read.csv("data/row_data.csv", row.names = 1)


####---------------------GR cistrome from Quagliarini et al. 2019-------------------------####
GEO <- "PRJNA428303"
metadata <- read_csv("data/PRJNA428303_SraRunTable.csv")

files <- paste0("data/", metadata$Run, "_quant/quant.sf")
names(files) <- metadata$Run

txi_quag <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                     countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE)

countsFromAbundance <- txi_quag$counts

txi_quag <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                     countsFromAbundance = "no", ignoreAfterBar = TRUE)

counts <- txi_quag$counts
lengths <- txi_quag$length
abundance <- txi_quag$abundance

countsFromAbundance <- countsFromAbundance[rowSums(counts)>0, ]
lengths <- lengths[rowSums(counts)>0, ]
abundance <- abundance[rowSums(counts)>0, ]
counts <- counts[rowSums(counts)>0, ]

rownames(counts) <- str_match(rownames(counts), "^ENSMUSG\\d+")[,1]
rownames(countsFromAbundance) <- str_match(rownames(countsFromAbundance), "^ENSMUSG\\d+")[,1]
rownames(lengths) <- str_match(rownames(lengths), "^ENSMUSG\\d+")[,1]
rownames(abundance) <- str_match(rownames(abundance), "^ENSMUSG\\d+")[,1]

data <- SummarizedExperiment(assays = list(counts = counts, countsFromAbundance = countsFromAbundance, 
                                           lengths = lengths, abundance = abundance),
                             colData = metadata,
                             rowData = row_data[rownames(counts), ])

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")


####---------------------Liver independence from Koronowski et al. 2019--------------------####
GEO <- "SRP153814"

metadata <- read_csv(paste0("data/", GEO, "_SraRunTable.csv"))

files <- paste0("data/", metadata$Run, "_quant/quant.sf")
names(files) <- metadata$Run

txi_koron <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                     countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE)

countsFromAbundance <- txi_koron$counts

txi_koron <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                      countsFromAbundance = "no", ignoreAfterBar = TRUE)

counts <- txi_koron$counts
lengths <- txi_koron$length
abundance <- txi_koron$abundance

countsFromAbundance <- countsFromAbundance[rowSums(counts)>0, ]
lengths <- lengths[rowSums(counts)>0, ]
abundance <- abundance[rowSums(counts)>0, ]
counts <- counts[rowSums(counts)>0, ]

rownames(counts) <- str_match(rownames(counts), "^ENSMUSG\\d+")[,1]
rownames(countsFromAbundance) <- str_match(rownames(countsFromAbundance), "^ENSMUSG\\d+")[,1]
rownames(lengths) <- str_match(rownames(lengths), "^ENSMUSG\\d+")[,1]
rownames(abundance) <- str_match(rownames(abundance), "^ENSMUSG\\d+")[,1]

data <- SummarizedExperiment(assays = list(counts = counts, countsFromAbundance = countsFromAbundance,
                                           abundance = abundance, lengths = lengths),
                             colData = metadata,
                             rowData = row_data[rownames(counts), ])

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")


####---------------------p66 redox from Pei et al. 2020----------------------####
GEO <- "PRJNA449625"

metadata <- read_csv(paste0("data/", GEO, "_SraRunTable.csv"))

files <- paste0("data/", metadata$Run, "_quant/quant.sf")
names(files) <- metadata$Run

txi_pei <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                      countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE)

countsFromAbundance <- txi_pei$counts

txi_pei <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                    countsFromAbundance = "no", ignoreAfterBar = TRUE)

counts <- txi_pei$counts
abundance <- txi_pei$abundance
lengths <- txi_pei$length

countsFromAbundance <- countsFromAbundance[rowSums(counts)>0, ]
lengths <- lengths[rowSums(counts)>0, ]
abundance <- abundance[rowSums(counts)>0, ]
counts <- counts[rowSums(counts)>0, ]

rownames(counts) <- str_match(rownames(counts), "^ENSMUSG\\d+")[,1]
rownames(countsFromAbundance) <- str_match(rownames(countsFromAbundance), "^ENSMUSG\\d+")[,1]
rownames(lengths) <- str_match(rownames(lengths), "^ENSMUSG\\d+")[,1]
rownames(abundance) <- str_match(rownames(abundance), "^ENSMUSG\\d+")[,1]

data <- SummarizedExperiment(assays = list(counts = counts, countsFromAbundance = countsFromAbundance, lengths = lengths, abundance = abundance),
                             colData = metadata,
                             rowData = row_data[rownames(counts), ])

saveRDS(data, file = paste0("data/", GEO, ".RDS"), compress = "gzip")


####----combining WT datasets (HFD, KD) from two different experiments------####

GEO1 <- "GSE52333"
geo1 <- get_geo(GEO1)
raw_file_name1 <- get_raw_files(GEO1)

GEO2 <- "GSE87425"
geo2 <- get_geo(GEO2)
raw_file_name2 <- get_raw_files(GEO2)

cel_files1 <- as.character(colData(geo1)[["supplementary_file"]])
cel_files2 <- as.character(colData(geo2)[["supplementary_file"]])

err_file_name <- cel_files1[grepl("GSM1263252", cel_files1)]
cel_files1[grepl("GSM1263252", cel_files1)] <- sub("pli[\\w\\-\\.]+", "CEL.gz", sub("bA", "bB", err_file_name), perl=TRUE)


cel_files1 <- sapply(cel_files1, basename)
cel_files2 <- sapply(cel_files2, basename)

cel_files1_full <- paste(dirname(raw_file_name1), cel_files1, sep="/")
cel_files2_full <- paste(dirname(raw_file_name2), cel_files2, sep="/")

raw_data <- oligo::read.celfiles(c(cel_files1_full, cel_files2_full),
                                 pkgname = "pd.mogene10st.mm.ensg")

eset <- oligo::rma(raw_data)
col_data1 <- colData(geo1) %>% subset(select = c(diet.ch1, tissue.ch1, harvest.timepoint.ch1))
col_data1$harvest.timepoint.ch1 <- as.numeric(sub("ZT", "", col_data1$harvest.timepoint.ch1))
col_data1$zt.ch1 <- col_data1$harvest.timepoint.ch1
col_data1$harvest.timepoint.ch1 <- NULL
col_data2 <- colData(geo2) %>% subset(select = c(diet.ch1, tissue.ch1, zt.ch1))
geo <- SummarizedExperiment(colData = rbind(col_data1, col_data2))

data <- update_se(geo, eset)
data <- subset(data, select = colData(data)$tissue.ch1 %in% c("liver", "LIVER") & colData(data)$diet.ch1 %in% c("normal chow", "Norma Chow"))
saveRDS(data, file = paste0("data/", "WT_compare", ".RDS"), compress = "gzip")
