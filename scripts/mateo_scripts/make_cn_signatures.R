#!/usr/bin/env Rscript


list.of.packages <- c("pheatmap", "tidyverse", "NMF", "flexmix", "QDNAseq", "corrplot", "GGally", "YAPSA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))
library(NMF)
library(flexmix)
library(QDNAseq)

library(corrplot)
library(GGally)
library(YAPSA)

args = commandArgs(trailingOnly=TRUE)

source("functions.R")
source("main_functions.R")

args[1] <- "../../Datasets/TCGA/PanCancer/TCGA_ASCAT_RAW_PVL/filteredAscat.rds" # ascat
args[2] <- "../../Datasets/TCGA/TCGA_survival_data_clean.txt"

cancerType <- "BLCA"

TCGA_BLCA_CIN_measures <- readRDS("TCGA_BLCA_CIN_measures.rds")


clinical <- read.delim2(args[2], sep = '\t') #get(load(args[2]))
clinical <- clinical %>% filter(type == "BLCA")

######### Ascat data analysis ################
cat("Reading Ascat Table\n")
ascat <- readRDS(args[1]) #get(load(args[1]))
test_ascat <- checkAscatData(ascat)
if(test_ascat$success == 1){        # everything good
  ascat <- test_ascat$ascat
} else if(test_ascat$success == 0){ # colnames do not match
  ascat <- changeAscatColnames(ascat)
} else {                            # other error
  stop("Stopping...Ascat table does not have correct format or order...")
}

cat("Filtering Chr 23 and 24 from Ascat...\n")
ascat <- ascat %>% filter(Chr != 23, Chr != 24)
#ascat$ID <- substr(ascat$ID, 3,8)

cat("Keeping rows with Abb.Cell.Frac < 0.9 and > 0.2 from Ascat...\n")
ascat <- ascat %>%
  filter(Aberrant.Cell.Fraction < 0.95) %>%
  filter(Aberrant.Cell.Fraction > 0.2)


tmp_ascat <- ascat %>% filter(cancer_type == cancerType)


colnames(tmp_ascat) <- c("ID", "chromosome", "start", "end", "nProbes", "cn","nA","nB","Ploidy","Aberrant.Cell.Fraction", "type") # TCGA have Type column
tmp_ascat$segVal <- log2((tmp_ascat$nA + tmp_ascat$nB)/tmp_ascat$Ploidy + 1)

list_of_segments <- list()
cnt = 1

for(id in unique(tmp_ascat$ID)) {
  tmp <- tmp_ascat %>%
    filter(ID == id)
  tmp <- tmp %>% dplyr::select("chromosome", "start", "end", "segVal")
  
  list_of_segments[[cnt]] <- tmp
  cnt = cnt + 1
  #if(cnt > 25) break
}
#"entire_TCGA_cnsigs"
folder <- "FINALoutput_cnsignatures_minprior_default_sigs_7"
output_dir <- file.path(getwd(), folder)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

cores <- 1
names(list_of_segments) <- unique(tmp_ascat$ID)#[1:25]

if(!file.exists(paste0(folder,"/BLCA_cn_features.rds"))){
  CN_features <- extractCopynumberFeatures(list_of_segments)
  saveRDS(CN_features, paste0(folder,"/BLCA_cn_features.rds"))
} else {
  CN_features <- readRDS(paste0(folder,"/BLCA_cn_features.rds"))
}

if(!file.exists(paste0(folder,"/BLCA_allcomponents.rds"))){
  all_components <- fitMixtureModels(CN_features, niter = 1000, cores = 1, nrep = 20, min_prior = 0, max_comp = 5, min_comp = 2) # this is optional
  cat("Checking for all components...\n")
  for(cmp in all_components){
    if(!cmp@converged){
      cat("Algorithm did not converged\n")
      print(cmp)
    }
    if(length(unique(cmp@cluster)) < 2){
      cat("# of Clusters is less than 2\n")
      print(cmp)
    }
  }
  saveRDS(all_components, paste0(folder,"/BLCA_allcomponents.rds"))
} else {
  all_components <- readRDS(paste0(folder, "/BLCA_allcomponents.rds"))
}

if(!file.exists(paste0(folder, "/sample_by_component.rds"))){
  sample_by_component <- generateSampleByComponentMatrix(CN_features = CN_features, all_components = all_components, cores = cores)
  saveRDS(sample_by_component, paste0(folder,"/sample_by_component.rds"))
} else {
  sample_by_component <- readRDS(paste0(folder,"/sample_by_component.rds"))
}

# if(!file.exists( paste0(folder, "/chooseNumberSignature.rds"))){
#   num <- chooseNumberSignatures(sample_by_component,min_sig = 3, max_sig = 9, cores = cores)
#   saveRDS(num, paste0(folder,"/chooseNumberSignature.rds"))
# } else {
#   num <- readRDS(paste0(folder,"/chooseNumberSignature.rds"))
# }

if(!file.exists(paste0(folder, "/BLCA_component_by_signature_sigs.rds"))){
  component_by_signature <- generateSignatures(sample_by_component, nsig = 7, cores = cores)
  pheatmap(component_by_signature@consensus)
  
  saveRDS(component_by_signature, paste0(folder, "/BLCA_component_by_signature_sigs.rds"))
} else {
  component_by_signature <- readRDS(paste0(folder,"/BLCA_component_by_signature_sigs.rds"))
}

if(!file.exists( paste0(folder,"/for_heatmap.rds"))){
  for_heatmap <- quantifySignatures(sample_by_component = sample_by_component, component_by_signature@fit@W)
  saveRDS(for_heatmap, paste0(folder,"/for_heatmap.rds"))
} else {
  for_heatmap <- as.data.frame((readRDS(paste0(folder,"/for_heatmap.rds"))))
}

pheatmap((for_heatmap), 
         #annotation_col = df_columns,
         cluster_rows = F,
         cutree_cols = 3, 
         clustering_method = "ward.D2",
         scale= "row"
)

