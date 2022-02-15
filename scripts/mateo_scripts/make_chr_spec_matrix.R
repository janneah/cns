setwd("~/GenomeDK_local/CancerEvolution/phd/Analysis/cnsignatures")
library(doParallel)
library(tidyverse)
library(gtools)
source("../cnsignatures/functions.R")


args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (cancer type)", call. = FALSE)
}

param_to_measure <- "segVal"
bin_size <- 1e6 # 0.5e7
args[1] <- "LUAD"
type <- args[1]

chrominfo_hg19 <- as.data.frame(chrominfo_hg19)
# Read the data and filter for downsampling
ascat <- readRDS(
  "../../Datasets/TCGA/PanCancer/TCGA_ASCAT_RAW_PVL/filteredAscatWithRaw.rds")
ascat$segval <- log2((ascat$nAraw + ascat$nBraw) / ascat$Ploidy)
ascat$segval[which(ascat$segVal == -Inf)] <- -max(ascat$segVal)
ascat$id <- as.character(ascat$ID)

output_pancan <- data.frame()
types <- args[1]
for (t in types) {
  cat(t, "\n")
  t_ascat <- ascat %>% filter(cancer_type == t)
  t_ascat <- t_ascat %>% filter(Chr != 23, Chr != 24)
  cat("Keeping rows with Abb.Cell.Frac < 0.9 and > 0.2 from Ascat...\n")
  t_ascat <- t_ascat %>%
    filter(Aberrant.Cell.Fraction < 0.99) %>%
    filter(Aberrant.Cell.Fraction > 0.2)
  seg_big <- t_ascat[(t_ascat[, 4] - t_ascat[, 3]) >= 3e6, ]
  seg_small <- t_ascat[(t_ascat[, 4] - t_ascat[, 3]) < 3e6, ]
  seg_shrink <- shrink.seg.ai.wrapper(seg_big)
  t_ascat <- rbind(seg_shrink, seg_small)
  mean(t_ascat$End - t_ascat$Start)
#ascat <- ascat[which(ascat$ID %in% TCGA_BLCA_CIN_measures$sample_id),] # nolint

  segVal_mat <- calculateChrEvents(t_ascat, chrominfo_hg19 = chrominfo_hg19,
                                 binSize = binSize, 
                                 paramToMeasure = "segVal")
  output_pancan <- rbind(output_pancan, segVal_mat)

}

id_num <- 1
tmp <- reshape::melt(output[id_num, ])
tmp$color <- ifelse(tmp$value == 0, "0", ifelse(tmp$value < 0, "Loss", "Gain"))
tmp$Chr <- as.numeric(unlist(lapply(tmp$variable, function(x) {
  split <- str_split(x, pattern = "_")
  substr(split[[1]][1], start = 4, 5)
})))
tmp$position <- 1: nrow(tmp)
tt <- tmp %>% 
  group_by(Chr) %>%
  summarise(n = n())
num <- 250
for (i in 2:22) {
  num <- num + tt$n[i]
  cat(num, " ")
}
p2 <- ggplot(tmp) +
  geom_point(aes(x = position, y = value, color = color)) +
  geom_col(aes(x = position, y = value, fill = color)) +
  geom_vline(xintercept = 250) + # 25
  geom_vline(xintercept = 494) + # 25
  geom_vline(xintercept = 693) + # 20
  geom_vline(xintercept = 885) + # 20
  geom_vline(xintercept = 1066) + # 19
  geom_vline(xintercept = 1238) + # 18
  geom_vline(xintercept = 1398) + # 19
  geom_vline(xintercept = 1542) + # 16
  geom_vline(xintercept = 1687) + # 15
  geom_vline(xintercept = 1823) + # 15
  geom_vline(xintercept = 1959) + # 14
  geom_vline(xintercept = 2093) + # 14
  geom_vline(xintercept = 2209) + # 14
  geom_vline(xintercept = 2317) + # 12
  geom_vline(xintercept = 2420) + # 11
  geom_vline(xintercept = 2511) + # 11  
  geom_vline(xintercept = 2593) + # 10
  geom_vline(xintercept = 2672) + # 9
  geom_vline(xintercept = 2732) + # 8
  geom_vline(xintercept = 2796) + # 6
  geom_vline(xintercept = 2845) + # 7
  geom_vline(xintercept = 2897) + # 5
  ggtitle(paste0("TCGA ID: ",rownames(output[id_num,]))) +
  xlab("Genome Position (1 Mb bin)")+ ylab("Mean segVal")+
  theme_minimal() + theme(axis.text.x = element_text(size = 5,angle = 90, hjust = 0.1))

p2
     
#TCGA_BLCA_CIN_measures <- signatures_matrix_df


#ggsave(ggarrange(p1,p2,ncol = 1,nrow = 2), filename = "~/Desktop/test.png",width = 20, height = 10)
output_pancan$sample_id <- rownames(output_pancan)
subsample_with_signatures <- merge(TCGA_BLCA_CIN_measures, output_pancan, by = "sample_id")

#both_signatures <- both_signatures #%>% filter(ajcc_pathologic_tumor_stage %in% late_stages)

output_combine <- combineSignatureWithChrSpec(subsample_with_signatures, segValColNames = colnames(output_pancan))


total_canc_spec_loss <- data.frame()
total_canc_spec_gain <- data.frame()
for(t in unique(subsample_with_signatures$type)){
  tt <- subsample_with_signatures %>% filter(type == t)
  output_combine <- combineSignatureWithChrSpec(tt, segValColNames = colnames(output_pancan))
  # loss
  canc_spec_loss <- t(output_combine$lossMatrix)
  rownames(canc_spec_loss) <- paste0(t,"_",rownames(canc_spec_loss) )
  total_canc_spec_loss <- rbind(total_canc_spec_loss, canc_spec_loss)
  # gains
  canc_spec_gain <- t(output_combine$gainMatrix)
  rownames(canc_spec_gain) <- paste0(t,"_",rownames(canc_spec_gain) )
  total_canc_spec_gain <- rbind(total_canc_spec_gain, canc_spec_gain)
  
}





ann <- data.frame(Chromosome = as.character(unlist(lapply(rownames(output_combine$lossMatrix), function(x){
  split<-str_split(x,pattern = "_")
  substr(split[[1]][2],start = 4,5)
}))))
rownames(ann) <- rownames(output_combine$lossMatrix)
ann_color = c("#FF0000", "#FF7F00", "#FFD400", "#FFFF00","#BFFF00","#6AFF00","#00EAFF","#0095FF", "#0040FF", "#AA00FF","#FF00AA","#EDB9B9", "#E7E9B9","#B9EDE0","#B9D7ED","#DCB9ED","#8F2323","#8F6A23","#4F8F23","#23628F","#000000","#737373")
names(ann_color) <- as.character(1:22)
anno_colors=list(Chromosome=ann_color)
ann_row <- data.frame(Cancer_Type = as.character(unlist(lapply(rownames(output_combine$lossMatrix), function(x){
  split<-str_split(x,pattern = "_")
  substr(split[[1]][1],start = 1,5)
}))))
rownames(ann_row) <- rownames(output_combine$lossMatrix)
heatmap_luad_loss <- pheatmap(abs(t(output_combine$lossMatrix)), 
                           clustering_method = "ward.D2", 
                           show_colnames = F,
                           show_rownames = T,
                           cluster_rows = F,
                           annotation_col = ann,
                           #annotation_row = ann_row,
                           annotation_colors = anno_colors,
                           cluster_cols = F,
                           scale = 'row',
                           color = c("#0115C7","#2C3CC7","#C7C60E","#C7B713","#C79E12","#CF8611","#B85A18", "#CF3711","#C41019", "#FF0000"),
                           breaks = seq(from = 0, to = 4, by = 0.3),
                           cutree_cols = 3,
                           fontsize = 10,
                           main = "TCGA LUAD - Loss", annotation_legend = F)


z <- do.call(gridExtra::grid.arrange, list( heatmap_luad_gain[[4]], heatmap_ov_gain[[4]], heatmap_luad_loss[[4]] , heatmap_ov_loss[[4]]))
plot(z)
ggsave(z,filename = "~/Desktop/blca_lusc.jpg", width = 12, height = 5)

f_plot <- reshape2::melt(abs(t(output_combine$lossMatrix)))
ggplot(f_plot, aes(x = Var2, y = Var1, size =value)) +
  geom_point()

library(umap)
umap_signatures <- umap(total_canc_spec_gain)
for_plot <- as.data.frame(umap_signatures$layout)
for_plot$type <- as.character(unlist(lapply(rownames(for_plot), function(x){
  split<-str_split(x,pattern = "_")
  substr(split[[1]][1],start = 1,5)
})))

colors <- c(ann_color, "#C133FF", "#3390FF", "#33FFA5")
names(colors) <- NA
p1 <- ggplot(for_plot, aes(x = V1, y=V2, color = type)) +
  geom_point() + xlab("UMAP Dim 1") + ylab("UMAP Dim 2") +
  geom_label(aes(label=type), label.size=0.1 ) + theme_minimal() + ggtitle("Pan-Cancer Copy Number Signatures: Gains ") + 
  theme(legend.position = "none") + theme(text = element_text(size=20)) 

ggsave(ggarrange(p1,p2), filename = "~/Desktop/umap.png", width = 20, height = 10 )

 #multiply_luad <- multiply
pc <- prcomp(both_signatures[,which(colnames(both_signatures) == "X1"):which(colnames(both_signatures) == "X14")], scale. = T, center = T)
summary(pc)
#library(ggbiplot)
# PC plot of signatures

set.seed(7)
k_means <- kmeans(x = pc$x[,1:6], 3, iter.max = 10000)
ggbiplot::ggbiplot(pc, ellipse=F, groups=as.factor(k_means$cluster)) + theme_minimal()
both_signatures$pca_kmeans_sigs <- k_means$cluster

surv_os <- Surv(as.numeric(as.character(both_signatures[,'OS.time']))/365, as.numeric(as.character(both_signatures[,'OS'])))
fit_os <- survfit(survplot::censor(surv_os, 5)~pca_kmeans_sigs, data = both_signatures)
p <- makeSurvPlot(fit_os, 
                  'TCGA Survival Analysis - CIN High Clustering', 
                  legend_labs = c("1","2","3"),
                  legen_title = "CIN", 
                  ylab = "OS")
p


# 
# tmp <- output_combine$meltedMergedGainLoss
# tmp$GainLoss <- substr(tmp$Var1,0,4)
# tmp$position <- substr(tmp$Var1,6,20)
# tmp$Chr <- unlist(lapply(tmp$Var1, function(x){
#   split<-str_split(x,pattern = "_")
#   substr(split[[1]][2],start = 4,5)
# }))
# p2 <- ggplot(tmp, aes(x = position, y = value, fill = Chr)) + 
#   geom_col() + ggtitle("LUSC") +
#   facet_wrap(~Var2, nrow = 2,ncol = 5) + scale_y_log10() + 
#   theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank())
# 
# ggsave(ggarrange(p1, p2), filename = "~/Desktop/gain_losses_blca_lusc.png", width = 20, height = 10)

