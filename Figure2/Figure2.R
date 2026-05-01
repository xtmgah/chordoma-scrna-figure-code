#Figure 2
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggsci)
library(tibble)
library(survival)
library(survminer)
library(cluster)

#Fig2a
# devtools::install_github("sqjin/CellChat")
library(Seurat)
library(CellChat)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
setwd("./sc/Analysis_CellChat/Analysis_CellChat/")#.
cellchat_object=readRDS("./sc/cellchat_object.rds")
zzm60colors1 <- c(
  '#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
  '#6d6fa0','#8d689d','#c8c7e1','#d25774',
  '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
  '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
  '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
  '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462'
)
cell_type_order <- c(
  "Chordoma", "Macrophage", "T cells", "Fibroblast", "Neutrophil", "Cycling",
  "B cells", "plasmablasts", "Endothelial", "Mast cells", "Osteoclast", "DCs",
  "Mural cells", "Plasma cell"
)

existing_types <- intersect(cell_type_order, levels(cellchat_object@idents))
cellchat_colors <- setNames(
  zzm60colors1[1:length(existing_types)],
  existing_types
)

groupSize <- as.numeric(table(cellchat_object@idents))
pdf("./cellchat_Num_all.pdf", height = 7, width = 12)
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat_object@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions",
                 color.use = cellchat_colors)

netVisual_circle(cellchat_object@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Interaction weights/strength",
                 color.use = cellchat_colors)

dev.off()

##Fig2e----
output_dir <- "./pathway_analysis_all/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("./pathway_analysis_all/")

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("./AllSignalHeatmap.pdf", width = 8, height = 6)
p1 <- netAnalysis_signalingRole_heatmap(cellchat_object, pattern = "outgoing",font.size = 5,width = 6,height = 9, color.use = cellchat_colors)
p2 <- netAnalysis_signalingRole_heatmap(cellchat_object, pattern = "incoming",font.size = 5,width = 6,height = 9, color.use = cellchat_colors)
p1+p2
dev.off()

##Fig2a2----
#PC cellcommu----
output_dir <- "./sc/Analysis_CellChat/pathway_analysis_PC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(".sc/Analysis_CellChat/pathway_analysis_PC/")

cellchat_PC=readRDS("./sc/cellchat_PC.rds")
pdf("./cellchat_Num_PC.pdf", height = 7, width = 12)
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat_PC@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions",
                 color.use = cellchat_colors)

netVisual_circle(cellchat_PC@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Interaction weights/strength",
                 color.use = cellchat_colors)

dev.off()


###Fig2b----
mat <- cellchat_PC@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  pdf(paste0(rownames(mat)[i], "_outgoing_signals_PC.pdf"), width = 6, height = 8)
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   edge.weight.max = max(mat), 
                   title.name = paste0(rownames(mat)[i], " -> others"),
                   color.use = cellchat_colors) 
  dev.off()
}

# Access all the signaling pathways showing significant communications----
##Fig2g----
pathways.show.all <- cellchat_PC@netP$pathways
levels(cellchat_PC@idents)
vertex.receiver = c(1:3) 
for (i in 1:length(pathways.show.all)) {
  pathway <- pathways.show.all[i]

  tryCatch({
    pdf(paste0(pathway, "_hierarchy_PC.pdf"), width = 8, height = 6)
    netVisual_aggregate(cellchat_PC, 
                        signaling = pathway,
                        vertex.receiver = vertex.receiver,
                        layout = "hierarchy",
                        color.use = cellchat_colors)
    dev.off()
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    while (dev.cur() > 1) {
      dev.off()
    }
  })
}


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores


####CC cellcommu----
library(CellChat)
options(stringsAsFactors = FALSE)
cellchat_CC<- readRDS("F:/Chordoma/Plots/cellchat_CC.rds")

output_dir <- "F:/Chordoma/Plots/Analysis_CellChat/pathway_analysis_CC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd("F:/Chordoma/Plots/Analysis_CellChat/pathway_analysis_CC/")
#Fig2b1----
pdf("./cellchat_Num_CC.pdf", height = 7, width = 12)
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat_CC@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Number of interactions",
                 color.use = cellchat_colors)

netVisual_circle(cellchat_CC@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = TRUE, 
                 label.edge = FALSE, 
                 title.name = "Interaction weights/strength",
                 color.use = cellchat_colors)

dev.off()

mat <- cellchat_CC@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  pdf(paste0(rownames(mat)[i], "_outgoing_signals_CC.pdf"), width = 6, height = 8)
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = TRUE, 
                   edge.weight.max = max(mat), 
                   title.name = paste0(rownames(mat)[i], " -> others"),
                   color.use = cellchat_colors)  
  dev.off()
}


# Access all the signaling pathways showing significant communications----
##Fig2f----
pathways.show.all <- cellchat_CC@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat_CC@idents)
vertex.receiver = c(1:3) 

for (i in 1:length(pathways.show.all)) {
  pathway <- pathways.show.all[i]
  
  tryCatch({
    pdf(paste0(pathway, "_hierarchy.pdf"), width = 8, height = 6)
    netVisual_aggregate(cellchat_CC, 
                        signaling = pathway,
                        vertex.receiver = vertex.receiver,
                        layout = "hierarchy",
                        color.use = cellchat_colors)
    dev.off()
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    while (dev.cur() > 1) {
      dev.off()
    }
  })
}


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
##FigS4
cellchat_CC<- netAnalysis_computeCentrality(cellchat_CC, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("/home/data2/zhangmy/sc/Analysis_CellChat/pathway_analysis_CC/AllSignalHeatmap_CC.pdf", width = 8, height = 6)
p1 <- netAnalysis_signalingRole_heatmap(cellchat_CC, pattern = "outgoing",font.size = 5,width = 6,height = 9, color.use = cellchat_colors)
p2 <- netAnalysis_signalingRole_heatmap(cellchat_CC, pattern = "incoming",font.size = 5,width = 6,height = 9, color.use = cellchat_colors)
p1+p2
dev.off()


object.list <- list(PC = cellchat_PC, CC = cellchat_CC)
cellchat_sample <- mergeCellChat(object.list, add.names = names(object.list)) 
pathway.show <- intersect(cellchat_PC@netP$pathways, cellchat_CC@netP$pathways)
patterns <- c("incoming", "outgoing", "all")
pathway.use <- c("MIF", "FN1","APP")

for (pattern in patterns) {
  ht_pc <- netAnalysis_signalingRole_heatmap(cellchat_PC, 
                                             pattern = pattern, 
                                             signaling = pathway.use, 
                                             title = paste("PC -", pattern),
                                             width = 5, 
                                             height = 6,
                                             color.heatmap = "Purples")
  
  ht_cc <- netAnalysis_signalingRole_heatmap(cellchat_CC, 
                                             pattern = pattern, 
                                             signaling = pathway.use, 
                                             title = paste("CC -", pattern),
                                             width = 5, 
                                             height = 6,
                                             color.heatmap = "Purples")
  
  
  pdf(paste0("/home/data2/zhangmy/sc/Analysis_CellChat/signaling_", pattern, "_comparison.pdf"), 
      width = 10, height = 8)
  draw(ht_pc + ht_cc,column_gap = unit(15, "mm"))
  dev.off()
  
}



#figure2c
sc_combined_clean=readRDS("./sc/sc_combined_clean_clinical_full.rds")
prop_table <- table(sc_combined_clean$cell_type, sc_combined_clean$orig.ident) %>%
  prop.table(margin = 2)
prop_matrix <- as.matrix(prop_table)
key_cells <- c("Chordoma", "Macrophage", "T cells")
prop_key <- prop_matrix[key_cells, ] %>% t()  
dist_mat <- dist(prop_key, method = "euclidean")
hclust_res <- hclust(dist_mat, method = "ward.D2")

sil_scores <- sapply(2:min(6, nrow(prop_key)-1), function(k) {
  clusters <- cutree(hclust_res, k)
  mean(silhouette(clusters, dist_mat)[, 3])
})
best_k <- which.max(sil_scores) + 1


sample_clusters <- data.frame(
  Sample = rownames(prop_key),
  Cluster = factor(cutree(hclust_res, k = 3))
)

cluster_groups <- cutree(hclust_res, k = 3) 
sample_clusters <- data.frame(
  Sample = names(cluster_groups),
  Cluster = factor(cluster_groups)  
)

anno_df <- sc_combined_clean@meta.data %>%
  select(orig.ident, clinical_group, gender, location, status, OS_event,PFS_event) %>%
  distinct() %>%
  remove_rownames() %>%
  column_to_rownames("orig.ident")

anno_df$Cluster <- sample_clusters$Cluster[match(rownames(anno_df), sample_clusters$Sample)]

library(RColorBrewer)

cluster_colors <- brewer.pal(3, "Set1")
names(cluster_colors) <- as.character(1:3)

ann_colors <- list(
  clinical_group = c("CC" = "#61bada", "PC" = "#72567a"),
  gender = c("female" = "#F4A582", "male" = "#92C5DE"),
  location = c("clivus" = "#66C2A5", "sacrum" = "#FC8D62", "mobile spine" = "#8DA0CB"),
  status = c("primary" = "#A6D854", "recurrence" = "#E78AC3"),
  OS_event = c("0" = "grey90", "1" = "black","/"="white"),
  PFS_event = c("0" = "grey90", "1" = "red","/"="white"),
  Cluster = cluster_colors
)

pdf("./sc/heatmap.pdf",width=12,height=9)
pheatmap(t(prop_key),  
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         annotation_col = anno_df,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,  # 细胞不聚类
         cluster_cols = hclust_res,  # 样本聚类
         show_colnames = TRUE,
         main = paste("Clustering based on", paste(key_cells, collapse = ", "), "(K =", 3, ")"))
dev.off()


##Fig2d----
clinical_df <- read.csv("./sc/clinical.csv", header=TRUE, stringsAsFactors = FALSE)
clinical_df=clinical_df[,1:13]

surv_df <- clinical_df %>%
  left_join(sample_clusters, by = "Sample") %>%
  mutate(
    OS_time = suppressWarnings(as.numeric(OS.Month.)),
    OS_status = suppressWarnings(as.numeric(as.character(OS_event))),
    PFS_time = suppressWarnings(as.numeric(PFS.Month.)),
    PFS_status = suppressWarnings(as.numeric(as.character(PFS_event))),
    age = suppressWarnings(as.numeric(age))
  ) %>%
  filter(!is.na(Cluster))

surv_df$Cluster <- as.factor(surv_df$Cluster)

cluster_levels <- levels(surv_df$Cluster)
cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
cluster_colors <- cluster_colors[seq_along(cluster_levels)]
names(cluster_colors) <- cluster_levels

age_median <- median(surv_df$age, na.rm = TRUE)

# -----------------------------
# OS: age-adjusted survival
# -----------------------------
os_df <- surv_df %>%
  filter(!is.na(OS_time), !is.na(OS_status), !is.na(age), !is.na(Cluster))

if (nrow(os_df) > 0 && length(unique(os_df$Cluster)) >= 2) {
  os_df$Cluster <- droplevels(os_df$Cluster)
  
  cox_os <- coxph(Surv(OS_time, OS_status) ~ Cluster + age, data = os_df)
  cox_os_null <- coxph(Surv(OS_time, OS_status) ~ age, data = os_df)
  
  os_compare <- anova(cox_os_null, cox_os, test = "Chisq")
  os_p <- NA_real_
  if (nrow(os_compare) >= 2) {
    p_col <- grep("P", colnames(os_compare), value = TRUE)
    if (length(p_col) > 0) {
      os_p <- suppressWarnings(as.numeric(os_compare[2, p_col[1]]))
    }
  }
  
  os_p_label <- if (is.na(os_p)) {
    "Adjusted Cox P = NA"
  } else if (os_p < 0.001) {
    "Adjusted Cox P < 0.001"
  } else {
    paste0("Adjusted Cox P = ", formatC(os_p, format = "f", digits = 3))
  }
  
  newdata_os <- data.frame(
    age = rep(age_median, length(levels(os_df$Cluster))),
    Cluster = factor(levels(os_df$Cluster), levels = levels(os_df$Cluster))
  )
  
  os_fit_adj <- survfit(cox_os, newdata = newdata_os)
  
  p_os <- ggsurvplot(
    os_fit_adj,
    data = newdata_os,
    conf.int = FALSE,
    risk.table = FALSE,
    censor = FALSE,
    break.time.by = 12,
    palette = cluster_colors[levels(os_df$Cluster)],
    title = paste0("Overall Survival by Cluster (age-adjusted, age = ", round(age_median, 1), ")"),
    xlab = "Time (months)",
    ylab = "Overall survival probability",
    legend.title = "Cluster",
    legend.labs = levels(os_df$Cluster),
    ggtheme = theme_bw(base_size = 14)
  )
  
  max_os_time <- max(os_df$OS_time, na.rm = TRUE)
  
  p_os$plot <- p_os$plot +
    annotate("text", x = max_os_time * 0.60, y = 0.15, label = os_p_label, size = 5) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.82, 0.82),
      legend.background = element_rect(fill = "white", color = "grey70")
    )
  
  pdf("./sc/OS_3cluster_age_adjusted.pdf", width = 12, height = 9)
  print(p_os)
  dev.off()
}

# -----------------------------
# PFS: age-adjusted survival
# -----------------------------
pfs_df <- surv_df %>%
  filter(!is.na(PFS_time), !is.na(PFS_status), !is.na(age), !is.na(Cluster))

if (nrow(pfs_df) > 0 && length(unique(pfs_df$Cluster)) >= 2) {
  pfs_df$Cluster <- droplevels(pfs_df$Cluster)
  
  cox_pfs <- coxph(Surv(PFS_time, PFS_status) ~ Cluster + age, data = pfs_df)
  cox_pfs_null <- coxph(Surv(PFS_time, PFS_status) ~ age, data = pfs_df)
  
  pfs_compare <- anova(cox_pfs_null, cox_pfs, test = "Chisq")
  pfs_p <- NA_real_
  if (nrow(pfs_compare) >= 2) {
    p_col <- grep("P", colnames(pfs_compare), value = TRUE)
    if (length(p_col) > 0) {
      pfs_p <- suppressWarnings(as.numeric(pfs_compare[2, p_col[1]]))
    }
  }
  
  pfs_p_label <- if (is.na(pfs_p)) {
    "Adjusted Cox P = NA"
  } else if (pfs_p < 0.001) {
    "Adjusted Cox P < 0.001"
  } else {
    paste0("Adjusted Cox P = ", formatC(pfs_p, format = "f", digits = 3))
  }
  
  newdata_pfs <- data.frame(
    age = rep(age_median, length(levels(pfs_df$Cluster))),
    Cluster = factor(levels(pfs_df$Cluster), levels = levels(pfs_df$Cluster))
  )
  
  pfs_fit_adj <- survfit(cox_pfs, newdata = newdata_pfs)
  
  p_pfs <- ggsurvplot(
    pfs_fit_adj,
    data = newdata_pfs,
    conf.int = FALSE,
    risk.table = FALSE,
    censor = FALSE,
    break.time.by = 12,
    palette = cluster_colors[levels(pfs_df$Cluster)],
    title = paste0("Progression-Free Survival by Cluster (age-adjusted, age = ", round(age_median, 1), ")"),
    xlab = "Time (months)",
    ylab = "Progression-free survival probability",
    legend.title = "Cluster",
    legend.labs = levels(pfs_df$Cluster),
    ggtheme = theme_bw(base_size = 14)
  )
  
  max_pfs_time <- max(pfs_df$PFS_time, na.rm = TRUE)
  
  p_pfs$plot <- p_pfs$plot +
    annotate("text", x = max_pfs_time * 0.60, y = 0.15, label = pfs_p_label, size = 5) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.82, 0.82),
      legend.background = element_rect(fill = "white", color = "grey70")
    )
  
  pdf("./sc/PFS_3cluster_age_adjusted.pdf", width = 12, height = 9)
  print(p_pfs)
  dev.off()
}

