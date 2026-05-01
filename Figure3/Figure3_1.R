library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Seurat)
zzm60colors2 <- c(
  '#76a2be','#4b6aa8','#c6adb0','#df5734','#6c408e',
  '#ac6894','#b7deea','#83ab8e','#d4c2db','#ece399',
  '#cbdaa9','#b95055','#2d3462','#bc9a7f','#e0cfda',
  '#e6b884','#d69a55','#64a776','#cc7f73','#927c9a',
  '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
  '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
  '#c49abc','#b05545','#405993','#9f8d89','#72567a',
  '#63a3b8','#c4daec','#3674a2','#537eb7','#e29eaf',
  '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
  '#408444','#9d3b62','#d5bb72','#d8a0c0','#61bada'
)

###subcluster tumor cells----
sc_tumor=readRDS("F:/Chordoma/Result/sc_tumor.rds")

tumor_meta <- sc_tumor@meta.data
cluster_col <- "RNA_snn_res.0.5" 
#print(paste("sc_tumor cell number:", ncol(sc_tumor)))

print(table(sc_tumor$clinical_group))
Idents(sc_tumor) <- cluster_col
n_clusters <- length(unique(Idents(sc_tumor)))
print(paste( n_clusters, "subclusters"))


p_umap_final <- DimPlot(sc_tumor, 
                        reduction = "umap", 
                        label = TRUE, 
                        label.size = 5, 
                        repel = TRUE, 
                        pt.size = 0.8, # 点大一点，看清楚分布
                        cols = zzm60colors2, 
                        raster = FALSE) +
  ggtitle(paste0("Tumor Sub-clustering (Res ", target_res, ")")) +
  theme(legend.position = "right")

# 保存
ggsave(paste0("F:/Chordoma/Plots/Tumor_Final_UMAP_Res", target_res, ".pdf"), 
       p_umap_final, width = 10, height = 8)

##barplot in cc/pc
library(ggplot2)
library(dplyr)
library(ggsci)

props <- sc_tumor@meta.data %>%
  group_by(!!sym(target_col), clinical_group) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# 画图
p_check <- ggplot(props, aes(x = !!sym(target_col), y = freq, fill = clinical_group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_fill_manual(values = c("CC"="#61bada", "PC"="#72567a")) + 
  theme_classic() +
  labs(title = paste0("Clinical Composition (Res ", target_res, ")"), 
       subtitle = "Look for clusters that are mostly PURPLE (PC) or BLUE (CC)",
       x = "Cluster ID", y = "Proportion") +
  geom_hline(yintercept = 0.5, linetype="dashed", color = "grey") +
  theme(axis.text.x = element_text(angle = 0, size = 10, face = "bold"))

# 保存并打印
ggsave(paste0("F:/Chordoma/Plots/Check_Res_", target_res, "_Composition.pdf"), 
       p_check, width = 12, height = 6)



###Tumor subtypes in 3 clusters plots###
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

props <- sc_tumor@meta.data %>%
  group_by(Cluster,RNA_snn_res.0.5) %>% # 按病人组和肿瘤亚群分组
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Cluster) %>%
  mutate(freq = n / sum(n))

p_bar <- ggplot(props, aes(x = Cluster, y = freq, fill = RNA_snn_res.0.5)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = zzm60colors2) +
  theme_classic() +
  labs(title = "Tumor Subtype Composition by Patient Cluster", 
       y = "Proportion", fill = "Tumor Subtype") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold"))

print(p_bar)
out_dir <- "F:/Chordoma/Plots/"
ggsave(file.path(out_dir, "PatientCluster_TumorComposition_Barplot.pdf"), p_bar, width = 8, height = 7)




# ==============================================================================
# 2. Butterfly Plot
# ==============================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel) 
deg_auc=read.csv("F:\\Chordoma\\Result\\DEGs_PCvsCC_AUC_Result.csv",header=TRUE,row.names=1)
fc_col <- ifelse("avg_log2FC" %in% colnames(deg_auc), "avg_log2FC", "avg_diff")
#  Top gene (high AUC)
top_pc <- deg_auc %>% filter(Group == "PC-Marker") %>% top_n(10, wt = final_AUC) %>% pull(gene)
top_cc <- deg_auc %>% filter(Group == "CC-Marker") %>% top_n(10, wt = final_AUC) %>% pull(gene)
top_fc_pc <- deg_auc %>% arrange(desc(!!sym(fc_col))) %>% head(2) %>% pull(gene)
top_fc_cc <- deg_auc %>% arrange(!!sym(fc_col)) %>% head(2) %>% pull(gene)

manual_genes <- c("TWIST1", "PTGES", "MMP9", "S100B", "TBXT", "KRT19")

labels_to_show <- unique(c(top_pc, top_cc, top_fc_pc, top_fc_cc, manual_genes))

print(length(labels_to_show))


labels_data <- deg_auc %>% filter(gene %in% labels_to_show)

p_scatter <- ggplot(deg_auc, aes(x = avg_log2FC, y = myAUC)) +
  geom_point(aes(color = myAUC, size = myAUC), alpha = 0.8) +
  scale_color_gradientn(colors = c("grey", "#FFEDA0", "#FEB24C", "#F03B20")) +
  scale_size(range = c(1, 4)) +
  
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  geom_text_repel(data = subset(deg_auc, gene %in% labels_to_show),
                  aes(label = gene),
                  size = 4, fontface = "bold",
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  

  theme_classic() +
  labs(title = "Marker Power: PC vs CC (AUC Test)",
       subtitle = "Y-axis: Predictive Power (AUC) | X-axis: Fold Change",
       x = "Log2 Fold Change",
       y = "Predictive Power (myAUC)") +
  theme(legend.position = "right")

print(p_scatter)

ggsave("F:\\Chordoma\\Plots\\butterfly_AUC_PC_vs_CC.pdf", 
       p_scatter, width = 10, height = 8)


####
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

pc_genes_122 <- deg_auc %>% 
  filter(myAUC > 0.7 & avg_log2FC > 0.25) %>% 
  pull(gene)

print(paste("提取到 PC 特异基因数量:", length(pc_genes_122)))

# ==============================================================================
# 2.  10  CC markers 
# ==============================================================================

cc_genes_10 <- deg_auc %>% 
  filter(myAUC < 0.3 & avg_log2FC < -0.25) %>% 
  pull(gene)

print(cc_genes_10) 


# ==============================================================================
# 3.  PC  122 markers_GO
# ==============================================================================
gene_convert <- bitr(pc_genes_122, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

if(nrow(gene_convert) > 0) {
  ego_pc <- enrichGO(gene = gene_convert$ENTREZID, 
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     readable = TRUE)
  
  if(!is.null(ego_pc)) {
    p_go <- dotplot(ego_pc, showCategory=20, title="Functions Gained in PC Group (122 Genes)")
    ggsave("F:\\Chordoma\\Plots\\GO_PC_122_Genes.pdf", p_go, width = 10, height = 16)
    print(p_go)
  }
}

library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)

Idents(sc_tumor) <- "tumor_functional_type" 
sc_pc_only <- subset(sc_tumor, subset = clinical_group == "PC")
sc_pc_only$tumor_functional_type <- Idents(sc_pc_only)

df <- FetchData(sc_pc_only, vars = c("TWIST1", "PTGES"))
df <- df[df$TWIST1 > 0 | df$PTGES > 0, ]


#TWIST1 cor PTGES boxplot----

DefaultAssay(sc_pc_only) <- "RNA"
expr_data <- FetchData(sc_pc_only, vars = c("TWIST1", "PTGES"))

# 2.  TWIST1 as Low, Mid, High 
#  (Quantile) group：0-33%, 33-66%, 66-100%
breaks <- quantile(expr_data$TWIST1[expr_data$TWIST1 > 0], probs = c(0.33, 0.66))

expr_data$Group <- case_when(
  expr_data$TWIST1 == 0 ~ "Neg",
  expr_data$TWIST1 <= breaks[1] ~ "Low",
  expr_data$TWIST1 <= breaks[2] ~ "Mid",
  TRUE ~ "High"
)


expr_data$Group <- factor(expr_data$Group, levels = c("Neg", "Low", "Mid", "High"))


p_box_trend <- ggplot(expr_data, aes(x = Group, y = PTGES, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.1, width = 0.2, size = 0.5) + 
  stat_compare_means(method = "kruskal.test", label.y = max(expr_data$PTGES) + 0.5) + # 总体差异P值
  scale_fill_brewer(palette = "Reds") +
  theme_classic() +
  labs(title = "PTGES increases with TWIST1 levels",
       x = "TWIST1 Expression Level", 
       y = "PTGES Expression")


ggsave("F:/Chordoma/Plots/TWIST1_PTGES_Trend_Boxplot.pdf", p_box_trend, width = 6, height = 6)
print(p_box_trend)


#########CC vs PC的co-correlation----
library(Seurat)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(corrr)
#saveRDS(graph_cc,"F:/Chordoma/Plots/graph_cc.rds")
graph_cc=readRDS("F:/Chordoma/Plots/graph_cc.rds")
# ==============================================================================
# plot cc net
# ==============================================================================
p_net_cc <- ggraph(graph_cc, layout = "nicely") + 
  
  geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight), color = weight), 
                 show.legend = TRUE) +
  
  scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  
  scale_edge_width(range = c(0.2, 1)) +
  
  geom_node_point(aes(color = type, size = degree), alpha = 0.9) +

  scale_color_manual(values = c("Core (Top 10)" = "#E41A1C", "Neighbor" = "#A6CEE3")) +
  
  scale_size(range = c(3, 8)) +
  
  geom_node_text(aes(label = ifelse(type == "Core (Top 10)" | degree > 0, name, "")), 
                 repel = TRUE, size = 3.5, fontface = "bold") +
  
  theme_void() +
  labs(title = "Co-expression Network in CC Tumor Cells",
       subtitle = paste0("Core Genes: CC-Specific Top 10 | Threshold: r > ", cor_threshold),
       color = "Gene Type")


ggsave("F:\\Chordoma\\Plots\\Network_CC_Top10_Neighbors.pdf", 
       p_net_cc, width = 12, height = 10)

#write.csv(edge_list,"F:\\Chordoma\\Result\\Network_CC_Neighbors.csv",quote = FALSE)
#write.csv(target_edges,"F:\\Chordoma\\Result\\Network_CC_Top10_Neighbors_target_edges.csv",quote = FALSE)
print(p_net_cc)

##PC co-correlation net----
library(Seurat)
library(tidyverse)
library(igraph)
library(ggraph)
library(corrr)
#saveRDS(graph_pc,"F:/Chordoma/Plots/graph_pc.rds")
graph_pc=readRDS("F:/Chordoma/Plots/graph_pc.rds")
p_net_pc <- ggraph(graph_pc, layout = "nicely") + 
  
  geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight), color = weight), 
                 show.legend = TRUE) +
  scale_edge_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_edge_width(range = c(0.2, 1)) +
  
  geom_node_point(aes(color = type, size = degree), alpha = 0.9) +
  
  scale_color_manual(values = c("Core (Top 10)" = "red", "Neighbor" = "#c5b0d5")) +
  
  scale_size(range = c(3, 8)) +

  geom_node_text(aes(label = ifelse(type == "Core (Top 10)" | degree > 0, name, "")), 
                 repel = TRUE, size = 3.5, fontface = "bold") +
  
  theme_void() +
  labs(title = "Co-expression Network in PC Tumor Cells",
       subtitle = paste0("Core Genes: PC-Specific Top 10 | Threshold: r > ", cor_threshold),
       color = "Gene Type")

print(p_net_pc)

ggsave("F:\\Chordoma\\Plots\\Network_PC_Top10_Neighbors.pdf", p_net, width = 12, height = 10)
#write.csv(edge_list, "F:\\Chordoma\\Result\\Tumor_PC\\Network_PC_All_Edges_Passed_Threshold.csv", quote = FALSE, row.names = FALSE)
#write.csv(target_edges, "F:\\Chordoma\\Result\\Tumor_PC\\Network_PC_Top10_Neighbors_Edges.csv", quote = FALSE, row.names = FALSE)

####Pseudotime analysis in tumor cells----
#devtools::install_github('cole-trapnell-lab/monocle3')
library(Seurat)
library(monocle3)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(patchwork)

cds=readRDS("F:\\Chordoma\\Plots\\cds.rds")
p1_custom <- plot_cells(cds, 
                        color_cells_by = target_col, 
                        label_cell_groups = TRUE, 
                        label_leaves = FALSE, 
                        label_branch_points = TRUE,
                        graph_label_size = 3,
                        cell_size = 0.8,
                        trajectory_graph_segment_size = 1.2) +
  
  scale_color_manual(values = zzm60colors2) +
  
  ggtitle(paste0("Trajectory (Root: Cluster ", root_cluster_id, ")")) +
  theme(legend.position = "right")


p2 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 cell_size = 0.8) +
  scale_color_viridis_c(option = "plasma") + 
  ggtitle("Pseudotime")


p_final_custom <- p1_custom | p2
print(p_final_custom)

ggsave("F:/Chordoma/Plots/Monocle3_Trajectory_Fixed_Color_Root12.pdf", 
       p_final_custom, width = 16, height = 7)



