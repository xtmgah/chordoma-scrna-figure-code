##Figure 1
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(rstatix)

sc_combined_clean=readRDS("./sc/sc_combined_clean.rds")
zzm60colors <- c(
  '#4b6aa8','#df5734','#6c408e','#ac6894','#d4c2db',
  '#b7deea','#83ab8e','#ece399','#61bada','#cc7f73',
  '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
  '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
  '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
  '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
  '#c49abc','#927c9a','#405993','#9f8d89','#72567a',
  '#63a3b8','#c4daec','#3674a2','#537eb7','#e29eaf',
  '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
  '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462'
)

target_resolution <- "RNA_snn_res.0.5" 
Idents(sc_combined_clean) <- target_resolution

p_cluster <- DimPlot(sc_combined_clean, 
                     reduction = "umap", 
                     label = TRUE,          
                     label.size = 5,        
                     repel = TRUE,          
                     pt.size = 0.1,         
                     cols = zzm60colors,    
                     raster = FALSE) +      
  ggtitle(paste0("UMAP Clustering (Resolution: ", target_resolution, ")")) +
  theme(legend.position = "right")    

print(p_cluster)

ggsave("./sc/clustering_resolution_0.5.pdf", 
       plot = combined_plot,
       width = 15, 
       height = 5,
       device = "pdf")



Idents(sc_combined_clean) <- "RNA_snn_res.0.5" # 设定默认身份

prop.table(table(Idents(sc_combined_clean), sc_combined_clean$orig.ident), margin = 1)
library(ggplot2)
library(dplyr) 
library(Seurat)

sample_count <- length(unique(sc_combined_clean$orig.ident))
if(sample_count > length(zzm60colors)){
  warning("样本数多于色盘颜色数，已自动扩展颜色。")
  zzm60colors <- colorRampPalette(zzm60colors)(sample_count)
}

cell_stats <- sc_combined_clean@meta.data %>%
  select(orig.ident, RNA_snn_res.0.5)

p_bar <- ggplot(cell_stats, aes(x = RNA_snn_res.0.5, fill = orig.ident)) +
  geom_bar(position = "fill", width = 0.8) +  # position="fill" 表示绘制百分比堆叠
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual(values = zzm60colors) +   # 应用你的自定义色盘
  theme_classic() +
  labs(
    title = "Sample Composition",
    x = "Cluster",
    y = "Proportion",
    fill = "Sample"
  ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "right" # 图例在右侧
  )

print(p_bar)
ggsave("./sc/p_bar.pdf", 
       plot = p_bar,
       width = 10, 
       height = 5,
       device = "pdf")

####
library(ggplot2)
zzm60colors1 <- c(
  '#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
  '#6d6fa0','#8d689d','#c8c7e1','#d25774',
  '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
  '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
  '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
  '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462'
)

p_anno <- DimPlot(sc_combined_clean, 
                  reduction = "umap", 
                  group.by = "cell_type",  
                  label = TRUE, 
                  label.size = 5,
                  repel = TRUE,
                  cols = zzm60colors1, 
                  raster = FALSE) + 
  ggtitle("Annotated Chordoma Atlas") +
  theme(legend.position = "right")

print(p_anno)

ggsave("./sc/anno_plot.pdf", 
       plot = p_anno,
       width = 15, 
       height = 15,
       device = "pdf")

p_sample <- DimPlot(sc_combined_clean, 
                    reduction = "umap", 
                    group.by = "orig.ident",  # 按样本着色
                    cols = zzm60colors,       # 使用你的色盘
                    pt.size = 0.05, 
                    shuffle = TRUE,           # 打乱点顺序，防止遮盖
                    raster = FALSE) + 
  ggtitle("Samples Distribution") +
  theme(legend.position = "right")

print(p_sample)
ggsave("./sc/p_batch.pdf", 
       plot = p_sample,
       width = 5, 
       height = 5,
       device = "pdf")

props <- sc_combined_clean@meta.data %>%
  group_by(orig.ident, cell_type) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n))

p_bar <- ggplot(props, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = zzm60colors1) +  # 使用你的色盘
  theme_classic() +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type", 
       title = "Cell Type Proportion per Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # X轴标签倾斜

print(p_bar)
ggsave("./sc/anno_bar_plot.pdf", 
       plot = p_bar,
       width = 15, 
       height = 15,
       device = "pdf")

# =======================================================
# 1. dotplot 
# =======================================================
level_order <- c(
  "Chordoma", 
  "Cycling",          # 增殖通常跟肿瘤放一起或单独放
  "Macrophage",
  "Osteoclast",       # 破骨跟巨噬系近
  "DCs",              # 髓系
  "Neutrophil",
  "Mast cells",
  "T cells", 
  "B cells", 
  "plasmablasts",     # B细胞谱系放一起
  "Plasma cell",
  "Fibroblast",
  "Endothelial",
  "Mural cells"
)

Idents(sc_combined_clean) <- factor(Idents(sc_combined_clean), levels = level_order)

markers_to_plot <- c(
  # 1. Chordoma
  "TBXT", "KRT19", "CD24",
  # 2. Cycling
  "MKI67", "TOP2A",
  # 3. Macrophage
  "CD68", "CD163", "C1QA","CD86",
  # 4. Osteoclast
  "CTSK", "ACP5",
  # 5. DCs
  "CLEC9A",
  # 6. Neutrophil
  "S100A8",
  # 7. Mast cells
  "TPSAB1",
  # 8. T cells
  "CD3D", "CD3E", 
  # 9. B cells 
  "MS4A1", "CD79A", 
  # 10. Plasma/plasmablasts 特有
  "JCHAIN", 
  # 11. Fibroblast
  "COL1A1", "DCN",
  # 12. Endothelial
  "PECAM1", "VWF",
  # 13. Mural cells
  "ACTA2"
)

# =======================================================
#  DotPlot
# =======================================================
sc_combined_clean=read.RDS("./sc/sc_combined_data.rds")
p_dot <- DotPlot(sc_combined_clean, 
                 features = markers_to_plot, 
                 cols = c("lightgrey", "#9d3b62"), 
                 dot.scale = 8,                   
                 cluster.idents = FALSE) +        
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold.italic", color = "black"), 
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(title = "Cell Type Markers")

print(p_dot)

ggsave("./sc/anno_dot_plot_fixed.pdf", 
       plot = p_dot,
       width = 14,  
       height = 8,   
       device = "pdf")

####Feature plot##
library(Seurat)
library(ggplot2)
library(patchwork) 

# =======================================================
# Marker 
# =======================================================

feats_to_plot <- c(
  "TBXT",    # Chordoma 
  "MKI67",   # Cycling 
  "COL1A1",  # Fibroblast 
  
  "CD3D",    # T cells 
  "MS4A1",   # B cells
 
  "JCHAIN",  # Plasma cells 
  "CD163",   # Macrophage 
  "CTSK",    # Osteoclast
  
  "S100A8",  # Neutrophil 
  "TPSAB1",  # Mast cell 
  "VWF"   ,   # Endothelial ,
  "ACTA2" ##Mural cells
)

# =======================================================
#  FeaturePlot (3列 * 4行)
# =======================================================
p_feat <- FeaturePlot(sc_combined_clean, 
                      features = feats_to_plot,
                      cols = c("#f0f0f0", "#da6f6d"), 
                      order = TRUE,                   
                      min.cutoff = 'q1',              
                      ncol = 3,                      
                      raster = FALSE) &               
  NoAxes() &                                         
  theme(plot.title = element_text(size = 12, face = "bold")) 


ggsave("./sc/anno_featureplot.pdf",
       plot = p_feat, 
       width = 12,    # 宽度适中
       height = 16,   # 高度设大，因为有4行
       limitsize = FALSE) # 防止图片过大报错

#print(p_feat)
#saveRDS(sc_combined_clean, "F:\\Chordoma\\sc_combined_clean.rds") 

# UMAP (CC vs PC)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

clinical_data <- read.csv("F:/Chordoma/clinical.csv", header=TRUE, stringsAsFactors = FALSE)
clinical_data=clinical_data[,1:13]
print(head(clinical_data))
meta <- sc_combined_clean@meta.data
meta$barcode <- rownames(meta) 
meta_df_merged <- meta %>%
  left_join(clinical_data, by = c("orig.ident" = "Sample"))

plot_data_grouped <- meta_df_merged %>%
  filter(!is.na(Subtype)) %>%
  
  # 2. 【核心差异】只按 组别(type) 和 细胞类型(cell_type) 分组
  # 这里不再包含 orig.ident，相当于把同组的所有样本合并成一个“大样本”
  group_by(Subtype, cell_type) %>%
  
  # 3. 统计数量
  summarise(count = n(), .groups = 'drop') %>%
  
  # 4. 计算在该组(PC或CC)内的总占比
  group_by(Subtype) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p_grouped_bar <- ggplot(plot_data_grouped, aes(x = Subtype, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) + 
  
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = zzm60colors1) + 
  theme_classic() +
  labs(x = "Group", y = "Proportion", fill = "Cell Type", 
       title = "Cell Type Proportion: PC vs CC (Aggregated)") +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", color = "black"), # 加大X轴字体
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

print(p_grouped_bar)

ggsave("./sc/PC_vs_CC_overall_barplot.pdf",
       plot = p_grouped_bar, 
       width = 16,    
       height = 12,  
       limitsize = FALSE)


type_map <- setNames(clinical_data$Subtype, clinical_data$Sample)
group_vector <- type_map[as.character(sc_combined_clean$orig.ident)]
sc_combined_clean$clinical_group <- unname(group_vector)

print("分组统计：")
print(table(sc_combined_clean$clinical_group))
p_split <- DimPlot(sc_combined_clean, 
                   reduction = "umap", 
                   group.by = "cell_type",      # 细胞类型
                   split.by = "clinical_group", # 按 PC/CC 分面
                   cols = zzm60colors1,          
                   label = TRUE,                
                   label.size = 3,
                   repel = TRUE,                
                   pt.size = 0.1,               
                   raster = FALSE) +            
  ggtitle("Cell Type Distribution: PC vs CC") +
  theme(legend.position = "right")

print(p_split)
ggsave("./sc/PC_vs_CC_dimplot.pdf",
       plot = p_split, 
       width = 20,    
       height = 15,  
       limitsize = FALSE)

#####cell composition compare test
sample_props_complete=readRDS("F:/Chordoma/Plots/sample_props_complete.rds")

p_box <- ggplot(sample_props_complete, aes(x = clinical_group, y = freq, fill = clinical_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5) + 
  
  facet_wrap(~cell_type, scales = "free_y", ncol = 4) +

  stat_pvalue_manual(stat.test, 
                     "p.adj.label",  
                     tip.length = 0.01,
                     size = 3.5) +     
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c("#61bada", "#72567a")) + 
  theme_bw() +
  labs(y = "Proportion", title = "Comparison of Cell Proportions (FDR Adjusted)") +
  theme(legend.position = "none", axis.text.x = element_text(face = "bold"))


print(p_box)
ggsave("./sc/PC_vs_CC_cell_comparison_boxplot.pdf",
       plot = p_box, 
       width = 16,    
       height = 16,  
       limitsize = FALSE)


###sort by pc/cc proportion in different cell type----
cell_type_props=readRDS("F:\\Chordoma\\Plots\\cell_type_props.rds")

my_colors <- c("CC" = "#61bada", "PC" = "#72567a")
p_bar_sorted <- ggplot(cell_type_props, aes(y = cell_type, x = freq, fill = clinical_group)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7, color = "black", size = 0.2) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    title = "Cell Type Composition (Sorted by PC Proportion)",
    y = "Cell Type",      
    x = "Proportion",     
    fill = "Group"
  ) +
  
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"), 
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_bar_sorted)
ggsave("./sc/Barplot_CellType_Composition_in_CC_PC_sorted.pdf", 
       plot = p_bar_sorted , 
       width = 10,
       height = 16)


