library(ComplexHeatmap)
library(readxl)
library(ggplot2)
library(ggfortify)
library(EnhancedVolcano)
library(stringr)

### Read data
segment = 3
hm <- read_excel(paste("Segment ", segment, "_Proteome Data_Für Alina.xlsx", sep=""))
hm <- data.frame(hm)
rownames(hm) <- hm$id
colnames(hm)

fc_cols <- c(paste("logFC_WT",segment,".over.Veh", segment, sep=""))
fc_cols_names <- c("WT vs Veh")

hm$id <- sapply(strsplit(hm$id,"_"), `[`, 1)
hm$id[which(hm$id=="")] = rownames(hm)[which(hm$id=="")]
hm$id = gsub("_", "", hm$id)

### Expression heatmap

expr_cols <- c(23:ncol(hm))
expr_cols <- colnames(hm[,expr_cols]) %>%
  str_subset("WT|Veh") %>%
  str_subset("treated", negate = TRUE)
expr_cols <- which(colnames(hm) %in% expr_cols)

anno_df <- data.frame(row.names = colnames(hm)[expr_cols], group = colnames(hm)[expr_cols], sample = colnames(hm)[expr_cols])
anno_df$group <- sapply(strsplit(anno_df$group,"_"), `[`, 1)
anno_df$sample <- sapply(strsplit(anno_df$sample,"_", fixed=T), `[`, 3)

# 0.05
hm.f <- hm[abs(hm[,paste("logFC_WT",segment,".over.Veh", segment, sep="")])>0.5,]
hm.f <- hm.f[hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]<=0.05,]

gcols=c("WT"="#808080", "Veh"="#ff8000", "treated"="#0000ff")

cha = HeatmapAnnotation(df = anno_df[,1,drop=F], col = list(group=gcols))
ht <- Heatmap(t(scale(t(hm.f[,expr_cols]))), name = "expression", column_labels = anno_df$sample, top_annotation = cha,
              row_labels = hm.f$id)

draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))

png(paste("sarah_new_heatmap_exprs_S", segment,"_005_05.png", sep=""), width = 6, height = 0.25*(nrow(hm.f)), units="in", res=100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
dev.off()

### PCA

pca <- prcomp(t(hm[,expr_cols]))

autoplot(pca, data = anno_df,
         colour = 'group',
         #label = TRUE, label.label = "sample", label.repel=T,
         frame=T, frame.type="norm",
         label.show.legend = F,
         frame.colour = 'group') +
  scale_color_manual(values = gcols) +
  scale_fill_manual(values = gcols)

ggsave(paste("sarah_new_pca_S",segment,".png", sep=""), width = 7, height = 6, scale = 0.8, dpi = 100)

### Volcano
EnhancedVolcano(
  hm,
  lab = hm$id,
  x = paste0("logFC_WT",segment,".over.Veh", segment),
  y = paste0("adj.P.Val_WT",segment,".over.Veh", segment),
  ylim = c(0, 3),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 5.0,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  legendLabSize = 12,
  legendIconSize = 4,
  legendPosition = "top"
) +
theme(
#  legend.position.inside    = c(0.01, 0.98),
  legend.box.spacing = unit(0, "pt"),
  legend.position      = "top",
  legend.justification = "left",
  legend.direction     = "horizontal",
  legend.box.margin    = margin(t = 0, r = 0, b = 0, l = -35, unit = "pt")
)

ggsave(paste("sarah_new_volcano_S",segment,".png", sep=""), width = 7, height = 7.5, scale = 0.8, dpi = 100)

### Heatmap for selected proteins

sp <- read_excel("220525 Protein list of interest_metabolism.xlsx", col_names = F)
colnames(sp) <- c("id", "pathway")
pathways <- levels(factor(sp$pathway))
i=1
hm.f <- hm[which(hm$id %in% sp[which(sp$pathway %in% pathways[i]),]$id),]

gcols=c("WT"="#808080", "Veh"="#ff8000", "treated"="#0000ff")
cha = HeatmapAnnotation(df = anno_df[,1,drop=F], col = list(group=gcols))
ht <- Heatmap(t(scale(t(hm.f[,expr_cols]))), column_title = pathways[i], heatmap_legend_param = list(title = "expression"),
              column_labels = anno_df$sample, top_annotation = cha,
              row_labels = hm.f$id)
ht

png(paste0("heatmap S", segment, " ", gsub("/", " ", pathways[i]),".png"), width = 6, height = 0.25*(nrow(hm.f)), units="in", res=100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
dev.off()
