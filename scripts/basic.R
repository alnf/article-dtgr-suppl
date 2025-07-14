library(ComplexHeatmap)
library(ggplot2)
library(ggfortify)
library(EnhancedVolcano)
library(stringr)
source("utils.R")

### Read data
segment = 2
hm <- readSegmentData(segment)

# Select needed columns
expr_cols <- c(23:ncol(hm))
expr_cols <- colnames(hm[,expr_cols]) %>%
  str_subset("WT|Veh") %>%
  str_subset("treated", negate = TRUE)
expr_cols <- which(colnames(hm) %in% expr_cols)
# Prepare annotation table
anno_df <- data.frame(row.names = colnames(hm)[expr_cols], group = colnames(hm)[expr_cols], sample = colnames(hm)[expr_cols])
anno_df$group <- sapply(strsplit(anno_df$group,"_"), `[`, 1)
anno_df$sample <- sapply(strsplit(anno_df$sample,"_", fixed=T), `[`, 3)
anno_df$group[which(anno_df$group=="Veh")] <- "dTGR"

### PCA
pca <- prcomp(t(hm[,expr_cols]))

autoplot(pca, data = anno_df,
         colour = 'group',
         #label = TRUE, label.label = "sample", label.repel=T,
         frame=T, frame.type="norm",
         label.show.legend = F,
         frame.colour = 'group') +
  scale_color_manual(values = gcols) +
  scale_fill_manual(values = gcols) +
  scale_y_continuous(
    breaks  = seq(-2, 2, by = 1),
    limits = c(-2.2,2.2)
  ) +
  scale_x_continuous(
    breaks  = seq(-2, 2, by = 1),
    limits = c(-2.2,2.2)
  ) +  
  coord_fixed(ratio = 1)

ggsave(paste("../plots/pca/pca_S",segment,".png", sep=""), width = 5.5, height = 5, scale = 0.8, dpi = 100)
ggsave(paste("../plots/pca/pca_S",segment,".eps", sep=""), device=cairo_ps, width = 5.5, height = 5, scale = 0.8, dpi = 100)

### Volcano
hm[,paste0("logFC_WT",segment,".over.Veh", segment)] <- -hm[,paste0("logFC_WT",segment,".over.Veh", segment)]

EnhancedVolcano(
  hm,
  lab = hm$id,
  x = paste0("logFC_WT",segment,".over.Veh", segment),
  y = paste0("adj.P.Val_WT",segment,".over.Veh", segment),
  ylim = c(0, 3),
  xlim = c(-6,6),
  pCutoff = 0.05,
  FCcutoff = 0.5,
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
  legend.box.spacing = unit(0, "pt"),
  legend.position      = "top",
  legend.justification = "left",
  legend.direction     = "horizontal",
  legend.box.margin    = margin(t = 0, r = 0, b = 0, l = -35, unit = "pt")
)

ggsave(paste0("../plots/volcano/volcano_S",segment,"_adjp005_lfc05.png"), width = 7, height = 7.5, scale = 0.8, dpi = 100)
ggsave(paste0("../plots/volcano/volcano_S",segment,"_adjp005_lfc05.eps"), device=cairo_ps, width = 7, height = 7.5, scale = 0.8, dpi = 100)

### DEGs heatmap
# 0.05
hm.f <- hm[abs(hm[,paste("logFC_WT",segment,".over.Veh", segment, sep="")])>0.5,]
hm.f <- hm.f[hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]<=0.05,]

gcols=c("WT"="#808080", "Veh"="#ff8000", "treated"="#0000ff", "dTGR"="#ff8000")

scaled_mat <- scale(t(hm.f[,expr_cols]))
cha_horizontal <- rowAnnotation(df = anno_df[,1,drop=F], col = list(group = gcols))
ht_horizontal <- Heatmap(scaled_mat, 
                         name = "expression",
                         row_labels = anno_df$sample,
                         column_labels = hm.f$id,
                         left_annotation = cha_horizontal)
draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
width = 0.23*(nrow(hm.f))
if (nrow(hm.f)<=20) {
  width = 0.24*(nrow(hm.f)+4.5)
}
if (nrow(hm.f)<=6) {
 width = 0.23*(nrow(hm.f)+4.5)
}
if (nrow(hm.f)<=3) {
  width = 0.23*(nrow(hm.f)+9)
  print(width)
}


postscript(paste0("../plots/degs_heatmaps/heatmap_degs_S", segment,"_adjp005_lfc05.eps"), width = width, height = 3.3,
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special")
draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
dev.off()

png(paste0("../plots/degs_heatmaps/heatmap_degs_S", segment,"_adjp005_lfc05.png"), width = width, height = 3.3, units="in", res=100)
draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
dev.off()

### Heatmap for selected proteins
sp <- read_excel("../data/220525 Protein list of interest_metabolism.xlsx", col_names = F)
colnames(sp) <- c("id", "pathway")
pathways <- levels(factor(sp$pathway))

for (i in 1:length(pathways)) {
  hm.f <- hm[which(hm$id %in% sp[which(sp$pathway %in% pathways[i]),]$id),]
  hm.f <- hm.f[hm.f[,paste("P.Value_WT",segment,".over.Veh", segment, sep="")]<=0.01,]
  
  if (nrow(hm.f)>0) {
    scaled_mat <- scale(t(hm.f[,expr_cols]))
    cha_horizontal <- rowAnnotation(df = anno_df[,1,drop=F], col = list(group = gcols))
    ht_horizontal <- Heatmap(scaled_mat, column_title = pathways[i], heatmap_legend_param = list(title = "expression"),
                             row_labels = anno_df$sample,
                             column_labels = hm.f$id,
                             left_annotation = cha_horizontal)
    width = 0.23*(nrow(hm.f))
    if (nrow(hm.f)<=11) {
      width = 0.25*(nrow(hm.f)+8)
    }
    png(paste0("../plots/selected_heatmaps/S", segment, "/heatmap S", segment, " ", gsub("/", " ", pathways[i]),".png"),
        width = width, height = 3.3, units="in", res=100)
    draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
    dev.off()
    postscript(paste0("../plots/selected_heatmaps/S", segment, "/heatmap S", segment, " ", gsub("/", " ", pathways[i]),".eps"), width = width, height = 3.3,
               horizontal = FALSE, 
               onefile = FALSE, 
               paper = "special")
    draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
    dev.off()    
  }
}
