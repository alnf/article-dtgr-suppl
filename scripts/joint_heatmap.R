library(RColorBrewer)

segment = 3
hm <- readSegmentData(segment)

hm.df <- read_excel("../data/220525 Protein list of interest with annotation.xlsx", col_names = T)

hname <- levels(factor(hm.df$heatmap_name))
hid = 1

hm.df <- hm.df[which(hm.df$heatmap_name==hname[hid]),]

hm.f <- hm[which(hm$id %in% hm.df$gene_name),]
hm.f <- hm.f[hm.f[,paste("P.Value_WT",segment,".over.Veh", segment, sep="")]<=0.01,]
hm.f <- hm.f[which(!duplicated(hm.f$id)),]

an.df <- hm.df[hm.df$gene_name %in% hm.f$id, ] 

complex_df <- as.data.frame(unique(an.df[, c("gene_name", "complex")]))
rownames(complex_df) <- complex_df$gene_name
complex_df <- complex_df[hm.f$id, "complex", drop = FALSE]
unique_complex <- unique(complex_df$complex)
complex_colors <- setNames(brewer.pal(length(unique_complex), "Set3"), unique_complex)


pathway_mat <- table(an.df$gene_name, an.df$pathway_name)
pathway_mat <- pathway_mat[hm.f$id, , drop = FALSE] # match heatmap gene order
pathway_mat[pathway_mat > 1] <- 1 
pathway_colors <- c("0" = "white", "1" = "black")


scaled_mat <- scale(t(hm.f[,expr_cols]))

complex_anno <- HeatmapAnnotation(
  complex = complex_df$complex,
  col = list(complex = complex_colors),
  annotation_name_side = "left"
)

col_list <- lapply(colnames(pathway_mat), function(x) pathway_colors)
names(col_list) <- colnames(pathway_mat)

pathway_anno <- HeatmapAnnotation(
  df = as.data.frame.matrix(pathway_mat),
  col = col_list,
  annotation_name_side = "left",
  show_legend = FALSE,
  #annotation_name_rot = 45,
  show_annotation_name = TRUE # display pathway names
)

cha_horizontal <- rowAnnotation(df = anno_df[,1,drop=F], col = list(group = gcols))
ht_horizontal <- Heatmap(scaled_mat, column_title = hname[hid], heatmap_legend_param = list(title = "expression"),
                         row_labels = anno_df$sample,
                         column_labels = hm.f$id,
                         left_annotation = cha_horizontal,
                         top_annotation = c(complex_anno, pathway_anno))

draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))

width = 0.3*(nrow(hm.f))
if (nrow(hm.f)<=11) {
  width = 0.31*(nrow(hm.f)+8)
}

png(paste0("../plots/selected_heatmaps/S", segment, "/heatmap S", segment, " ", hname[hid], ".png"),
    width = width, height = 5, units="in", res=100)
draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

 
