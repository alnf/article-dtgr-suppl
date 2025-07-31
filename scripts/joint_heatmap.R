library(RColorBrewer)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)

segment = 3
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

hm.df <- read_excel("../data/220525 Protein list of interest with annotation.xlsx", col_names = T)
hname <- levels(factor(hm.df$heatmap_name))
levels(factor(hm.df$complex))

# Complexes plot
hm.f <- hm[which(hm$id %in% hm.df$gene_name),]
hm.f <- hm.f[hm.f[,paste("P.Value_WT",segment,".over.Veh", segment, sep="")]<=0.01,]
hm.f <- hm.f[which(!duplicated(hm.f$id)),]

hid = 2
an.df <- hm.df[which(hm.df$heatmap_name==hname[hid]),]
an.df <- an.df[an.df$gene_name %in% hm.f$id, ] 
an.df <- an.df[which(!duplicated(an.df$gene_name)),] 

complex_df <- unique(an.df[, c("gene_name", "complex")])
complex_df <- as.data.frame(complex_df)
rownames(complex_df) <- complex_df$gene_name

# Filter expression matrix
hm.f <- hm.f[hm.f$id %in% rownames(complex_df),]

# Sort genes by complex
corder <- c("I", "II", "Coenzym Q", "III", "Fe-S cluster (I, II, III)", "I, III, IV, V", "IV", "V","ETC associated", "electron transfer flavoprotein")
# Create a factor with levels in your desired order
complex_factor <- factor(complex_df$complex, levels = corder, ordered = TRUE)
# Sort complex_df by this factor
sorted_indices <- order(complex_factor, na.last = TRUE)
# Get sorted dataframe
complex_df <- complex_df[sorted_indices, ]

rownames(hm.f) <- hm.f$id
hm.f <- hm.f[rownames(complex_df), ]
expr_mat   <- hm.f[, expr_cols]              # genes in the right order
scaled_mat <- scale(t(expr_mat))
  
# Row annotation
cha_horizontal <- rowAnnotation(df = anno_df[,1,drop=F], col = list(group = gcols))

splits <- factor(complex_df$complex, levels = unique(complex_df$complex))

ht <- Heatmap(
  scaled_mat,
  name            = "expression",       # legend title + internal heatmap ID
  cluster_columns = FALSE,
  show_column_dend= FALSE,
  
  column_split    = splits,
  #column_gap      = unit(4, "mm"),     # give room for titles
  column_title    = levels(splits),    # <-- vector of split‐titles
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  #column_title_rot= 90,                # rotate them if you like
  
  row_labels      = anno_df$sample,
  column_labels   = hm.f$id,
  left_annotation = cha_horizontal
)

draw(ht, merge_legend = TRUE, padding = unit(c(8,9,2,5), "mm"))

# Saving the heatmap
width = 0.3 * (nrow(hm.f))
if (nrow(hm.f) <= 11) {
  width = 0.31 * (nrow(hm.f) + 8)
}

png(paste0("../plots/joint_heatmaps/S", segment, "/heatmap S", segment, " ", hname[hid], " sorted by complex.png"),
    width = width, height = 3.1, units = "in", res = 100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

postscript(paste0("../plots/joint_heatmaps/S", segment, "/heatmap S", segment, " ", hname[hid], " sorted by complex.eps"),
    width = width, height = 3.1,
    horizontal = FALSE, 
    onefile = FALSE, 
    paper = "special")
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

# Pathways plot
hm.f <- hm[which(hm$id %in% hm.df$gene_name),]
hm.f <- hm.f[hm.f[,paste("P.Value_WT",segment,".over.Veh", segment, sep="")]<=0.01,]
hm.f <- hm.f[which(!duplicated(hm.f$id)),]

hid = 1
an.df <- hm.df[which(hm.df$heatmap_name==hname[hid]),]
an.df <- an.df[an.df$gene_name %in% hm.f$id, ] 

# Filter expression matrix
hm.f <- hm.f[hm.f$id %in% an.df$gene_name,]

pathway_mat <- table(an.df$gene_name, an.df$pathway_name)
pathway_mat <- pathway_mat[hm.f$id, , drop = FALSE] # match heatmap gene order
pathway_mat[pathway_mat > 1] <- 1 
pathway_colors <- c("0" = "white", "1" = "black")

scaled_mat <- scale(t(hm.f[,expr_cols]))

col_list <- lapply(colnames(pathway_mat), function(x) pathway_colors)
names(col_list) <- colnames(pathway_mat)

# Cluster genes based on pathway membership to get ordering
gene_dist <- dist(pathway_mat, method = "binary")  # distance based on binary membership
gene_clust <- hclust(gene_dist)
gene_order <- gene_clust$order

# Reorder hm.f and pathway_mat by clustered gene order
hm.f <- hm.f[gene_order, ]
pathway_mat <- pathway_mat[gene_order, ]
#pathway_mat <- rbind(pathway_mat[-7, , drop = FALSE],
#                     pathway_mat[ 7, , drop = FALSE])

# Now create pathway_anno and heatmap as before
pathway_anno <- HeatmapAnnotation(
  df = as.data.frame.matrix(pathway_mat),
  col = col_list,
  annotation_name_side = "left",
  show_legend = FALSE,
  show_annotation_name = TRUE
)

scaled_mat <- scale(t(hm.f[,expr_cols]))

ht <- Heatmap(scaled_mat, column_title = hname[hid], heatmap_legend_param = list(title = "expression"),
                         row_labels = anno_df$sample,
                         column_labels = hm.f$id,
                         left_annotation = cha_horizontal,
                         top_annotation = pathway_anno,
                         cluster_columns = F)

draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))

# Saving the heatmap
width = 0.37 * (nrow(hm.f))
if (nrow(hm.f) <= 11) {
  width = 0.38 * (nrow(hm.f) + 8)
}

png(paste0("../plots/joint_heatmaps/S", segment, "/heatmap S", segment, " ", hname[hid], " sorted by pathway.png"),
    width = width, height = 4.1, units = "in", res = 100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

postscript(paste0("../plots/joint_heatmaps/S", segment, "/heatmap S", segment, " ", hname[hid], " sorted by pathway.eps"),
    width = width, height = 4.1,
    horizontal = FALSE, 
    onefile = FALSE, 
    paper = "special")
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()
