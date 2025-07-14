library(RColorBrewer)
library(readxl)
library(dplyr)
library(fgsea)
library(ggplot2)
library(GSVA)
library(limma)
library(ComplexHeatmap)
library(circlize)

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

### GSEA
ranks <- setNames(-hm[,paste0("logFC_WT",segment,".over.Veh", segment)], hm$id)
dup <- which(duplicated(names(ranks)))
ranks <- ranks[-dup]
ranks <- sort(ranks, decreasing=TRUE)
pathways <- fgsea::gmtPathways("../pathways/m2.cp.reactome.v2024.1.Mm.symbols.gmt")

fg <- fgsea(pathways=pathways,
            stats   =ranks,
            nperm   =10000)

top <- fg[order(fg$padj), ][1:20, ]
top$pathway <- str_remove(top$pathway, "^[^_]*_")

ggplot(top, aes(x=NES, 
                y=reorder(pathway, NES), 
                size=size, 
                color=padj)) +
  geom_point() +
  scale_color_viridis_c(trans="log10", direction=-1) +
  labs(x="Normalized Enrichment Score (NES)",
       y=NULL,
       color="adj. p-value",
       size="Genes in set") +
  theme_minimal()

ggsave(paste("../plots/gsea/gsea_S",segment,".png", sep=""), width = 15, height = 6, scale = 0.8, dpi = 100)

### ssGSEA

hm <- hm[-which(duplicated(hm$id)),]
expr.mat <- hm[,expr_cols]
rownames(expr.mat) <- hm$id
pathways <- fgsea::gmtPathways("../pathways/m2.cp.reactome.v2024.1.Mm.symbols.gmt")

params <- ssgseaParam(
  exprData  = as.matrix(expr.mat),      # <-- was `expr` before
  geneSets  = pathways,      # <-- was `gset.idx.list` before
  assay     = NA,            # usually ignored for plain matrices
  minSize   = 5,             # filter out very small gene‐sets, if you like
  maxSize   = Inf,
  alpha     = 0.75,          # <-- was `tau`; exponent for the running‐sum
  normalize = F           # <-- was `ssgsea.norm`
)

ssgsea.scores <- gsva(params, verbose = FALSE)


anno_df$group <- factor(anno_df$group)
design <- model.matrix(~ 0 + group, data=anno_df)
colnames(design) <- levels(anno_df$group)  # “Veh” “WT”

# Fit limma
fit   <- lmFit(ssgsea.scores, design)
cont  <- makeContrasts(WTvsdTGR = WT - dTGR, levels=design)
fit2  <- contrasts.fit(fit, cont)
fit2  <- eBayes(fit2)

# Extract results
path.res <- topTable(fit2,
                     coef        = "WTvsdTGR",
                     number      = Inf,
                     sort.by     = "P",
                     adjust.method = "BH")


pathwaySizes <- lengths(pathways)
names(pathways)
path.res$setSize <- pathwaySizes[rownames(path.res)]

path.res <- path.res[path.res$adj.P.Val<=0.1,]
topN <- head(path.res[order(path.res$adj.P.Val), ], 30)
topN <- head(path.res[order(path.res$logFC), ], 10)
topN <- path.res
  
ggplot(topN, aes(
  x     = logFC,
  y     = reorder(rownames(topN), logFC),
  size  = setSize,
  color = adj.P.Val
)) +
  geom_point() +
  scale_color_viridis_c(trans="log10", direction=-1) +
  labs(
    x     = "ssGSEA enrichment (WT – dTGR)",
    y     = NULL,
    size  = "Genes in set",
    color = "adj. p-value"
  ) +
  theme_minimal(base_size=14)


# Plotting scores

vars <- apply(ssgsea.scores, 1, var)
top20 <- names(sort(vars, decreasing = TRUE))[1:20]
mat <- ssgsea.scores[top20, ]

# 3) Scale each pathway (row) to Z‐scores for better contrast
mat_z <- t(scale(t(mat)))

# 4) Column annotation: add group info
#    Make sure anno_df has rownames matching colnames(ssgsea.scores)
ha_col <- HeatmapAnnotation(
  group = anno_df[colnames(mat_z), "group"],
  col = list(group = gcols)
)

# 5) Define a diverging color function
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# 6) Draw the heatmap
Heatmap(
  mat_z,
  name = "ssGSEA Z-score",
  col = col_fun,
  top_annotation = ha_col,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_title = "Pathways",
  column_title = "Samples",
#  clustering_distance_rows    = "euclidean",
#  clustering_distance_columns = "euclidean",
#  clustering_method           = "complete",
  heatmap_legend_param = list(
    title = "Z-score",
    at    = c(-2, -1, 0, 1, 2),
    labels= c("-2","-1","0","1","2")
  )
)
