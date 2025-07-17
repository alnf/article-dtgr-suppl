library(RColorBrewer)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)

### Reading files
mets <- read.table("../data/Heatmap_Seg3_short_average.csv", sep=",", header = T)
mets <- mets[c(1:116),-1]
mets_names <- read.table("../data/MTX_Pathway_Tag.csv", sep=";", header = T)
mets_names <- mets_names[,c(1:2)]

### Plotting with hierarchical clustering
mat_fc <- mets %>%
  dplyr::select(short_name_prod, Treatment, averagezscore) %>%
  pivot_wider(names_from = Treatment, values_from = averagezscore) %>%
  column_to_rownames("short_name_prod") %>%
  as.matrix()

head(mat_fc)

col_meta <- data.frame(
  Treatment = colnames(mat_fc),
  stringsAsFactors = FALSE
)
rownames(col_meta) <- colnames(mat_fc)

treatment_cols <- c(
  "WT vehicle"   = "#808080",
  "dTGR vehicle" = "#ff8000"
)

ha <- HeatmapAnnotation(
  df = col_meta,
  col = list(Treatment = treatment_cols),
  annotation_name_side = "left",
  show_annotation_name = FALSE
)

binary_df <- mets_names %>% 
  mutate(Present = 1L) %>% 
  pivot_wider(
    names_from   = Pathway,
    values_from  = Present,
    values_fill  = list(Present = 0L)
  ) %>% 
  column_to_rownames("short_name_prod")

missing_mets <- setdiff(rownames(mat_fc), rownames(binary_df))

if (length(missing_mets) > 0) {
  fill_block <- matrix(
    0L,
    nrow   = length(missing_mets),
    ncol   = ncol(binary_df),
    dimnames = list(missing_mets, colnames(binary_df))
  )
  binary_df <- rbind(binary_df, fill_block)
}

binary_df <- binary_df[rownames(mat_fc), , drop = FALSE]

ann_cols <- lapply(colnames(binary_logical), function(.col) {
  c("0" = "white", "1" = "black")
})
names(ann_cols) <- colnames(binary_logical)

ra_all <- rowAnnotation(
  df  = as.data.frame.matrix(binary_df), 
  col = ann_cols,
  annotation_name_side = "top",
  show_legend = FALSE,
  show_annotation_name = TRUE
)

ht <- Heatmap(
  mat_fc,
  name               = "avg z-scores",
  top_annotation     = ha,
  left_annotation    = ra_all,
  show_column_names  = FALSE
)

segment = 3
png(paste0("../plots/joint_heatmaps/S", segment, "/metabolites heatmap S", segment, ".png"),
    width = 7, height = 10, units = "in", res = 100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

### Ordered by pathways and z-scores together

pathway_df <- mets_names %>%
  column_to_rownames("short_name_prod")
row_pathway <- pathway_df[ rownames(mat_fc), "Pathway" ]
row_pathway[is.na(row_pathway)] <- "None"
row_pathway <- factor(row_pathway,
                      levels = c("Glycolysis", "PPP", "TCA", "Lipid synthesis", 
                                 "Ketone bodies", "Polyol", "AA","None"))
row_pathway

ordings <- lapply(levels(row_pathway), function(pw) {
  idx <- which(row_pathway == pw)
  idx[order(mat_fc[idx,1], decreasing=TRUE)]
})
names(ordings) <- levels(row_pathway)
global_ord <- unlist(ordings)

ht <- Heatmap(
  mat_fc,
  name         = "avg z-scores",
  top_annotation= ha,
  left_annotation = ra_all,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  row_title    = NULL,
  row_split    = row_pathway,
  row_order    = global_ord
)

draw(ht, merge_legend = TRUE)

segment = 3
png(paste0("../plots/joint_heatmaps/S", segment, "/metabolites heatmap S", segment, " sorted.png"),
    width = 6.5, height = 10, units = "in", res = 100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 9, 2, 5), "mm"))
dev.off()

### Transposed

mat_fc_t <- t(mat_fc)
ra_treat <- rowAnnotation(
  df = data.frame(Treatment = rownames(mat_fc_t)),
  col = list(Treatment = treatment_cols),
  show_annotation_name = FALSE
)

ann_cols <- lapply(colnames(binary_df), function(.col) {
  c("0" = "white", "1" = "black")
})
names(ann_cols) <- colnames(binary_df)

ca_met <- columnAnnotation(
  df = as.data.frame(binary_df),
  col = ann_cols,
  show_legend         = FALSE,
  show_annotation_name= TRUE
)

col_split  <- row_pathway
col_order  <- global_ord

ht_horiz <- Heatmap(
  mat_fc_t,
  name               = "avg z-scores",
  top_annotation     = ca_met,
  right_annotation    = ra_treat,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  column_split       = col_split,
  column_order       = col_order,
  column_title       = NULL
)

draw(ht_horiz, merge_legend = TRUE, padding = unit(c(2,9,2,5), "mm"))

segment <- 3
png(paste0("../plots/joint_heatmaps/S", segment, 
           "/metabolites heatmap S", segment, " sorted horiz.png"),
    width = 15, height = 5, units = "in", res = 100)
draw(ht_horiz, merge_legend = TRUE, padding = unit(c(2,9,2,5), "mm"))
dev.off()


### Exact order
mat_fc <- mets %>%
  dplyr::select(short_name_prod, Treatment, averagezscore) %>%
  pivot_wider(names_from  = Treatment,
              values_from = averagezscore) %>%
  column_to_rownames("short_name_prod") %>%
  .[ mets_names$short_name_prod, , drop = FALSE ] %>%
  as.matrix()

binary_df <- mets_names %>%
  mutate(Present = 1L) %>%
  pivot_wider(names_from   = Pathway,
              values_from  = Present,
              values_fill  = list(Present = 0L)) %>%
  column_to_rownames("short_name_prod") %>%
  .[ rownames(mat_fc), , drop = FALSE ]

mat_fc_t <- t(mat_fc)

treatment_cols <- c("WT vehicle"   = "#808080",
                    "dTGR vehicle" = "#ff8000")
ra_treat <- rowAnnotation(
  Treatment = rownames(mat_fc_t),
  col       = list(Treatment = treatment_cols),
  show_annotation_name = FALSE
)

ann_cols <- binary_df %>%
  names() %>%
  set_names() %>%
  map(~ c("0" = "white", "1" = "black"))

ca_met <- columnAnnotation(
  df                   = as.data.frame(binary_df),
  col                  = ann_cols,
  show_legend          = FALSE,
  show_annotation_name = TRUE
)

ht_horiz <- Heatmap(
  mat_fc_t,
  name               = "avg z-scores",
  top_annotation     = ca_met,
  right_annotation   = ra_treat,
  show_row_names     = TRUE,
  show_column_names  = TRUE,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE
)

draw(ht_horiz, merge_legend = TRUE, 
     padding = unit(c(2, 9, 2, 5), "mm"))

png(paste0("../plots/joint_heatmaps/S", segment, 
           "/metabolites heatmap S", segment, " exact order horiz.png"),
    width = 15, height = 5, units = "in", res = 100)
draw(ht_horiz, merge_legend = TRUE,
     padding = unit(c(2,9,2,5), "mm"))
dev.off()
