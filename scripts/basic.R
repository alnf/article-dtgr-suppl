library(ComplexHeatmap)
library(readxl)
library(broom)
library(ggplot2)
library(ggfortify)
library(EnhancedVolcano)

segment = 1
hm <- read_excel(paste("Segment ", segment, "_Proteome Data_Für Alina.xlsx", sep=""))
hm <- data.frame(hm)
rownames(hm) <- hm$id
colnames(hm)

fc_cols <- c(paste("logFC_WT",segment,".over.Veh", segment, sep=""),
             paste("logFC_WT",segment,".over.treated", segment, sep=""),
             paste("logFC_Veh",segment,".over.treated", segment, sep=""))

fc_cols_names <- c("WT vs Veh", "WT vs treated", "Veh vs treated")
hm$id <- sapply(strsplit(hm$id,"_"), `[`, 1)
hm$id[which(hm$id=="")] = rownames(hm)[which(hm$id=="")]
hm$id = gsub("_", "", hm$id)

### Fold changes heatmap

# cha = HeatmapAnnotation(comparisons = anno_block(labels = fc_cols_names[c(3, 1, 2)], gp = gpar(fill = 2:4)))
# ht <- Heatmap(hm[,fc_cols], name = "logFC", column_labels = fc_cols_names, top_annotation = cha, column_split=3,
#         show_column_names = F, row_labels = hm$id)
# draw(ht, padding = unit(c(2, 2, 2, 5), "mm"))
# 
# png("sarah_heatmap_logfc.png", width = 7, height = 8, units="in", res=100)
# draw(ht, padding = unit(c(2, 2, 2, 5), "mm"))
# dev.off()

### Expression heatmap

expr_cols <- c(23:ncol(hm))
anno_df <- data.frame(row.names = colnames(hm)[expr_cols], group = colnames(hm)[expr_cols], sample = colnames(hm)[expr_cols])
anno_df$group <- sapply(strsplit(anno_df$group,"_"), `[`, 1)
anno_df$sample <- sapply(strsplit(anno_df$sample,".", fixed=T), `[`, 2)

#hm <- hm[hm$P.Value_WT3.over.Veh3<=0.01,]
#hm <- hm[hm$P.Value_Veh3.over.treated3<=0.001,]
#hm <- hm[hm$adj.P.Val_Veh3.over.treated3<=0.05,]
#hm <- hm[hm$adj.P.Val_WT3.over.treated3<=0.05,]

# 0.05
hm.f <- hm[abs(hm[,paste("logFC_WT",segment,".over.Veh", segment, sep="")])>0.5,]
hm.f <- hm.f[hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]<=0.05,]

gcols=c("WT"="#808080", "Veh"="#ff8000", "treated"="#0000ff")

cha = HeatmapAnnotation(df = anno_df[,1,drop=F], col = list(group=gcols))
ht <- Heatmap(t(scale(t(hm.f[,expr_cols]))), name = "expression", column_labels = anno_df$sample, top_annotation = cha,
              row_labels = hm.f$id)

draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))

png(paste("sarah_heatmap_exprs_S", segment,"_005.png", sep=""), width = 8, height = 0.25*nrow(hm.f), units="in", res=100)
draw(ht, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))
dev.off()

# 0.1
hm.f <- hm[abs(hm[,paste("logFC_WT",segment,".over.Veh", segment, sep="")])>0.5,]
hm.f <- hm.f[hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]<=0.1,]
ht <- Heatmap(t(scale(t(hm.f[,expr_cols]))), name = "expression", column_labels = anno_df$sample, top_annotation = cha,
              row_labels = hm.f$id)
png(paste("sarah_heatmap_exprs_S", segment,"_01.png", sep=""), width = 8, height = 0.25*nrow(hm.f), units="in", res=100)
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

ggsave(paste("sarah_pca_S",segment,".png", sep=""), width = 9, height = 7, scale = 0.8, dpi = 100)


## Correlation

funcdata <- read_excel("GLS für Alina.xlsx")
funcdata <- as.data.frame(funcdata)
exprdata <- hm[,expr_cols]
colnames(exprdata) <- sapply(strsplit(colnames(exprdata),"_"), `[`, 3)
colnames(exprdata) <- sapply(strsplit(colnames(exprdata),".", fixed = T), `[`, 1)
anno_df$sampleID <- sapply(strsplit(anno_df$sample,"_"), `[`, 2)
rownames(funcdata) <- funcdata$ID
funcdata <- funcdata[colnames(exprdata),-1]

### Tidyverse solution

df_expr <- exprdata %>%
  as.data.frame() %>%
  rownames_to_column("protein") %>%
  pivot_longer(
    cols = -protein,
    names_to = "sample",
    values_to = "expression"
  )

df_func <- funcdata %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    cols = -sample,
    names_to = "parameter",
    values_to = "value"
  )

df_joined <- df_expr %>%
  inner_join(df_func, by = "sample")

df_corr <- df_joined %>%
  group_by(protein, parameter) %>%
  summarize(
    cor_est = cor(expression, value, method = "spearman",  use = "pairwise.complete.obs"),
    p_val   = cor.test(expression, value, method = "spearman")$p.value,
    .groups = "drop"
  )

df_corr %>%
  arrange(desc(cor_est)) %>%
  slice_head(n = 10)

df_corr %>%
  arrange(cor_est) %>%
  slice_head(n = 10)

write.table(df_corr, "df_corr_for_Sarah.tsv", sep="\t", row.names = F)


### Linear regression

df_regress_results <- df_joined %>%
  group_by(protein, parameter) %>%
  do(
    # Fit a linear regression model for each group:
    mod = lm(value ~ expression, data = .)
  ) %>%
  # Extract a tidy summary of model coefficients:
  mutate(
    tidied = list(tidy(mod)),   # slope, intercept, p-value, ...
    glanced = list(glance(mod)) # r.squared, adj.r.squared, etc.
  )

df_regress_results

df_coefficients <- df_regress_results %>%
  select(protein, parameter, tidied) %>%
  unnest(cols = c(tidied))

df_modelstats <- df_regress_results %>%
  select(protein, parameter, glanced) %>%
  unnest(cols = c(glanced))

df_slope <- df_coefficients %>%
  filter(term == "expression") %>%
  select(protein, parameter, estimate, p.value)

df_filtered <- df_joined %>%
  inner_join(df_modelstats %>% filter(r.squared > 0.7) %>% select(protein, parameter),
             by = c("protein", "parameter"))

df_annotation <- df_modelstats %>%
  mutate(label = paste0("R² = ", round(r.squared, 2))) %>% filter(r.squared > 0.7)

df_filtered <- df_filtered %>%
  left_join(anno_df[,c("group", "sampleID")], by = c("sample" = "sampleID"))

p1 <- ggplot(df_filtered %>% filter(parameter == "GLS"), aes(x = expression, y = value)) +
  geom_point(aes(color=group)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ protein, scales = "free") +
  # Add the R² label to each facet
  geom_text(
    data = df_annotation %>% filter(parameter == "GLS"),
    aes(label = label),
    x = Inf,         # Place at right edge
    y = Inf,         # Place at top edge
    hjust = 1.1,     # Adjust horizontal position slightly
    vjust = 1.5,     # Adjust vertical position slightly
    inherit.aes = FALSE
  ) +
  ggtitle("GLS") +
  theme(legend.position="none", plot.title = element_text(size=13, face="bold"))

p2 <- ggplot(df_filtered %>% filter(parameter == "Seg3"), aes(x = expression, y = value)) +
  geom_point(aes(color=group)) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ protein, scales = "free") +
  # Add the R² label to each facet
  geom_text(
    data = df_annotation %>% filter(parameter == "Seg3"),
    aes(label = label),
    x = Inf,         # Place at right edge
    y = Inf,         # Place at top edge
    hjust = 1.1,     # Adjust horizontal position slightly
    vjust = 1.5,     # Adjust vertical position slightly
    inherit.aes = FALSE
  ) +
  ggtitle("Seg3") +
  theme(plot.title = element_text(size=13, face="bold"))

wp <- wrap_plots(list(p1, p2), ncol=2)
wp

ggsave("proteins_correlations.png", width = 21, height = 7, scale = 0.8, dpi = 100)



####### New request
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

expr_cols <- c(23:31)
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
