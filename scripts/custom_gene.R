library(AnnotationDbi)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(homologene)

but_genes <- read.table("../data/butyrate_genes.tsv", sep="\t", header = T)

human_syms <- but_genes$symbol   # ← replace with your list

# 2) Map human SYMBOL → human Entrez
human_ent <- mapIds(
  org.Hs.eg.db,
  keys       = human_syms,
  column     = "ENTREZID",
  keytype    = "SYMBOL",
  multiVals  = "first"
)
# drop any that failed
human_ent <- human_ent[!is.na(human_ent)]

# 3) Use homologene() to get rat Entrez IDs
#    inTaxon = 9606 (Homo sapiens), outTaxon = 10116 (Rattus norvegicus)
homologs <- homologene::homologene(
  unique(human_ent),
  inTax   = 9606,
  outTax  = 10116
)


### Read data
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

hm.f <- hm[which(toupper(hm$id) %in% but_genes$symbol),]
scaled_mat <- scale(t(hm.f[,expr_cols]))
cha_horizontal <- rowAnnotation(df = anno_df[,1,drop=F], col = list(group = gcols))
ht_horizontal <- Heatmap(scaled_mat, column_title = "Butyrate genes", heatmap_legend_param = list(title = "expression"),
                         row_labels = anno_df$sample,
                         column_labels = hm.f$id,
                         left_annotation = cha_horizontal)

width = 0.23*(nrow(hm.f))
if (nrow(hm.f)<=11) {
  width = 0.25*(nrow(hm.f)+8)
}

draw(ht_horizontal, merge_legend = TRUE, padding = unit(c(2, 2, 2, 5), "mm"))


## Violin plot

genename = "BDH1"
genename = "OXCT1"
genename = "ACAT1"
hm.f <- hm[which(toupper(hm$id) == genename),]
hm.f[,expr_cols]
anno_df


expr_df <- hm.f[,expr_cols] %>%
  as.data.frame() %>% 
  rownames_to_column("feature") %>% 
#  filter(feature == "Bdh1_A0A0G2JSH2") %>%
  dplyr::select(-feature) %>% 
  pivot_longer(
    cols      = everything(),
    names_to  = "sample",
    values_to = "expression"
  )

anno_df$sample <- rownames(anno_df)
# 2. Join to your annotation to get group info
plot_df <- expr_df %>%
  left_join(anno_df, by = "sample")

# 3. Violin plot
gcols=c("WT"="#808080", "dTGR"="#ff8000")

p.adj <- hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]

ymax   <- max(plot_df$expression)
y.label <- ymax + 0.05 * (ymax - min(plot_df$expression))

ggplot(plot_df, aes(x = group, y = expression, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +   # overlay a narrow boxplot
  geom_jitter(width = 0.15, size = 1, alpha = 0.8) + # show individual points
  scale_fill_manual(values = gcols) +
  annotate(
    "text",
    x     = 1.5,                                           # between the two groups
    y     = y.label+0.5,
    label = paste0("adj. pval = ", signif(p.adj, 2)),
    size  = 5,
    fontface = "italic"
  ) +
  labs(
    title = paste0(genename, " expression by group"),
    x     = "Group",
    y     = "log2 expression"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

ggsave(paste("../plots/S",segment,"_",genename,".png", sep=""), width = 5.5, height = 5, scale = 0.8, dpi = 100)
