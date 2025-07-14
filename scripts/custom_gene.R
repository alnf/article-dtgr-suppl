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

