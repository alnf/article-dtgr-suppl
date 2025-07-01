library(AnnotationDbi)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(homologene)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(openxlsx)
source("utils.R")

### Read data
segment = 6
hm <- readSegmentData(segment)

### Get ids mappings 
rat_uniprots <- sapply(strsplit(rownames(hm),"_", fixed=T), `[`, 2)
rat_entrez <- mapIds(
  x         = org.Rn.eg.db,
  keys      = rat_uniprots,
  column    = "ENTREZID",
  keytype   = "UNIPROT",
  multiVals = "first"
)
rat_entrez <- data.frame(rat_uniprot=names(rat_entrez), rat_entrez=rat_entrez, row.names=NULL)

homo_df <- homologene(rat_entrez$rat_entrez, inTax=10116, outTax=9606)
colnames(homo_df) <- c("rat_symbol", "human_symbol", "rat_entrez", "human_entrez")
mappings <- merge(rat_entrez, homo_df, by="rat_entrez")

### Perform enrichment
#gmt_file <- "../pathways/SysMyo_Muscle_Gene_Sets.txt"
files <- c("../pathways/GO_Cellular_Component_2025.txt", "../pathways/Reactome_Pathways_2024.txt",
           "../pathways/MGI_Mammalian_Phenotype_Level_4_2024.txt", "../pathways/GO_Biological_Process_2025.txt")

for (gmt_file in files) {

  gmt_df <- read.gmt(gmt_file)
  
  # 0.01
  hm.f <- hm[abs(hm[,paste("logFC_WT",segment,".over.Veh", segment, sep="")])>0.5,]
  hm.f <- hm.f[hm.f[,paste("P.Value_WT",segment,".over.Veh", segment, sep="")]<=0.01,]
  #hm.f <- hm.f[hm.f[,paste("adj.P.Val_WT",segment,".over.Veh", segment, sep="")]<=0.05,]
  gene_list <- sapply(strsplit(rownames(hm.f),"_", fixed=T), `[`, 2)
  gene_list <- mappings[mappings$rat_uniprot %in% gene_list,]
  
  ora_res <- enricher(
    gene         = gene_list$human_symbol,
    TERM2GENE    = gmt_df,
    pvalueCutoff = 0.05,
    pAdjustMethod= "BH"
  )
  
  if (nrow(ora_res) > 0) {
    print(nrow(ora_res))
    head(ora_res, n=20)
    
    pname = sapply(strsplit(gmt_file,"/"), `[`,3)
    pname = sapply(strsplit(pname,".txt"), `[`, 1)
    
    # Dot plot
    dotplot(ora_res, showCategory=10, title=gsub("_", " ", pname))
    
    ggsave(paste0("../plots/ora/S",segment,"/ora_dot_S",segment,"_",pname,"_p001_lfc05.png"), width = 8, height = 8, scale = 0.8, dpi = 100)
    ggsave(paste0("../plots/ora/S",segment,"/ora_dot_S",segment,"_",pname,"_p001_lfc05.eps"), device=cairo_ps, width = 8, height = 8, scale = 0.8, dpi = 100)
    
    enrich_df <- as.data.frame(ora_res)
    write.xlsx(enrich_df, file = paste0("../plots/ora/S",segment,"/ora_S",segment,"_",pname,"_p001_lfc05.xlsx"), rowNames = FALSE)
    
    # Heatmap
    hm.f$rat_uniprot <- sapply(strsplit(rownames(hm.f),"_", fixed=T), `[`, 2)
    glist <- hm.f[hm.f$rat_uniprot %in% gene_list$rat_uniprot,]
    glist <- glist[,c("rat_uniprot",paste("logFC_WT",segment,".over.Veh", segment, sep=""))]
    glist <- merge(glist, gene_list, by="rat_uniprot")
    foldChange <- glist[,2]
    names(foldChange) <- glist$human_symbol
    
    heatplot(ora_res, foldChange=-foldChange, showCategory=10) +
      scale_fill_gradient2(
        low      = "#327eba",
        mid      = "white",
        high     = "#e06663",
        midpoint = 0
      ) + ggtitle(gsub("_", " ", pname))
    
    width = 0.07*(length(foldChange))
    if (length(foldChange)<=6) {
      width = 0.07*(length(foldChange)+4.5)
    }
    ggsave(paste0("../plots/ora/S",segment,"/ora_heat_S",segment,"_",pname,"_p001_lfc05.png"), width = width, height = 6, scale = 0.8, dpi = 100)
    ggsave(paste0("../plots/ora/S",segment,"/ora_heat_S",segment,"_",pname,"_p001_lfc05.eps"), device=cairo_ps, width = width, height = 6, scale = 0.8, dpi = 100)
  }
}