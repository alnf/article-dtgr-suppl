library(readxl)
readSegmentData <- function(segment) {
  hm <- read_excel(paste0("../data/Segment ", segment, "_Proteome Data_Für Alina.xlsx"))
  hm <- data.frame(hm)
  rownames(hm) <- hm$id
  hm$id <- sapply(strsplit(hm$id,"_"), `[`, 1)
  hm$id[which(hm$id=="")] = rownames(hm)[which(hm$id=="")]
  hm$id = gsub("_", "", hm$id)
  return(hm)
}