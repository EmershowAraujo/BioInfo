library(miRNAtap)
library(miRNAtap.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)


#X = miRNA//db = NÃO ALTERAR//spc = Especie
tmp <- g2Predicts(x = c("hsa-miR-155",
                        "Hsa-miR-103a",
                        "Hsa-miR-29a",
                        "Hsa-miR-30b",
                        "Hsa-miR-125",
                        "Hsa-miR-143",
                        "Hsa-miR-145"))

geneName <- AnnotationDbi::select(x = org.Hs.eg.db,
                                  keys = tmp$targetID,
                                  columns = "SYMBOL",
                                  keytype = ("ENTREZID"))

tmp %>% 
  mutate(symbol = geneName$SYMBOL) %>% 
  write_tsv(file = "predictEmershow.tsv")

tmp %>% 
  mutate(symbol = geneName$SYMBOL) %>% 
  group_by(miRNA) %>% 
  count()

#Genes interesse
tmp %>% 
  mutate(symbol = geneName$SYMBOL) %>% 
  filter(str_detect(string = symbol,pattern = "RUNX2|FOXO|AKT"))


miRNA_tags <- tmp %>% 
  pivot_wider(names_from = miRNA,
              values_from = targetID)
#KEGG
KEGG <- list()


for(i in 1:ncol(miRNA_tags)){
  KEGG[[i]] <- enrichKEGG(gene = miRNA_tags[[i]][[1]],
                          organism     = 'hsa',
                          #nPerm        = 1000,
                          minGSSize    = 3,
                          maxGSSize    = 800,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH", ##no p-value correction or BH
                          keyType       = "ncbi-geneid")
}
names(KEGG) <- colnames(miRNA_tags)

for(i in 1:length(KEGG)){
  nm <- names(KEGG)[i]
  KEGG[[i]]@result %>% 
    write_tsv(file = paste(nm,"KEGG",".tsv",sep = ""))
  
}

#GO
GO <- list()


for(i in 1:ncol(miRNA_tags)){
  GO[[i]] <- enrichGO(gene = miRNA_tags[[i]][[1]],
                      ont = "ALL",
                      OrgDb = org.Hs.eg.db,
                      #nPerm        = 500,
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      keyType       = "ENTREZID")
}
names(GO) <- colnames(miRNA_tags)

for(i in 1:length(GO)){
  nm <- names(GO)[i]
  GO[[i]]@result %>% 
    write_tsv(file = paste(nm,"GO",".tsv",sep = ""))
  
}
