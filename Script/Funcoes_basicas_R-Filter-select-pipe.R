#importa e ñ salva
read_tsv(file = "Data/Hsa-miR-103aGO.tsv" )
#import e salva na variavel(objeto)
miR103aGO <- read_tsv(file = "Data/Hsa-miR-103aGO.tsv")
#filtro column
miR103aGO$Description
#filtro linha/column
miR103aGO[1:10,"GeneRatio"]
 library(tidyverse)
#Pipe/filtro/seleçao
miR103aGO %>%
  select(pvalue, ID, Description) %>%
  filter(pvalue<0.05)
?select

