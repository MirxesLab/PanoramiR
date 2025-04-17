# AIM ----------
# update the miRNA to gene relationship using miRTarBase v10
# Resources: https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/index.php 
  # MTI.csv

# ENVIRONMENT ------
library(tidyverse)
setwd('~/Documents/GitHub/PanoramiR/resources/')

file_out = "Predicted_Targets_Human_Filtered_Rearranged_GeneID.txt"
mirBase = read_csv('table/hsa_miRTarBase_MTI_v10.csv')
mirBase$`Support Type` %>% table() %>% as.data.frame()
# 1            Functional MTI   18312
# 2     Functional MTI (Weak) 3984724
# 3        Non-Functional MTI    4282
# 4 Non-Functional MTI (Weak)     430

mirBase_strong = mirBase %>% 
  dplyr::filter(`Support Type` %in% c('Functional MTI', 'Non-Functional MTI')) %>%
  dplyr::select(2, 4, 5) %>%
  unique()
colnames(mirBase_strong) = c('miRNA', 'Gene', 'Gene ID')


library(org.Hs.eg.db)
gene = mirBase_strong$Gene %>% unique()
ensembl = mapIds(org.Hs.eg.db, keys = gene, 
                 column = 'ENSEMBL', keytype = 'SYMBOL')
df.ensembl = data.frame(Gene = names(ensembl), 'Ensembl ID' = ensembl) %>% unique()
newmap = mirBase_strong %>% left_join(df.ensembl, by = 'Gene') # 9305
newmap = unique(newmap) 

# Save results
write.table(newmap, file = file.path("table", file_out),quote = FALSE, sep = '\t')


