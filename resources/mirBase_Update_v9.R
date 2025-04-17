library(readxl)
library(tidyverse)
library(ggVennDiagram)

setwd('/Users/plateau/Documents/GitHub/PanoramiR/resources/')

# ----- read in files -----
mirBase = read_xlsx('table/hsa_MTI.xlsx')
oldmap  = read_tsv('table/old.txt')
colnames(oldmap)[1:2] = c('miRNA', 'Gene')
colnames(mirBase)[4] = 'Gene'

# ----- Check data -----
df = table(mirBase$`Support Type`) %>% as.data.frame()
tmp = data.frame(Var1 = "No Recording",
                 Freq = sum(is.na(mirBase$`Support Type`)))
df = rbind(tmp, df)
df$type = c('NA', 'Strong', 'Weak', 'Strong', 'Weak')
ggplot(df, aes(x = Var1, y = Freq, fill = type)) +
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))

mirBase.support = mirBase[!is.na(mirBase$`Support Type`), ]
mirBase.strong  = mirBase.support[!str_detect(mirBase.support$`Support Type`, 'Weak'), ]

# miRNA check
# 1. Number of miRNA
df.miRNA = data.frame(Map = factor(c('miRTarBase_all', 'miRTarBase_support', 'miRTarBase_strong', 'Old Map'),
                                   levels = c('miRTarBase_all', 'miRTarBase_support', 'miRTarBase_strong', 'Old Map')),
                      miRNA = c(length(unique(mirBase$miRNA)),
                                length(unique(mirBase.support$miRNA)),
                                length(unique(mirBase.strong$miRNA)),
                                length(unique(oldmap$miRNA))))

ggplot(df.miRNA, aes(x = Map, y = miRNA, fill = Map)) +
    geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))
# 2. Common miRNA and Unique miRNA 
ggVennDiagram(list(Old = unique(oldmap$miRNA),
                               miRTarBase = unique(mirBase.strong$miRNA))) +
    scale_fill_distiller(palette = 'YlGnBu')

# 3. Relationship between miRNA and Gene
idx = oldmap$miRNA %in% mirBase.strong$miRNA
common.old = oldmap[idx, ]
unique.old   = oldmap[!idx, ]
idx = mirBase.strong$miRNA %in% oldmap$miRNA
common.new   = mirBase.strong[idx, ]
unique.new   = mirBase.strong[!idx, ]

fun.gene = function(map, type) {
    tmp = group_split(map, by = miRNA)
    num = lapply(tmp, function(x) nrow(x)) %>% unlist()
    num.name = lapply(tmp, function(x) x$miRNA[1]) %>% unlist()
    df = data.frame(miRNA = num.name,
                    number_gene = num,
                    type = type)

    return(df)
}

old.common.gene = fun.gene(common.old, 'old_common')
new.common.gene = fun.gene(common.new, 'new_common')
old.unique.gene = fun.gene(unique.old, 'old_unique')
new.unique.gene = fun.gene(unique.new, 'new_unique')

df = rbind(old.common.gene,
           new.common.gene,
           old.unique.gene,
           new.unique.gene) 
ggplot(df, aes(x = type, y = number_gene)) +
    geom_boxplot() +
    theme_classic()


df.gene.summary = tibble(summary(old.common.gene$number_gene),
       summary(old.unique.gene$number_gene),
       summary(new.common.gene$number_gene),
       summary(new.unique.gene$number_gene))
colnames(df.gene.summary)= c('old_common', 'old_unique',
                             'new_common', 'new_unique')
df.gene.summary = df.gene.summary %>% as.data.frame()
rownames(df.gene.summary) = c('Min', 
                              '1st Qu',
                              'Median',
                              'Mean',
                              '3rd Qu',
                              'Max')
df.gene.summary = round(df.gene.summary, 2)

common.miRNA = common.old$miRNA %>% unique()
fun.commongene = function(id) {
    mirna = common.miRNA[id]
    tmp.old = common.old[common.old$miRNA == mirna, ]
    tmp.new = common.new[common.new$miRNA == mirna, ]
    num.old = nrow(tmp.old)
    num.new = nrow(tmp.new)
    num.common = sum(tmp.new$Gene %in% tmp.old$Gene)
    num.unique.new = num.new - num.common
    num.unique.old = num.old - num.common
    
    df = data.frame(miRNA = mirna,
                    Old_all = num.old,
                    Old_unique = num.unique.old,
                    New_all = num.new,
                    New_unique = num.unique.new,
                    Common  = num.common)
    return(df)
}

for (i in 1:length(common.miRNA)) {
    if (i == 1) {
        df.commongene = fun.commongene(i)
    } else {
        tmp = fun.commongene(i)
        df.commongene = rbind(df.commongene, tmp)
    }
}

df.commongene$Perc_New_common = round(df.commongene$Common / df.commongene$New_all, 2)
ggplot(df.commongene, aes(x = Perc_New_common)) +
    geom_density()



# ----- Save new map -----
library(org.Hs.eg.db)
newmap = mirBase.strong[, c(2, 4, 5)] %>% unique() 
colnames(newmap) = colnames(oldmap)[-4]
gene = newmap$Gene
ensembl = mapIds(org.Hs.eg.db, keys = gene, 
                 column = 'ENSEMBL', keytype = 'SYMBOL')
df.ensembl = data.frame(Gene = names(ensembl),
                        'Ensembl ID' = ensembl) %>% 
    unique()
newmap = left_join(newmap, df.ensembl, by = 'Gene') # 90393
newmap = unique(newmap) # 8489

# newmap$miRNA %>% unique () %>% length() # 740
# newmap$Gene %>% unique() %>% length() # 2849
# paste(newmap$miRNA, newmap$Gene) %>% unique() %>% length() # 8489
# 
# oldmap$Gene %>% unique() %>% length()
# paste(oldmap$miRNA, oldmap$Gene) %>% unique() %>% length()
# unique(oldmap) %>% nrow()

# Final Check
miRNA_Gene = fun.gene(newmap, 'Final_Map')
mean(miRNA_Gene$number_gene)
median(miRNA_Gene$number_gene)
bar = ggplot(miRNA_Gene, aes(x = number_gene)) +
    geom_histogram(bins = 50)
violin = ggplot(miRNA_Gene, aes(x = type, y = number_gene)) +
    geom_violin() +
    coord_flip() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = 1))
ggarrange(bar, violin,  nrow = 2)
summary(miRNA_Gene$number_gene)
# Save results
write.table(newmap, file = 'table/Predicted_Targets_Human_Filtered_Rearranged_GeneID.txt',quote = FALSE, sep = '\t')
