# =========================================================================== #
# Step 4.1 Overview Information
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.GlobalNorm
    # ls.compare.group

# output:
    # fun.pca(comp, cols)
    # fun.heatmap(df, main, anno.df, anno.col)
    
    # df.anno.all
    # color.anno.all
    # allmiRNA [for overview info heat map]

    # sample.dist
    # sample.dist.matrix
# =========================================================================== #

# Prepare Global Normalization Data
# df.input.data.GlobalNorm.tmp [Samples, miRNA ...]
df.input.data.GlobalNorm.tmp = t(df.input.data.GlobalNorm[, -1])
colnames(df.input.data.GlobalNorm.tmp) = df.input.data.GlobalNorm$miRNA

df.input.data.GlobalNorm.tmp = as.data.frame(df.input.data.GlobalNorm.tmp) %>%
    dplyr::mutate(Samples = colnames(df.input.data.GlobalNorm)[-1]) %>%
    dplyr::select(Samples, df.input.data.GlobalNorm$miRNA)

# PCA Plot
# ------------------------------------ #
fun.pca = function(comp, cols) {
    # comp is the comparison group: eg 'Comparison 1'
    # cols is the color representing group: col.compare.1
    
    # Extract Comparison Group
    df.compare = ls.compare.group[[comp]]
    
    # Prepare the Data Frame for PCA
    df.tmp = dplyr::inner_join(df.compare, df.input.data.GlobalNorm.tmp, by = 'Samples')
    
    # PCA Analysis
    pca.res = prcomp(df.tmp[, -c(1:2)])
    pca.var = pca.res$sdev^2
    pca.var = data.frame(PC = factor(1:length(pca.var)),
                         variance = pca.var) %>%
        dplyr::mutate(pct = round(variance/sum(variance) * 100, 2)) %>%
        dplyr::mutate(pct_cum = cumsum(pct))
    
    # Prepare PCA data frame for Plot
    df.pca = df.tmp[, 1:2] %>% cbind(pca.res$x[, 1:2])
    colnames(df.pca) = c('sample', 'comp', 'PC1', 'PC2')
    df.pca$sample = paste0('sample_', df.pca$sample)
    df.pca$comp = factor(df.pca$comp)
    
    p.pca = ggplot(df.pca,
                   aes(x = PC1, y = PC2, color = comp, label = sample)) +
        geom_point(size = 5) +
        xlab(paste("PC1(", pca.var$pct[1], "%)")) +
        ylab(paste("PC2(", pca.var$pct[2], "%)")) +
        scale_color_manual(values = cols) +
        theme_classic() +
        ggtitle(paste("PCA for",comp))
    p.pca$labels$colour = comp
    
    return(p.pca)
}

# Hierarchical Cluster Analysis
# ------------------------------------ #
colors.heatmap = colorRampPalette(col.heatmap)(255)
fun.heatmap = function(df, main, anno.df, anno.col) {
    # df: data frame for heatmap [column:sample, row: miRNA]
    # main: the Caption of the heatmap
    # anno.df: the annotation data frame
    # anno.col: the annotation color
    p.heatmap = pheatmap(df,
                         scale = 'row',
                         show_rownames = ifelse(nrow(df) > 20, FALSE, TRUE),
                         border_color = ifelse(nrow(df) > 20, FALSE, TRUE),
                         col = colors.heatmap,
                         annotation_col = anno.df,
                         annotation_names_col = FALSE,
                         annotation_colors = anno.col,
                         main = paste("Clustergram analysis for", main))
    
    return(p.heatmap)
}

# Heatmap For All miRNA
# ------------------------------------- #
# Prepare the data frame for heatmap
allmiRNA = df.input.data.GlobalNorm[, -1] %>% as.data.frame()
index = which(apply(allmiRNA, 1, sd) == 0)
if(length(index) != 0) {
    allmiRNA = allmiRNA[-index, ]
}


rownames(allmiRNA) = df.input.data.GlobalNorm$miRNA


# Prepare annotation data frame for heatmap
df.anno.all = df.samplesheet[, c('Unique Sample ID', comparisons[comparisons != "Comparison 0"])]
df.anno.all = df.anno.all[df.anno.all$`Unique Sample ID` %in% colnames(df.input.data.GlobalNorm)[-1], ] %>%
    tibble::column_to_rownames('Unique Sample ID')


# Prepare annotation color for heatmap
color.anno.all = list('Comparison 1' = c(A = col.compare.1[1], B = col.compare.1[2]),
                      'Comparison 2' = c(A = col.compare.2[1], B = col.compare.2[2]),
                      'Comparison 3' = c(A = col.compare.3[1], B = col.compare.3[2]))
color.anno.all = color.anno.all[comparisons[comparisons != "Comparison 0"]]

# Heatmap for sample-sample distance
# ------------------------------------- #
sample.dist = dist(df.input.data.GlobalNorm.tmp[, -1])
sample.dist.matrix = as.matrix(sample.dist) 



