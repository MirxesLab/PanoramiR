# =========================================================================== #
# Step 8 Global Normalization
# --------------------------------------------------------------------------- #
# Input:
    # df.input.data.noNA [miRNA, 1, 2, ...]
    
# Output:
    # df.input.data.GlobalNorm      -> Save
    # fun.plot.violin 
    
# =========================================================================== #

# Global Normalization
# -------------------------- #
# Calculate the mean of each sample
mean.sample = apply(df.input.data.noNA[, -1], 2, mean, na.rm = TRUE)

# Calculate the mean of all samples
mean.all = mean(mean.sample)

# Calculate the Global Norm Factor
GN.Factor = mean.sample - mean.all

# Global Normalization
df.input.data.GlobalNorm = df.input.data.noNA
num.sample = ncol(df.input.data.noNA) - 1
for (i in 1:num.sample) {
    df.input.data.GlobalNorm[, i+1] = df.input.data.noNA[, i+1] - GN.Factor[i]
}


# Visualization: miRNA Distribution
# ---------------------------------- #
fun.plot.violin = function(comp, cols) {
    ls.violin = list()
    # Data frame recording the Sample and Comparison Group
    df.compare = ls.compare.group[[comp]]
    
    groupA.sample = df.compare$Samples[df.compare[, 2] == "A"]
    groupB.sample = df.compare$Samples[df.compare[, 2] == "B"]
    
    # Preparation data frame of SP normalization CT values 
    df.input.data.GlobalNorm.tmp = t(df.input.data.GlobalNorm[, -1])
    colnames(df.input.data.GlobalNorm.tmp) = df.input.data.GlobalNorm$miRNA
    df.input.data.GlobalNorm.tmp = as.data.frame(df.input.data.GlobalNorm.tmp) %>%
        dplyr::mutate(Samples = colnames(df.input.data.GlobalNorm)[-1]) %>%
        dplyr::select(Samples, df.input.data.GlobalNorm$miRNA)
    
    if(class(df.compare$Samples) != class(df.input.data.GlobalNorm.tmp$Samples)) {
        df.compare$Samples = as.character(df.compare$Samples)
    }
    
    df.tmp = dplyr::inner_join(df.compare, df.input.data.GlobalNorm.tmp,by = "Samples")
    
    # Reorder the miRNA based on the mean of miRNA across samples
    df.tmp.t = t(df.tmp[, -c(1:2)][, order(apply(df.tmp[, -c(1:2)],
                                                 2,
                                                 mean,
                                                 na.rm = TRUE),
                                           decreasing = FALSE)])
    colnames(df.tmp.t) = df.tmp$Samples 
    
    # Prepare data frame for Violin plot
    df.tmp.melt = reshape2::melt(df.tmp.t[-c(1,2), ])
    df.tmp.melt = df.tmp.melt %>% 
        dplyr::mutate(group =  ifelse(Var2 %in% groupA.sample, 'A', 'B'))
    colnames(df.tmp.melt) = c('miRNA', 'sample', 'Ct', 'group')
    
    df.tmp.melt$sample = factor(df.tmp.melt$sample,
                                levels = unique(df.tmp.melt$sample))
    df.tmp.melt$group = factor(df.tmp.melt$group)
    
    # Violin Plot
    p.violin = ggplot(df.tmp.melt,
                      aes(x = sample, y = Ct, fill = group)) +
        geom_boxplot(width = 0.1) +
        #geom_jitter(shape = 16, position = position_jitter(0.2), color = NA) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +
        scale_fill_manual(values = alpha(cols, 0.5)) +
        ylab("Ct Values") +
        xlab("sample") +
        ggtitle(paste0("miRNA Expression Level ( ", comp, " )"))
    if(nrow(df.tmp) > 25) {
      p.violin = p.violin + geom_violin(trim = F, color = NA)
    } else {
      p.violin = p.violin + geom_violin(trim = F)
    }
  
    
    fig.violin = df.tmp.melt %>%
        plot_ly(x = ~ sample,
                y = ~ Ct,
                split = ~ sample,
                showlegend = FALSE,
                color = ~ group,
                colors = c(A = cols[1], B = cols[2]),
                box = list(visible = T),
                meanline = list(visible = T),
                type = 'violin') %>%
        layout(title = list(text = paste('miRNA expression level in each sample (',
                                         comp,
                                         ')')),
               xaxis = list(title = 'sample', tickangle = 90),
               yaxis = list(title = 'Ct Values', zeroline = F))
    ls.violin$ggplot = p.violin
    ls.violin$plotly = fig.violin
    return(ls.violin)
}

# Save Results
# -------------------------- #
write.csv(df.input.data.GlobalNorm, file.path(dir.out.tbl, 'Data GlobalNorm.csv'))



