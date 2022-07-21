# ============================================================================= #
# Input:
    # df.input.data
    # cutoff.max
    # cutoff.min

# output:
    # fig.filter.miRNA
# ============================================================================= #
index.sp = which(df.input.data$miRNA %in% name.sp)
index.ipc = which(df.input.data$miRNA %in% name.ipc)


# Create reordered input data: row is miRNA, column is sample
df.input.data.reordered = df.input.data[-c(index.ipc, index.sp), ]
df.input.data.reordered = df.input.data.reordered[order(apply(df.input.data.reordered[, -1],
                                                              1,
                                                              mean,
                                                              na.rm = TRUE),
                                                        decreasing = FALSE),]


# Overview info, box plot
x = reshape2::melt(df.input.data.reordered)
colnames(x) = c('miRNA', 'sample', 'Ct')
x$sample = paste0('sample_', x$sample)
x$miRNA = factor(x$miRNA, levels = unique(x$miRNA))

mean.miRNA = apply(df.input.data.reordered[, -1], 1, mean, na.rm = TRUE)
df.mean.miRNA = data.frame(miRNA = df.input.data.reordered$miRNA,
                           mean.ct = mean.miRNA,
                           order = 1:length(mean.miRNA))
df.mean.miRNA$miRNA = factor(df.mean.miRNA$miRNA, levels = unique(df.mean.miRNA$miRNA))
rownames(df.mean.miRNA) = NULL
x$order = rep(1:length(mean.miRNA), num.sample)


p.filtermiRNA = ggplot() +
    geom_boxplot(data = x,
                 aes(x = miRNA, y = Ct),
                 color = col.grey,
                 outlier.color = col.others[1]) +
    geom_point(data = df.mean.miRNA,
               aes(x = miRNA, y = mean.ct),
               color = col.others[2],
               size = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    geom_hline(yintercept = c(cutoff.max), 
               linetype = "dashed",
               color = col.compare.2) +
    xlab("miRNA") +
    ylab("Ct value") +
    ggtitle("Global Ct values")

hline = function(y = 0, color = 'black') {
    list(type = 'line',
         x0 = 0,
         x1 = 1,
         xref = "paper",
         y0 = y,
         y1 = y,
         line = list(color = color, dash = 'dot'))
}
fig.filtermiRNA = plot_ly(data = x,
                            x = ~ miRNA, y = ~ Ct, split = ~ miRNA, text = ~ sample,
                            type = 'box',
                            marker = list(color = col.others[1]),
                            line = list(color = col.grey, width = 1.5),
                            showlegend = FALSE) %>%
    add_trace(data = df.mean.miRNA,
              x = ~ miRNA, y = ~ mean.ct,
              type = 'scatter',
              mode = 'markers+line',
              line = list(color = col.others[2]),
              marker = list(color = col.others[2], size = 2),
              showlegend = FALSE,
              inherit = F) %>%
    add_annotations(x = df.mean.miRNA$miRNA[round(length(mean.miRNA)/2)], 
                    y = min(x$Ct, na.rm = TRUE)-10,
                    showarrow = FALSE,
                    text = "miRNA",
                    font = list(size = 14)) %>%
    layout(title = "Global Ct values",
           shapes = list(hline(cutoff.max, col.compare.2[2])),
           xaxis = list(categoryarray = ~order, 
                        visible = FALSE),
           yaxis = list(title = "Ct value",
                        zeroline = FALSE))
