# =========================================================================== #
# Step 4 Filter and SP normalization
# --------------------------------------------------------------------------- #
# input
    # df.input.data / df.input.data.ipc
    # cutoff.max
    # cutoff.min
    # file.ref.miRthld [only for PanoramiR]

# Output:
    # min.up.limit [Only for PanoramiR]
    # df.input.data.Filt1
    # fig.filtermiRNA
# =========================================================================== #

                   # -------------------------------- #
                   #        First Round Filter        #
                   # -------------------------------- #
# Prepare, PanoramiR do not have Verification check step
if(arg.pipeline == 'PanoramiR') {df.input.data.ipc = df.input.data}

# for PanoramiR, cutoff.max is min.up.limit ~ cutoff.max (by default: 33)
if(arg.pipeline == 'PanoramiR') {
    threshold = read_excel(file.ref.miRthld)
    min.up.limit = min(threshold$Threshold, na.rm = TRUE)
    
    # Make sure the row order of miRNA is same as the one of df.input.data
    threshold = as.data.frame(threshold)
    rownames(threshold) = threshold$miRNA
    threshold = threshold[df.input.data.ipc$miRNA, ]
    
    # Change the cutoff.max based on the reference file
    threshold[, 2] = sapply(threshold[, 2], 
                            function(x) ifelse(x > cutoff.max, cutoff.max, x))
} else {
    threshold = data.frame(miRNA = df.input.data.ipc$miRNA,
                           Threshold = cutoff.max)
}

# Filter based on cutoff.max, cutoff.min
df.input.data.Filt1 = df.input.data.ipc
num.sample = ncol(df.input.data.ipc) - 2
for (i in 1:nrow(df.input.data.ipc)) {
    for (j in 1:num.sample) {
        bloon.1 = is.na(df.input.data.Filt1[i, j+2])
        bloon.2 = df.input.data.Filt1[i,j+2] > threshold[i, 2]
        bloon.3 = df.input.data.Filt1[i,j+2] < cutoff.min

        if (bloon.1[1] | bloon.2[1] | bloon.3[1] ) {
            df.input.data.Filt1[i,j+2] = NA
        } 
    }
}

                # -------------------------------- #
                #  Visualization (Cancer/Biofluid) #
                # -------------------------------- #
# This plot do not convey much information. Can be removed later
if(arg.pipeline != 'PanoramiR') {
    
    index.sp = which(df.input.data.ipc$miRNA %in% name.sp)
    index.ipc = which(df.input.data.ipc$miRNA %in% name.ipc)
    
    # Create reordered input data
    df.input.data.reordered = df.input.data.ipc[-c(index.ipc, index.sp), -1]
    df.input.data.reordered = df.input.data.reordered[order(apply(df.input.data.reordered[, -c(1,2)],
                                                                  1,
                                                                  mean,
                                                                  na.rm = TRUE),
                                                            decreasing = FALSE),]
    
    
    # Overview info, box plot
    x = reshape2::melt(df.input.data.reordered, id.vars = 'miRNA')
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
    
    # This Plot is for saving
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
                   color = col.others[1]) +
        geom_hline(yintercept = c(cutoff.sp), 
                   linetype = "dashed",
                   color = col.others[1]) +
        xlab("miRNA") +
        ylab("Ct value") +
        ggtitle("Global Ct values")
    
    # This Plot is for displaying
    hline = function(y = 0, color = 'black', cutoff = cutoff.max) {
        list(type = 'line',
             x0 = 0,
             x1 = 1,
             xref = "paper",
             y0 = y,
             y1 = y,
             line = list(color = color, dash = 'dot'),
             title = cutoff)
    }
    fig.filtermiRNA = plot_ly(data = x,
                              x = ~ miRNA, y = ~ Ct, split = ~ miRNA, text = ~ sample,
                              type = 'box',
                              boxpints = 'outliers',
                              marker = list(color = col.others[1]),
                              line = list(color = col.grey, width = 1.5),
                              showlegend = FALSE) %>%
        add_trace(data = df.mean.miRNA,
                  x = ~ miRNA, y = ~ mean.ct,
                  type = 'scatter',
                  mode = 'markers',
                  line = list(color = col.others[2]),
                  marker = list(color = col.others[2], size = 2),
                  showlegend = FALSE,
                  inherit = F) %>%
        add_annotations(x = df.mean.miRNA$miRNA[round(length(mean.miRNA)/2)], 
                        y = min(x$Ct, na.rm = TRUE)-10,
                        showarrow = FALSE,
                        text = "miRNA",
                        font = list(size = 14)) %>%
        add_annotations(x = df.mean.miRNA$miRNA[50],
                        y = cutoff.max + 1,
                        showarrow = F,
                        text = paste0('Round1 Filter: ', cutoff.max),
                        font = list(size = 12)) %>%
        add_annotations(x = df.mean.miRNA$miRNA[50],
                        y = cutoff.sp - 1,
                        showarrow = F,
                        text = paste0('Round2 Filter: ', cutoff.sp),
                        font = list(size = 12)) %>%
        layout(title = "Global Ct values",
               shapes = list(hline(cutoff.max, col.others[1], cutoff.max),
                             hline(cutoff.sp, col.others[1]), cutoff.sp),
               xaxis = list(categoryarray = ~order, 
                            visible = FALSE),
               yaxis = list(title = "Ct value",
                            zeroline = FALSE))
}






