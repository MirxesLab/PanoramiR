# =========================================================================== #
# Step 7 Imputation
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.Filt2
    # threshold.impute [by default: 0.1]

# output:
    # p.num.miRNA
    # df.input.data.noNA       
    # df.miss
    # df.exclude.miRNA       -> Save
# =========================================================================== #

                    # --------------------------------- #
                    #   Plot: Number of detected miRNA  #
                    # --------------------------------- #
df.input.data.Filt2 = df.input.data.Filt2[, -1] # Remove the 'Group' Column
num.miRNA           = nrow(df.input.data.Filt2) # All miRNA
num.NA.miRNA        = apply(df.input.data.Filt2[, -1], 
                            2, 
                            function(x) sum(is.na(x)))

df.detect = data.frame(sample    = paste0('sample_', colnames(df.input.data.Filt2)[-1]),
                       num.miRNA = num.NA.miRNA,
                       status    = 'Undetected') %>% 
    rbind(data.frame(sample = paste0('sample_', colnames(df.input.data.Filt2)[-1]),
                     num.miRNA = num.miRNA - num.NA.miRNA,
                     status = 'Detected'))

df.detect$sample = factor(df.detect$sample, levels = unique(df.detect$sample))
df.detect$status = factor(df.detect$status, levels = c("Undetected", "Detected"))

p.num.miRNA = ggplot(data = df.detect,
                     aes(x = sample, y = num.miRNA, fill = status)) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = c(col.grey, col.others[2])) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1)) +
    ggtitle('The number of miRNA detected') +
    labs(y = 'Number of miRNA')



                    # --------------------------------- #
                    #            Imputation             #
                    # --------------------------------- #

# Remove miRNA with the number of missing values more than 10% number of sample
index.remove = which(apply(df.input.data.Filt2,
                           1, 
                           function(x) sum(is.na(x))) > threshold.impute * num.sample)

if (sum(index.remove) != 0) {
    df.input.data.noNA = df.input.data.Filt2[-index.remove, ]
} else {
    df.input.data.noNA = df.input.data.Filt2
}


# Impute missing value based on the greatest Ct value of miRNA across samples
df.input.data.noNA[, -1] = apply(df.input.data.noNA[, -1], 
                                 1, 
                                 function(x) ifelse(is.na(x), 
                                                    max(x, na.rm = TRUE), 
                                                    x)) %>% t

                    # ---------------------------------- #
                    # Table: Included and Excluded miRNA #
                    # ---------------------------------- #
# Visualization
num.include = nrow(df.input.data.noNA)
num.exclude = num.miRNA - num.include
threshold.sample = floor(num.sample * threshold.impute)

df.miss = data.frame('Missing values across samples' = paste(c('<=', '>'),
                                                             threshold.sample),
                     'Number of miRNA' = c(num.include, num.exclude),
                     'Percent' = paste0(round(c(num.include, num.exclude)/num.miRNA * 100, 2),
                                        '%'),
                     status = c('include', 'exclude'))
rownames(df.miss) = NULL


                    # ---------------------------------- #
                    #           Save results             #
                    # ---------------------------------- #
df.exclude.miRNA = df.input.data.Filt2[index.remove, ]
tmp.colname      = colnames(df.exclude.miRNA)[-1]
tmp.colname      = paste('sample', tmp.colname)

df.exclude.miRNA = df.exclude.miRNA %>%
    dplyr::mutate('Missing values' = apply(df.exclude.miRNA,
                                           1,
                                           function(x) sum(is.na(x)))) %>%
    dplyr::mutate('Threshold' = threshold.sample)
colnames(df.exclude.miRNA) = c('miRNA', tmp.colname, 'Missing Values', 'Threshold')

write.csv(df.exclude.miRNA, file.path(dir.out.tbl,
                                       'Excluded miRNA (before Global Norm).csv'))

