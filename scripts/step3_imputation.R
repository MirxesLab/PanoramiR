# =========================================================================== #
# Step 3 Impute Missing Value (<= 10%, analyze, NA = minimum)
# --------------------------------------------------------------------------- #
# Input: 
    # df.input.data.SPnorm     [Step 2]

# Params:
    # threshold.impute    = 0.1

# Output:
    # df.input.data.noNA
    # num.miRNA.noNA
    # p.num.noNA
# =========================================================================== #

# Remove miRNA with the number of missing values more than 10% number of sample
index.remove = which(apply(df.input.data.SPnorm, 
                           1, 
                           function(x) length(which(is.na(x)))) > threshold.impute * num.sample)

df.input.data.noNA = df.input.data.SPnorm[-index.remove, ]

# Impute missing value based on the minimum value of miRNA across samples
min.miRNA = as.numeric(apply(df.input.data.noNA, 1, min, na.rm = TRUE))
num.miRNA.noNA = length(min.miRNA)
for (i in 1:num.miRNA.noNA) {
    df.input.data.noNA[i, -1] = replace(df.input.data.noNA[i, -1],
                                      which(is.na(df.input.data.noNA[i, -1])),
                                      min.miRNA[i])
}


# Visualization
num.miRNA = nrow(df.input.data.SPnorm)
num.NA = apply(df.input.data.SPnorm[, -1], 1, function(x) sum(is.na(x)))
num.include = sum(num.NA <= num.sample * threshold.impute)
num.exclude = sum(num.NA > num.sample * threshold.impute)
num.summary = c(num.include, num.exclude)
df.miss = data.frame('Missing values across samples' = c(paste('<=',floor(num.sample * threshold.impute)),
                                                         paste('>', floor(num.sample * threshold.impute))),
                     'Number of miRNA' = num.summary,
                     'Percent' = paste0(round(num.summary/num.miRNA *100, 2), '%'),
                     status = c('include', 'exclude'))
rownames(df.miss) = NULL




# num.noNA = num.sample - apply(df.input.data.SPnorm[, -1], 1, function(x) sum(is.na(x)))
# names(num.noNA) = df.input.data.SPnorm$miRNA
# num.noNA = num.noNA[order(num.noNA, decreasing = TRUE)]
# df.num.noNA = data.frame(miRNA = factor(names(num.noNA), levels = names(num.noNA)),
#                          nonMissing = num.noNA,
#                          status = ifelse(num.noNA >= num.sample* (1 - threshold.impute),'include', 'exclude'))
# p.num.noNA = ggplot(df.num.noNA,
#                     aes(x = miRNA, 
#                         y = nonMissing,
#                         fill = status)) +
#     geom_col(width = 1) +
#     scale_fill_manual(values = c(col.grey, col.others[2])) +
#     theme_classic() +
#     theme(axis.text.x = element_blank()) +
#     ggtitle('Number of non-missing values of miRNA across samples') +
#     ylab('Number of non-missing values of each miRNA')

