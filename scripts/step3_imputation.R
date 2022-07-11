# =========================================================================== #
# Step 3 Impute Missing Value (<= 10%, analyze, NA = minimum)
# --------------------------------------------------------------------------- #
# Input: 
    # df.input.data.SPnorm

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
num.noNA = num.sample - apply(df.input.data.SPnorm[, -1], 1, function(x) sum(is.na(x)))
names(num.noNA) = df.input.data.SPnorm$miRNA
num.noNA = num.noNA[order(num.noNA, decreasing = TRUE)]
df.num.noNA = data.frame(miRNA = factor(names(num.noNA), levels = names(num.noNA)),
                         nonMissing = num.noNA,
                         status = ifelse(num.noNA >= num.sample* (1 - threshold.impute),'analyze', 'not analyze'))
p.num.noNA = ggplot(df.num.noNA,
                    aes(x = miRNA, 
                        y = nonMissing,
                        fill = status)) +
    geom_col(width = 1) +
    scale_fill_manual(values = c(col.others[2], col.grey)) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    ggtitle('Number of non-missing values of miRNA across samples') +
    ylab('Number of non-missing values of each miRNA')

