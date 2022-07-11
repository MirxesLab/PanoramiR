# =========================================================================== #
# Step 3 Global Normalization
# --------------------------------------------------------------------------- #
# Input:
    # df.input.data.noNA [row: 376 miRNA, column: miRNA, 1, 2, ...]
    
# Output:
    # df.input.data.GlobalNorm [row: 376 miRNA, column: miRNA, 1, 2, ...]
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
for (i in 1:num.sample) {
    df.input.data.GlobalNorm[, i+1] = df.input.data.noNA[, i+1] - GN.Factor[i]
}


# Save Results
# -------------------------- #
df.input.data.GlobalNorm.tmp = df.input.data.GlobalNorm
colnames(df.input.data.GlobalNorm.tmp) = c('miRNA', 
                                           paste0('sample_', 
                                                  colnames(df.input.data.GlobalNorm.tmp)[-1]))
write.table(df.input.data.GlobalNorm.tmp,
            file = file.path(dir.out.tbl, "Global.Norm.tsv"),
            sep = '\t',
            row.names = FALSE)
