# =========================================================================== #
# Step 2 Spike-in Normalization(Only for PanoramiR)
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.Filt1 [Step 1]    
    # file.ref.RNAvol     [config file]
    # ls.compare.group    [Step 1.1]

# params:
    # cutoff.sp [32]
    # is.RTsp

# output:
    # df.SP.factor [A1.SP, B1.SP, C1.SP ...]
    # df.input.data.SPnorm [miRNA, 1, 2, ... n, ps: no spike-in rows]
    # p.num.miRNA (Plot)
    # fun.plot.violin(comp, cols)
    # num.miRNA
# =========================================================================== #

# Read in reference files
# ---------------------------------- #
RNAvol = read_excel(file.ref.RNAvol)    # [sample, RNA input volume]

# Check RNA input volume
# ---------------------------------- #
if(nrow(RNAvol) != num.sample) {stop('Please check ther number of sample in the files')}

# Define function to calculate Factor
# ---------------------------------- #
fun.factor = function(df.x) {
    mean.sample = apply(df.x[, -c(1,2)], 2, mean, na.rm = TRUE)
    mean.all    = mean(mean.sample)
    factor  = mean.sample - mean.all
    return(factor)
}

fun.df.factor = function(df) {
    # df is the data frame recording SP or IPC Ct values
   df.factor = data.frame()
   for (i in c('A', 'B', 'C', 'D')) {
       df.tmp = df[df$Group == i, ]
       df.factor.tmp = fun.factor(df.tmp)
       if (i == 'A') {
           df.factor = df.factor.tmp
       } else {
           df.factor = rbind(df.factor, df.factor.tmp)
       }
   }
   df.factor = as.data.frame(df.factor) %>%
       dplyr::mutate(Group = c('A', 'B', 'C', 'D')) %>%
       dplyr::select('Group', colnames(df.factor))
   rownames(df.factor) = NULL
   
   return(df.factor)
}

# Calculate IPC Factor
# ---------------------------------- #
miRList = read_excel(file.ref.miRList) %>%  # [Group, Position, miRBASE v22 Accession, miRNA ID]
    dplyr::select('RT Group', 'Position', 'miRNA ID')
colnames(miRList) = c('Group', 'Position', 'miRNA')
miRList = miRList[miRList$miRNA == df.input.data.Filt1$miRNA, ]
df.input.data.Filt1.tmp = cbind(miRList[, -2], df.input.data.Filt1[, -1])

# record the IPC and SP index
index.sp = which(df.input.data.Filt1.tmp$miRNA %in% name.sp)
index.ipc = which(df.input.data.Filt1.tmp$miRNA %in% name.ipc)

# Extract IPC Ct values
df.ipc = df.input.data.Filt1.tmp[index.ipc, ]
df.sp  = df.input.data.Filt1.tmp[index.sp, ]

# Calculate the IPC factor
df.ipc.factor = fun.df.factor(df.ipc)

# IPC Normalization
df.input.data.IPC.norm = df.input.data.Filt1.tmp
for (i in c('A', 'B', 'C', 'D')) {
    df.input.data.IPC.tmp = df.input.data.IPC.norm[df.input.data.IPC.norm$Group == i,]
    df.ipc.tmp = df.ipc.factor[df.ipc.factor$Group == i, ]
    for (j in 1:num.sample) {
        df.input.data.IPC.tmp[, j+2] = df.input.data.IPC.tmp[, j+2] - df.ipc.tmp[, j+1]
    }
    df.input.data.IPC.norm[df.input.data.IPC.norm$Group == i,] = df.input.data.IPC.tmp
}

# Calculate SP Factor
# ---------------------------------- #
# If it is isolation spike-in, Spike-in should be normalized based on RNA input volume
if(!is.RTsp) {
    # Calculate Input Scaling Factor
    RNAvol = RNAvol %>% dplyr::mutate(log2Vol = log2(RNAvol$`RNA input volume`))
    mean.log2Vol = mean(RNAvol$log2Vol)
    RNAvol = RNAvol %>% dplyr::mutate(InputSF = RNAvol$log2Vol - mean.log2Vol)
    
    # Calculate RNA Volume Normalized SP
    df.SP.volNorm = df.sp
    for (i in 1:num.sample) {
        df.SP.volNorm[, i+2] = df.SP.volNorm[, i+2] - RNAvol$InputSF[i]
    }
    
    # Calculate SP factor
    df.sp.factor = fun.df.factor(df.SP.volNorm)

} else{
    # Calculate SP factor directly
    df.sp.factor = fun.df.factor(df.sp)
}

# Spike-in Normalization
# ---------------------------------- #
df.input.data.SPnorm = df.input.data.IPC.norm
for (i in c('A', 'B', 'C', 'D')) {
    df.input.data.SP.tmp = df.input.data.SPnorm[df.input.data.SPnorm$Group == i,]
    df.sp.tmp = df.sp.factor[df.sp.factor$Group == i, ]
    for (j in 1:num.sample) {
        df.input.data.SP.tmp[, j+2] = df.input.data.SP.tmp[, j+2] - df.sp.tmp[, j+1]
    }
    df.input.data.SPnorm[df.input.data.SPnorm$Group == i,] = df.input.data.SP.tmp
}

# Remove IPC & Spike-in Rows
# ---------------------------------- #
df.input.data.SPnorm = df.input.data.SPnorm[-c(index.ipc, index.sp), -1]

# Round 2 Filter After Spike-in Normalization
df.input.data.SPnorm[, -1] = apply(df.input.data.SPnorm[, -1], 
                                 2, 
                                 function(x) ifelse(x > cutoff.sp, NA, x))


# Visualization: Number of miRNA detected
# ---------------------------------- #
num.miRNA = nrow(df.input.data.SPnorm)
num.NA.miRNA = apply(df.input.data.SPnorm[, -1], 
                     2, 
                     function(x) length(which(is.na(x))))
df.miRNA = data.frame(sample = paste0('sample_', colnames(df.input.data.SPnorm)[-1]),
                      num.miRNA = num.NA.miRNA,
                      status = 'Undetected') %>%
    rbind(data.frame(sample = paste0('sample_', colnames(df.input.data.SPnorm)[-1]),
                     num.miRNA = num.miRNA - num.NA.miRNA,
                     status = 'Detected'))
df.miRNA$sample = factor(unique(df.miRNA$sample), levels = unique(df.miRNA$sample))
df.miRNA$status = factor(df.miRNA$status, levels = c("Undetected", "Detected"))

p.num.miRNA = ggplot(data = df.miRNA,
                     aes(x = sample, y = num.miRNA, fill = status)) +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = c(col.grey, col.others[2])) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 1)) +
    ggtitle('The number of miRNA detected') +
    labs(y = 'Number of miRNA')



# Save Files
# ---------------------------------- #
df.input.data.SPnorm.tmp = df.input.data.SPnorm
colnames(df.input.data.SPnorm.tmp)[-1] = paste0('sample_', colnames(df.input.data.SPnorm.tmp)[-1])
write.table(df.input.data.SPnorm.tmp,
            file = file.path(dir.out.tbl, 'SpikeIn.Norm.tsv'),
            sep = '\t',
            row.names = FALSE)

