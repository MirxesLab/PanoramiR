# =========================================================================== #
# Step 2 Spike-in Normalization(Only for PanoramiR)
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.Filt1 [Step 1]
    # file.ref.miRsp      [config file]
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
miRsp  = read_excel(file.ref.miRsp)     # [Position, miRNA ID, SP, Group]
RNAvol = read_excel(file.ref.RNAvol)    # [sample, RNA input volume]

# Check RNA input volume
# ---------------------------------- #
if(nrow(RNAvol) != num.sample) {stop()}

# Calculate SP Factor
# ---------------------------------- #
df.SP.raw = df.input.data.Filt1[grep('*.SP', df.input.data.Filt1$miRNA), ]

# If it is isolation spike-in, Spike-in should be normalized based on RNA input volume
if(!is.RTsp) {
    # Calculate Input Scaling Factor
    RNAvol = RNAvol %>% dplyr::mutate(log2Vol = log2(RNAvol$`RNA input volume`))
    mean.log2Vol = mean(RNAvol$log2Vol)
    RNAvol = RNAvol %>% dplyr::mutate(InputSF = RNAvol$log2Vol - mean.log2Vol)
    
    # Calculate RNA Volume Normalized SP
    df.SP.volNorm = df.SP.raw
    for (i in 1:num.sample) {
        df.SP.volNorm[, i+1] = df.SP.volNorm[, i+1] - RNAvol$InputSF[i]
    }
    mean.SP = apply(df.SP.volNorm[, -1], 1, mean, na.rm = TRUE)
    
    # Calculate SP factor
    df.SP.factor = apply(df.SP.volNorm[, -1], 
                         1, 
                         function(x) x - mean(x)) %>% as.data.frame()
    colnames(df.SP.factor) = df.SP.raw$miRNA

} else{
    # Calculate SP factor
    df.SP.factor = apply(df.SP.raw[, -1], 
                         1, 
                         function(x) x - mean(x)) %>% as.data.frame()
    
    colnames(df.SP.factor) = df.SP.raw$miRNA
}

# Spike-in Normalization
# ---------------------------------- #
# Prepare data frame
miRsp = miRsp[-grep('*.SP', miRsp$`miRNA ID`), ]
colnames(miRsp) = c("Posation", "miRNA", "SP", "Group")
miRsp = miRsp %>% dplyr::select('miRNA', 'SP')

# Remove Spike-in Rows
df.input.data.Filt1 = df.input.data.Filt1[-grep('*.SP', df.input.data.Filt1$miRNA), ]

# Spike-in Normalization
df.input.data.Filt1.tmp = inner_join(miRsp, df.input.data.Filt1, by = "miRNA")
for (i in 1:num.sample) {
    for (j in colnames(df.SP.factor)) {
        SF = df.SP.factor[i, j]
        Ct = df.input.data.Filt1.tmp[df.input.data.Filt1.tmp$SP==j, i+2]
        df.input.data.Filt1.tmp[df.input.data.Filt1.tmp$SP==j, i+2] = Ct - SF
    }
}
df.input.data.SPnorm = df.input.data.Filt1.tmp[, c('miRNA', as.character(1:num.sample))]

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

# Visualization: miRNA Distribution
# ---------------------------------- #
fun.plot.violin = function(comp, cols) {
    # Data frame recording the Sample and Comparison Group
    df.compare = ls.compare.group[[comp]]
    groupA.sample = df.compare$Samples[df.compare[, 2] == "A"]
    
    # Preparation data frame of SP normalization CT values 
    df.input.data.SPnorm.tmp = t(df.input.data.SPnorm[, -1])
    colnames(df.input.data.SPnorm.tmp) = df.input.data.SPnorm$miRNA
    df.input.data.SPnorm.tmp = as.data.frame(df.input.data.SPnorm.tmp) %>%
        dplyr::mutate(Samples = as.numeric(colnames(df.input.data.SPnorm)[-1])) %>%
        dplyr::select(Samples, df.input.data.SPnorm$miRNA)
    
    df.tmp = dplyr::inner_join(df.compare, df.input.data.SPnorm.tmp,by = "Samples")
    
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
    colnames(df.tmp.melt) = c('miRNA', 'sample', 'log2cp', 'group')
    
    
    df.tmp.melt$sample = paste0('sample_', df.tmp.melt$sample)
    df.tmp.melt$sample = factor(df.tmp.melt$sample,
                                levels = unique(df.tmp.melt$sample))
    df.tmp.melt$group = factor(df.tmp.melt$group)
    
    # Violin Plot
    p.violin = ggplot(df.tmp.melt,
                      aes(x = sample, y = log2cp, fill = group)) +
        geom_violin(trim = TRUE) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +
        scale_fill_manual(values = alpha(cols, 0.5)) +
        ylab("Ct Values") +
        ggtitle(paste0("miRNA Expression Levels in Each Samples ( ", comp, " )"))
    
    return(p.violin)
}


# Save Files
# ---------------------------------- #
df.input.data.SPnorm.tmp = df.input.data.SPnorm
colnames(df.input.data.SPnorm.tmp)[-1] = paste0('sample_', colnames(df.input.data.SPnorm.tmp)[-1])
write.table(df.input.data.SPnorm.tmp,
            file = file.path(dir.out.tbl, 'SpikeIn.Norm.tsv'),
            sep = '\t',
            row.names = FALSE)

df.SP.factor.tmp = df.SP.factor %>%
    dplyr::mutate(sample = paste0('sample_', 1:num.sample)) %>%
    dplyr::select(sample, colnames(df.SP.factor))
write.table(df.SP.factor.tmp,
            file = file.path(dir.out.tbl, 'SpikeIn.Factor.tsv'),
            sep = '\t',
            row.names = FALSE)
