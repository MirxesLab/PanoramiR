# =========================================================================== #
# Step 5 Differential Expression Analysis
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.GlobalNorm
    # ls.compare.group
# params:
    # threshold.DE.pvalue = 0.05   
    # threshold.DE.dCt = 1   
# output:
    # fun.DEanalysis(comp)
        # groupA     = df.miRNA.A [miRNA, samples]
        # groupB     = df.miRNA.B [miRNA, samples]
        # res.ori    = df.res.ori 
        # res.order  = df.res.order  [Save] [Ordered by p value]
        # res.filter = df.res.filter [Save] 

    # fun.plot.dCt(ls.diff, comp) # based on pvalue
    # fun.plot.scatter(ls.diff, comp) # based on dCt
    # fun.plot.volcano(ls.diff, comp) # based on both pvalue and dCt
# =========================================================================== #

# Generate data frame of differential expression analysis
fun.DEanalysis = function(comp) {
    # comp is the comparison group. eg: 'Comparison 1'
    
    ls.diff = list()
    
    # Group A is the test sample, Group B is the control sample
    df.compare = ls.compare.group[[comp]]
    groupA.sample = df.compare$Samples[df.compare[comp] == 'A']
    groupB.sample = df.compare$Samples[df.compare[comp] == 'B']
    
    df.miRNA.A = df.input.data.GlobalNorm[, c('miRNA', groupA.sample)]
    df.miRNA.B = df.input.data.GlobalNorm[, c('miRNA', groupB.sample)]
    
    mean.A = apply(df.miRNA.A[, -1], 1, mean)
    mean.B = apply(df.miRNA.B[, -1], 1, mean)
    
    sd.A = apply(df.miRNA.A[, -1], 1, sd)
    sd.B = apply(df.miRNA.A[, -1], 1, sd)
    
    pvalue = c() # Student's T test
    for (i in 1:num.miRNA.noNA) {
        pvalue.tmp = t.test(df.miRNA.A[i, -1], df.miRNA.B[i, -1])$p.value
        pvalue = append(pvalue, pvalue.tmp)
    }
    
    df.res.ori = data.frame(miRNA = df.input.data.GlobalNorm$miRNA,
                            mean.A = round(mean.A, 3),
                            sd.A = round(sd.A, 3),
                            mean.B = round(mean.B, 3),
                            sd.B = round(sd.B, 3),
                            dCt = round(mean.A - mean.B, 3),
                            abs.dCt = round(abs(mean.A - mean.B), 3),
                            pvalue = round(pvalue, 5))
    
    df.res.order = df.res.ori[order(df.res.ori$pvalue, decreasing = FALSE), ] %>%
        dplyr::select(miRNA, mean.A, sd.A, mean.B, sd.B, dCt, pvalue)
    rownames(df.res.order) = NULL
    
    # Filter based on p value and dCt
    df.res.filter = df.res.ori %>%
        dplyr::filter(pvalue <= threshold.DE.pvalue & abs.dCt >= threshold.DE.dCt) %>%
        dplyr::select(miRNA, mean.A, sd.A, mean.B, sd.B, dCt, pvalue)
    df.res.filter = df.res.filter[order(df.res.filter$pvalue, 
                                        decreasing = FALSE), ]
    rownames(df.res.filter) = NULL
    
    ls.diff$groupA     = df.miRNA.A
    ls.diff$groupB     = df.miRNA.B
    ls.diff$res.ori    = df.res.ori
    ls.diff$res.order  = df.res.order
    ls.diff$res.filter = df.res.filter
    
    return(ls.diff)
}

# Generate Figures
cols.sig.DE = c("sig.up" = col.others[1], "sig.down" = col.others[2], "ns" = col.grey)
cols.DE = c("up" = col.others[1], "down" = col.others[2], "ns" = col.grey)

## bar plot of dCt: filter based on pvalue
fun.plot.dCt    = function(ls.diff, comp) {
    df.res.order = ls.diff$res.order
    
    df.tmp = df.res.order %>%
        dplyr::filter(pvalue <= threshold.DE.pvalue ) %>%
        dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt ~ 'sig.up',
                                       dCt <= -threshold.DE.dCt ~ 'sig.down',
                                       TRUE ~ 'ns'))
    
    df.tmp = df.tmp[order(df.tmp$dCt, decreasing = FALSE), ]
    df.tmp$miRNA = factor(df.tmp$miRNA, levels = df.tmp$miRNA)
    
    p.dCt = ggplot(df.tmp, 
                   aes(x = miRNA, y = dCt, fill = type)) +
        geom_bar(stat = 'identity', width = 0.5) +
        scale_fill_manual(values = cols.sig.DE,
                           breaks = c('sig.up', 'sig.down', 'ns')) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5)) +
        ylab("dCt Values") +
        geom_hline(yintercept = c(threshold.DE.dCt, -threshold.DE.dCt), linetype = 'dashed') +
        ggtitle(paste0("Differentially expressed miRNA ( ", comp, " )"))
    
    return(p.dCt)    
}

## Scatter plot: Filter based on dCt
fun.plot.scatter = function (ls.diff, comp) {
    df.res.ori = ls.diff$res.ori
    df.tmp = df.res.ori %>% 
        dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt  ~ 'up',
                                       dCt <= -threshold.DE.dCt ~ 'down',
                                       TRUE ~ 'ns'))
    
    p.scatter = ggplot(data = df.tmp,
                       aes(x = mean.B, y = mean.A, color = type, label = miRNA)) +
        geom_point(size = 1) +
        geom_abline(slope = 1, intercept = 0, color = col.grey) + 
        geom_abline(slope = 1, intercept = threshold.DE.dCt, color = col.others[1]) +
        geom_abline(slope = 1, intercept = -threshold.DE.dCt, color = col.others[2]) +
        theme_classic() +
        scale_color_manual(values = cols.DE,
                           breaks = c("up", "down", "ns")) +
        labs(x = "Control Ct value", 
             y = "Treatment Ct value",
             title = paste0(comp, ', dCt cut-off = ', threshold.DE.dCt)) +
        theme_classic() +
        theme(legend.position = "top")
    
    return(p.scatter)
}

## volcano plot, filter based on both p value and dCt
fun.plot.volcano = function(ls.diff, comp) {
    df.res.ori = ls.diff$res.ori
    
    df.tmp = df.res.ori %>% 
        dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt & pvalue < threshold.DE.pvalue ~ 'sig.up',
                                       dCt <= -threshold.DE.dCt & pvalue < threshold.DE.pvalue ~ 'sig.down',
                                       TRUE ~ 'ns'))
    
    p.volcano = ggplot(df.tmp,
                       aes(x = dCt, y = -log10(pvalue), col = type, label = miRNA)) +
        geom_point(size = 3) +
        scale_color_manual(values = cols.sig.DE,
                           breaks = c("sig.up", "sig.down", "ns")) +
        geom_hline(yintercept = -log10(threshold.DE.pvalue), linetype = "dashed") +
        geom_vline(xintercept = c(-threshold.DE.dCt, threshold.DE.dCt), linetype = "dashed") +
        labs(x = "dCt value", 
             y = "- Log10 (p value)",
             title = comp) +
        theme_classic() +
        theme(legend.position = "top")
    return(p.volcano)
}


    

