# =========================================================================== #
# Step 10 Differential Expression Analysis
# --------------------------------------------------------------------------- #
# input:
    # df.input.data.GlobalNorm

# params:
    # threshold.DE.pvalue = 0.05   
    # threshold.DE.dCt = 1   
# output:
    # fun.DEanalysis(comp)
        # groupA     = df.miRNA.A [miRNA, samples]
        # groupB     = df.miRNA.B [miRNA, samples]
        # res.order  = df.res.order  [Save] [Ordered by p value]
        # res.filter = df.res.filter [Save]
        # fdr        = FALSE / TRUE

    # fun.plot.dCt(ls.diff, comp) # based on pvalue
    # fun.plot.scatter(ls.diff, comp) # based on dCt
    # fun.plot.volcano(ls.diff, comp) # based on both pvalue and dCt
# =========================================================================== #

# update ls.compare.group
ls.compare.group <- list()
n.comp.group = length(comparisons[comparisons != 'Comparison 0'])
for (i in 1:sum(n.comp.group)) {
    df.tmp = df.samplesheet[,c(1,(1+i))]  # sample, comparison
    df.tmp = df.tmp[df.tmp$Samples %in% colnames(df.input.data.GlobalNorm), ]
    df.tmp = df.tmp %>%
        dplyr::filter(!is.na(df.tmp[2]))
    
    ls.compare.group[[i]] <- df.tmp    
}
names(ls.compare.group) <- comparisons[comparisons != 'Comparison 0']


# Generate data frame of differential expression analysis
fun.summary = function(df.group, group) {
    mean   = apply(df.group[, -1],1, mean, na.rm = TRUE)
    sd     = apply(df.group[, -1], 1, sd, na.rm = TRUE)
    
    df.summary = data.frame(miRNA = df.group$miRNA,
                            mean = round(mean,2),
                            sd = round(sd,2))
    
    colnames(df.summary) = c('miRNA', 
                             paste0(c('mean.','sd.'), group))
    
    return(df.summary)
}

# DE analysis
fun.DEanalysis.Ttest = function(comp) {
    # comp is the comparison group. eg: 'Comparison 1'
    
    ls.diff = list()
    
    # Group A is the test sample, Group B is the control sample
    # dCt = ctrl - exp dCt < 0: downregulated; dCt > 0: upregulated
    df.compare = ls.compare.group[[comp]]
    groupA.sample = df.compare$Samples[df.compare[comp] == 'A']
    groupB.sample = df.compare$Samples[df.compare[comp] == 'B']
    
    df.miRNA.A = df.input.data.GlobalNorm[, c('miRNA', groupA.sample)]
    df.miRNA.B = df.input.data.GlobalNorm[, c('miRNA', groupB.sample)]
    
    df.summary.A = fun.summary(df.miRNA.A, 'A')
    df.summary.B = fun.summary(df.miRNA.B, 'B')
    df.summary   = dplyr::inner_join(df.summary.A,
                                     df.summary.B,
                                     by = 'miRNA')
    dCt = df.summary$mean.B - df.summary$mean.A
    
    num.miRNA = nrow(df.input.data.GlobalNorm)
    pvalue = c() # Student's T test
    for (i in 1:num.miRNA) {
        pvalue.tmp = t.test(df.miRNA.A[i, -1] %>% as.numeric(), 
                            df.miRNA.B[i, -1] %>% as.numeric())$p.value
        pvalue = append(pvalue, pvalue.tmp)
    }
    names(pvalue) = df.input.data.GlobalNorm$miRNA
    adj.pvalue = p.adjust(pvalue, method = 'BH')
    
    df.res.ori = data.frame(miRNA = df.input.data.GlobalNorm$miRNA ,
                            dCt = dCt,
                            pvalue = pvalue,
                            adj.pvalue = adj.pvalue) %>%
        cbind(df.summary[, -1])
    
    df.res.order = df.res.ori[order(df.res.ori$adj.pvalue, df.res.ori$pvalue, decreasing = FALSE), ] 
    rownames(df.res.order) = NULL
    
    # Filter based on p value and dCt
    df.res.filter = df.res.order %>%
        dplyr::filter(adj.pvalue <= threshold.DE.pvalue & abs(dCt) >= threshold.DE.dCt)
    fdr = TRUE
    
    if(nrow(df.res.filter) == 0) {
        df.res.filter = df.res.order %>%
            dplyr::filter(pvalue <= threshold.DE.pvalue & abs(dCt) >= threshold.DE.dCt)
        fdr = FALSE
    }
    rownames(df.res.filter) = NULL
    
    
    ls.diff$groupA     = df.miRNA.A
    ls.diff$groupB     = df.miRNA.B
    ls.diff$res.order  = df.res.order
    ls.diff$res.filter = df.res.filter
    ls.diff$fdr        = fdr
    
    return(ls.diff)
}

fun.DEanalysis.limma = function(comp) {
    ls.diff = list()
    
    df.compare = ls.compare.group[[comp]]
    df.compare = ls.compare.group[[comp]]
    groupA.sample = df.compare$Samples[df.compare[comp] == 'A']
    groupB.sample = df.compare$Samples[df.compare[comp] == 'B']
    
    df.miRNA.A = df.input.data.GlobalNorm[, c('miRNA', groupA.sample)]
    df.miRNA.B = df.input.data.GlobalNorm[, c('miRNA', groupB.sample)]
    
    df.summary.A = fun.summary(df.miRNA.A, 'A')
    df.summary.B = fun.summary(df.miRNA.B, 'B')
    df.summary  = inner_join(df.summary.A, df.summary.B, by = 'miRNA')
    
    # metadata 
    meta = df.compare
    colnames(meta) = c('sample', 'condition')
    meta$condition = factor(meta$condition, levels = c('A', 'B'))
    
    # contrast matrix
    design = model.matrix(~ 0 + meta$condition)
    colnames(design) = c('A', 'B')
    contrast_martix = makeContrasts(compare = B - A,
                                    levels = design)
    df.miRNA = df.input.data.GlobalNorm %>%
        tibble::column_to_rownames('miRNA')
    
    # DE analysis with limma
    fit = lmFit(df.miRNA, design)
    fit2 = contrasts.fit(fit, contrasts = contrast_martix)
    ebayes = eBayes(fit2)
    
    # result
    name.sum = colnames(df.summary)
    df.res.order = topTable(ebayes, number = Inf) %>% 
        tibble::rownames_to_column('miRNA') %>%
        dplyr::select('miRNA', 'logFC', 'P.Value', 'adj.P.Val')
    df.res.order = inner_join(df.res.order, df.summary, by = 'miRNA')
    colnames(df.res.order) = c('miRNA', 'dCt', 'pvalue', 'adj.pvalue', name.sum[-1])

    df.res.filter = df.res.order %>%
        dplyr::filter(adj.pvalue <= threshold.DE.pvalue & abs(dCt) >= threshold.DE.dCt)
    rownames(df.res.filter) = NULL
    
    fdr = TRUE
    
    ls.diff$groupA     = df.miRNA.A
    ls.diff$groupB     = df.miRNA.B
    ls.diff$res.order  = df.res.order
    ls.diff$res.filter = df.res.filter
    ls.diff$fdr        = fdr
    
    return(ls.diff)
}

# Generate Figures
cols.sig.DE = c("sig.up" = col.others[1], "sig.down" = col.others[2], "ns" = col.grey)
cols.DE = c("up" = col.others[1], "down" = col.others[2], "ns" = col.grey)

## bar plot of dCt: filter based on pvalue
fun.plot.dCt    = function(ls.diff, comp) {
    df.res.order = ls.diff$res.order
    fdr = ls.diff$fdr
    
    if(fdr) {
        df.tmp = df.res.order %>%
            dplyr::filter(adj.pvalue <= threshold.DE.pvalue ) %>%
            dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt ~ 'sig.up',
                                           dCt <= -threshold.DE.dCt ~ 'sig.down',
                                           TRUE ~ 'ns'))
    } else {
        df.tmp = df.res.order %>%
            dplyr::filter(pvalue <= threshold.DE.pvalue ) %>%
            dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt ~ 'sig.up',
                                           dCt <= -threshold.DE.dCt ~ 'sig.down',
                                           TRUE ~ 'ns'))
    }

    
    df.tmp = df.tmp[order(df.tmp$dCt, decreasing = TRUE), ]
    df.tmp$miRNA = factor(df.tmp$miRNA, levels = df.tmp$miRNA)
    
    if(length(df.tmp$miRNA) > 40) {
        p.dCt = ggplot(df.tmp, 
                       aes(x = miRNA, y = dCt, fill = type)) +
            geom_bar(stat = 'identity', width = 0.5) +
            scale_fill_manual(values = cols.sig.DE,
                              breaks = c('sig.up', 'sig.down', 'ns')) +
            theme_classic() +
            theme(axis.text.x =element_blank())+
            ylab("dCt Values") +
            geom_hline(yintercept = c(threshold.DE.dCt, -threshold.DE.dCt), linetype = 'dashed') +
            ggtitle(paste0("Differentially expressed miRNA ( ", comp, " )")) 
    } else {
        p.dCt = ggplot(df.tmp, 
                       aes(x = miRNA, y = dCt, fill = type)) +
            geom_bar(stat = 'identity', width = 0.5) +
            scale_fill_manual(values = cols.sig.DE,
                              breaks = c('sig.up', 'sig.down', 'ns')) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 270, hjust = 1, vjust = 0.5))+
            ylab("dCt Values") +
            geom_hline(yintercept = c(threshold.DE.dCt, -threshold.DE.dCt), linetype = 'dashed') +
            ggtitle(paste0("Differentially expressed miRNA ( ", comp, " )")) 
    }
    return(p.dCt)    
}

## Scatter plot: Filter based on dCt
fun.plot.scatter = function (ls.diff, comp) {
    df.order = ls.diff$res.order
    df.tmp = df.order %>% 
        dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt  ~ 'up',
                                       dCt <= -threshold.DE.dCt ~ 'down',
                                       TRUE ~ 'ns'))
    
    p.scatter = ggplot(data = df.tmp,
                       aes(x = mean.A, y = mean.B, color = type, label = miRNA)) +
        geom_point(size = 1) +
        geom_abline(slope = 1, intercept = 0, color = col.grey) + 
        geom_abline(slope = 1, intercept = threshold.DE.dCt, color = col.others[1]) +
        geom_abline(slope = 1, intercept = -threshold.DE.dCt, color = col.others[2]) +
        theme_classic() +
        scale_color_manual(values = cols.DE,
                           breaks = c("up", "down", "ns")) +
        labs(x = "Treatment Ct value", 
             y = "Control Ct value",
             title = paste0(comp, ', dCt cut-off = ', threshold.DE.dCt)) +
        theme_classic() +
        theme(legend.position = "top")
    
    return(p.scatter)
}

## volcano plot, filter based on both p value and dCt
fun.plot.volcano = function(ls.diff, comp) {
    df.res.order = ls.diff$res.order
    fdr          = ls.diff$fdr
    if(fdr) {
        df.tmp = df.res.order%>% 
            dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt & adj.pvalue < threshold.DE.pvalue ~ 'sig.up',
                                           dCt <= -threshold.DE.dCt & adj.pvalue < threshold.DE.pvalue ~ 'sig.down',
                                           TRUE ~ 'ns'))
    } else {
        df.tmp = df.res.order%>% 
            dplyr::mutate(type = case_when(dCt >= threshold.DE.dCt & pvalue < threshold.DE.pvalue ~ 'sig.up',
                                           dCt <= -threshold.DE.dCt & pvalue < threshold.DE.pvalue ~ 'sig.down',
                                           TRUE ~ 'ns'))
    }
    
    
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


# # check
# ls.diff.1.ttest = fun.DEanalysis.Ttest('Comparison 1')
# ls.diff.1.limma = fun.DEanalysis.limma('Comparison 1')
# 
# ls.diff.1.ttest$res.filter %>% write_xlsx('tmp/ttest.PanoramiR.xlsx')
# ls.diff.1.limma$res.filter %>% write_xlsx('tmp/limma.PanoramiR.xlsx')
# 
# all.ttest = ls.diff.1.ttest$res.order
# all.limma = ls.diff.1.limma$res.order
# res.ttest = ls.diff.1.ttest$res.filter
# res.limma = ls.diff.1.limma$res.filter
# 
# miRNA = df.input.data.GlobalNorm[df.input.data.GlobalNorm$miRNA %in% res.limma$miRNA, ]
# 
# # diff.ttest = all.ttest[all.ttest$miRNA %in% diff$miRNA, ]
# # 
# # 
# # test.limma = all.limma[all.limma$miRNA %in% a, ] %>% write_xlsx('tmp/test.xlsx')
# # 
# # all.ttest[all.ttest$miRNA == 'hsa-miR-181c-3p', ]
# # A = df.input.data.GlobalNorm[df.input.data.GlobalNorm$miRNA == 'hsa-miR-495-3p', c(2:4)] %>% as.numeric()
# # B = df.input.data.GlobalNorm[df.input.data.GlobalNorm$miRNA == 'hsa-miR-495-3p', c(5:7)] %>% as.numeric()
# # mean(A) - mean(B)
# # 
#  filter.1 = res.ttest$sd.A - res.ttest$sd.B
#  filter.2 = res.limma$sd.A - res.limma$sd.B
# 
# df.sd = data.frame(sd = filter.1, group = 't-test') %>%
#     rbind(data.frame(sd = filter.2, group = 'limma'))
# ggplot(df.sd, aes(x = group, y = sd, fill = group))+
#     geom_boxplot()
