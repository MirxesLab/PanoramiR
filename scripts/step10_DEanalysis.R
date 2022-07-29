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
    dCt = df.summary$mean.A - df.summary$mean.B
    
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
fun.DEanalysis.Utest = function(comp) {
    ls.diff = list()
    
    # Group A is the test sample, Group B is the control sample
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
    
    dCt = df.summary$mean.A - df.summary$mean.B
    
    
    # test.norm.dis.A = c()
    # test.norm.dis.B = c()
    # 
    # for (i in miRNA) {
    #     pvalue = shapiro.test(df.miRNA.A[df.miRNA.A$miRNA == i, -1] %>%as.numeric())
    #     pvalue = pvalue$p.value
    #     test.norm.dis.A = append(test.norm.dis.A, pvalue)
    # }
    # 
    # for (i in miRNA) {
    #     pvalue = shapiro.test(df.miRNA.B[df.miRNA.B$miRNA == i, -1] %>%as.numeric())
    #     pvalue = pvalue$p.value
    #     test.norm.dis.B = append(test.norm.dis.B, pvalue)
    # }
    # index.A = which(test.norm.dis.A<0.05)
    # index.B = which(test.norm.dis.B<0.05)
    
    miRNA = df.input.data.GlobalNorm$miRNA
    pvalue = c()
    for(i in miRNA) {
        A = df.miRNA.A[df.miRNA.A$miRNA == i, -1] %>% as.numeric()
        B = df.miRNA.B[df.miRNA.B$miRNA == i, -1] %>% as.numeric()
        p.tmp = wilcox.test(A, B, paired = FALSE, exact = TRUE)
        p.tmp = p.tmp$p.value
        names(p.tmp) = i
        pvalue = append(pvalue, p.tmp)
    }
    
    adj.pvalue = p.adjust(pvalue, method = 'BH')
    
    df.res.ori = df.summary %>%
        dplyr::mutate(dCt    = dCt %>% round(2),
                      pvalue = pvalue[df.summary$miRNA] %>% round(5),
                      adj.pvalue = adj.pvalue[df.summary$miRNA] %>% round(5))

    df.res.order = df.res.ori[order(df.res.ori$adj.pvalue, 
                                    df.res.ori$pvalue, 
                                    decreasing = FALSE), ]
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
} # remove
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
    contrast_martix = makeContrasts(compare = A - B,
                                    levels = design)
    
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

    
    df.tmp = df.tmp[order(df.tmp$dCt, decreasing = FALSE), ]
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


# Check with Jeremy
# excel.result = read_excel('tmp/Cancer384_v3.1.xlsx', 'Results', skip = 1)
# colnames(excel.result)
# 
# excel.result = excel.result[, c('miRNA ID', 'p value')]
# colnames(excel.result) = c('miRNA', 'excel.value')
# 
# check.excel = inner_join(df.res.ori, excel.result, by = 'miRNA')
# check.excel
# diff.1 = mean(check.excel$pvalue - check.excel$excel.value, na.rm = TRUE)
# diff.2 = mean(check.excel$adj.pvalue - check.excel$excel.value, na.rm = TRUE)
# write_xlsx(check.excel, 'tmp/checkpvalue.xlsx')
