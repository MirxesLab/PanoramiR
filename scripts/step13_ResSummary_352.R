# ============================================================================= #
# Input:
    # ls.diff 
    # ls.miRNA2gene [NA or length(ls.miRNA2gene) = 3]
    # ls.GOsum      [NA or length(ls.miRNA2gene) = 2]

# output:
    # fun.plot.summary(ls.diff, ls.miRNA2gene, ls.GOsum, comp)
# ============================================================================= #
fun.plot.summary = function(ls.diff, ls.miRNA2gene, ls.GOsum, comp) {
    # ls.diff storing the results of DE analysis
    # ls.miRNA2gene storing the results of target gene analysis
    # ls.GOsum storing the results of GO enrichment
    # comp is the 'Comparison n'
    
    # Drawing Box
    num.miRNA.all = nrow(df.input.data.Filt1)
    num.miR = num.miRNA.all - length(index.ipc) - length(index.sp)
    allmiRNA      = boxGrob(glue("ID3EAL Knowledge panel",
                                 "{num} miRNAs",
                                 num = txtInt(num.miR),
                                 .sep = '\n'))
    
    num.sample = nrow(df.samplesheet)
    
    num.miRNA.noNA = nrow(df.input.data.GlobalNorm)
    noNAmiRNA      = boxGrob(glue("{num} miRNAs detected in all samples",
                                  "sample = {sample}",
                                  num = txtInt(num.miRNA.noNA),
                                  sample = txtInt(num.sample),
                                  .sep = '\n'))
    
    num.miRNA.NA = num.miR - num.miRNA.noNA
    NAmiRNA      = boxGrob(glue("{num} miRNAs",
                                "are not detected in all samples",
                                num = txtInt(num.miRNA.NA),
                                .sep = '\n'))
    
    num.DE.miRNA = nrow(ls.diff$res.filter)
    
    num.noDE.miRNA = num.miRNA.noNA - num.DE.miRNA
    noDEmiRNA = boxGrob(glue("{num} miRNAs",
                             "do not have significant expression difference",
                             num = txtInt(num.noDE.miRNA),
                             .sep = '\n'))
    if (num.DE.miRNA != 0 ) {
        num.DE.up = sum(ls.diff$res.filter$dCt > 0)
        num.DE.down = sum(ls.diff$res.filter$dCt < 0)
        
        DEmiRNA = boxGrob(glue("{num} miRNAs have significant expression difference",
                               'p <= {pvalue}, abs(dCt) >= {dCt}',
                               '{up} miRNAs upregulated',
                               '{down} miRNAs downregulated',
                               num = txtInt(num.DE.miRNA),
                               pvalue = threshold.DE.pvalue,
                               dCt = threshold.DE.dCt,
                               up = txtInt(num.DE.up),
                               down = txtInt(num.DE.down),
                               .sep = '\n'))
    } else {
        DEmiRNA = boxGrob(glue("{num} miRNAs have significant expression difference",
                               'p < {pvalue}, abs(dCt) < {dCt}',
                               num = txtInt(num.DE.miRNA),
                               pvalue = threshold.DE.pvalue,
                               dCt = threshold.DE.dCt,
                               .sep = '\n'))
    }
    
    if(length(ls.miRNA2gene) != 1) {
        num.gene = nrow(ls.miRNA2gene$miRNA2gene.sig)
        gene = boxGrob(glue("{num} target genes are found",
                            num = txtInt(num.gene)))
    }
    
    if(length(ls.GOsum) != 1) {
        num.GO = nrow(ls.GOsum$res.save)
        GO = boxGrob(glue("{num} significantly enriched GO hits",
                          num = txtInt(num.GO)))
    }
    
    # Plot flow chart
    png(file = file.path(dir.out.fig, paste0('ResSum_', comp, '.png')),
        width = 1500, height = 1200, res = 150)
    
    grid::grid.newpage()
    if(length(ls.miRNA2gene) != 1) {
        if(length(ls.GOsum) != 1) {
            vert = spreadVertical(allmiRNA = allmiRNA,
                                  noNAmiRNA = noNAmiRNA,
                                  DEmiRNA = DEmiRNA,
                                  gene = gene,
                                  GO = GO)
        } else {
            vert = spreadVertical(allmiRNA = allmiRNA,
                                  noNAmiRNA = noNAmiRNA,
                                  DEmiRNA = DEmiRNA,
                                  gene = gene)
        }
    } else {
        vert = spreadVertical(allmiRNA = allmiRNA,
                              noNAmiRNA = noNAmiRNA,
                              DEmiRNA = DEmiRNA)
    }
    NAmiRNA = moveBox(NAmiRNA, 
                      x = 0.8,
                      y = coords(vert$noNAmiRNA)$top + distance(vert$allmiRNA, 
                                                                vert$noNAmiRNA,
                                                                half = TRUE,
                                                                center = FALSE))
    noDEmiRNA = moveBox(noDEmiRNA,
                        x = 0.8,
                        y = coords(vert$DEmiRNA)$top + distance(vert$noNAmiRNA,
                                                                vert$DEmiRNA,
                                                                half = TRUE,
                                                                center = FALSE))
    for (i in 1:(length(vert) - 1)) {
        connectGrob(vert[[i]], vert[[i + 1]], type = 'vert') %>% print()
    }
    
    print(connectGrob(vert$allmiRNA, NAmiRNA, type = 'L'))
    print(connectGrob(vert$noNAmiRNA, noDEmiRNA, type = 'L' ))
    
    print(vert)
    print(NAmiRNA)
    print(noDEmiRNA)
    
    dev.off()

}




