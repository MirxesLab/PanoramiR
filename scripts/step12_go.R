# ============================================================================= #
# Input:
    # ls.miRNA2gene [symbol.pvalue]
    # threshold.GO.pvalue [by default 0.01]
    # GO.nodesize         [by default 10]
# Output:
    # fun.GOenrich(ls.miRNA2gene, ont) [used in fun.GOsum]
        # ls.GOres [res.save; res.show] or [NA] 
    # fun.GOsum (ls.miRNA2gene)
        # ls.GOsum [res.save; res.show] or [NA]
# ============================================================================= #

fun.GOenrich = function(ls.miRNA2gene, ont) {
    # ls.miRNA2gene is the list saving the results of target gene analysis
    # ont is the GO ontology, 'CC', 'MF', 'BP'
    
    ls.GOres = list()
    
    # Create gene list
    genelist = ls.miRNA2gene$symbol.pvalue$art.pvalue
    names(genelist) = ls.miRNA2gene$symbol.pvalue$symbol
    
    # Create gene select function
    topDiffGenes = function(genelist) {return(genelist < 0.01)}
    
    # Create topGOdata object
    GOdata = new('topGOdata',
                 ontology = ont,
                 allGenes = genelist,
                 geneSel = topDiffGenes,
                 nodeSize = GO.nodesize,
                 annot = annFUN.org, mapping = 'org.Hs.eg.db', ID = 'symbol')
    
    # Enrichment test
    classicfisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
    pvalue = score(classicfisher)
    
    # Data frame in report (part of the GO terms)
    num.node = length(pvalue[pvalue < threshold.GO.pvalue])
    
    if (num.node == 0) {
        ls.GOres = NA       # based on the threshold, no GO term is found
    } else if (num.node < 10) {
        res.show = GenTable(GOdata,
                            classicFisher = classicfisher,
                            topNodes = num.node)
        res.show$`GO Term` = rep(ont, num.node)
        colnames(res.show)[6] = 'pvalue'
        
    } else {
        res.show = GenTable(GOdata,
                            classicFisher = classicfisher,
                            topNodes = 10)
        res.show$`GO Term` = rep(ont, 10)
        colnames(res.show)[6] = 'pvalue'
    }
    
    # Data frame for output (all GO terms)
    if (num.node == 0) {
        ls.GOres = NA 
    } else {
        res.save = GenTable(GOdata,
                            pvalue = classicfisher,
                            topNodes = num.node)
        res.save$`GO Term` = rep(ont, num.node)
        colnames(res.save)[6] = 'pvalue'
        
        
        ls.GOres$res.save = res.save
        ls.GOres$res.show = res.show
    }
    
    return(ls.GOres)
}

fun.GOsum = function(ls.miRNA2gene) {
    
    ls.GOsum = list()
    
    # GO enrichment
    goRes.bp = fun.GOenrich(ls.miRNA2gene, 'BP')
    goRes.cc = fun.GOenrich(ls.miRNA2gene, 'CC')
    goRes.mf = fun.GOenrich(ls.miRNA2gene, 'MF')
    
    df.res.show = data.frame()
    df.res.save = data.frame()
    
    # Summary the results (p value < 0.01)
    if(length(goRes.bp) != 1) {
        df.res.show = rbind(df.res.show, goRes.bp$res.show)
        df.res.save = rbind(df.res.save, goRes.bp$res.save)
    }
    if(length(goRes.cc) != 1) {
        df.res.show = rbind(df.res.show, goRes.cc$res.show)
        df.res.save = rbind(df.res.save, goRes.cc$res.save)
    }
    if(length(goRes.mf) != 1) {
        df.res.show = rbind(df.res.show, goRes.mf$res.show)
        df.res.save = rbind(df.res.save, goRes.mf$res.save)
    }
    
    # nrow == 0 means no significantly enriched GO hits
    if(nrow(df.res.show) != 0) {
        ls.GOsum$res.show = df.res.show
        ls.GOsum$res.save = df.res.save
    } else {
        ls.GOsum = NA
    }

    return(ls.GOsum)
}

