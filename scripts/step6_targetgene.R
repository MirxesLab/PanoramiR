# ============================================================================= #
# Input:
    # ls.diff
    # file.ref.target
# Output:
    # fun.targetGene(ls.diff) [Find the target gene]
        # ls.miRNA2gene: miRNA2gene.sig, symbol.pvalue, miRNA2gene.show
# ============================================================================= #

# Read in Reference file
ref.target = read.csv(file.ref.target, sep = '\t')
colnames(ref.target) = c('miRNA', 'GeneSymbol', 'GeneID', 'EnsemblID')

# Function: find target gene
fun.targetGene = function(ls.diff) {
    # ls.diff is the list recording the results of DE analysis
    
    ls.miRNA2gene = list()
    
    df.res.ori = ls.diff$res.ori         # Contain all miRNA
    df.res.filter = ls.diff$res.filter   # Contain miRNA with pvalue <= 0.05, abs(dCt) >= 1
    
    if(nrow(df.res.filter) == 0) {      
        ls.miRNA2gene = NA               # No miRNA meets the requirement, length(ls.miRNA2gene) = 1
    } else {
        miRNA2gene.all = ref.target %>%
            dplyr::inner_join(df.res.ori, by = 'miRNA') %>%
            dplyr::select(miRNA, GeneSymbol, GeneID, EnsemblID) %>%
            unique()
        
        miRNA2gene.sig = ref.target %>%
            dplyr::inner_join(df.res.filter, by = 'miRNA') %>%
            dplyr::select(miRNA, GeneSymbol, GeneID, EnsemblID) %>%
            unique()
    }
    
    if (nrow(miRNA2gene.sig) == 0) {
        ls.miRNA2gene = NA             # No target gene is found, length(ls.miRNA2gene) = 1
    } else {
        # generate gene list for go analysis
        # the artificial pvalue of gene targeted by selected miRNA is 0.001
        # the artificial pvalue of other genes is 1
        miRNA2gene.all = miRNA2gene.all %>%
            dplyr::mutate(art.pvalue = ifelse(miRNA2gene.all$GeneSymbol %in% miRNA2gene.sig$GeneSymbol,
                                              0.001,
                                              1))
        
        # Gene Symbol with artificial p value for GO analysis
        symbol.pvalue = data.frame(symbol = miRNA2gene.all$GeneSymbol,
                                   art.pvalue = miRNA2gene.all$art.pvalue) %>%
            unique()
        
        # Generate data frame for target gene analysis
        miRNA.select = unique(miRNA2gene.sig$miRNA)
        df.miRNA2gene.show = data.frame() 
        for (i in miRNA.select) {
            # Selected 5 target genes randomly
            index = sample(which(miRNA2gene.sig$miRNA == i),
                           size = 5)
            df.miRNA2gene.show = rbind(df.miRNA2gene.show,
                                       miRNA2gene.sig[index, ])
        }
        rownames(df.miRNA2gene.show) = NULL
        
        ls.miRNA2gene$miRNA2gene.sig = miRNA2gene.sig # Save
        ls.miRNA2gene$symbol.pvalue  = symbol.pvalue  # For GO enrichment
        ls.miRNA2gene$miRNA2gene.show = df.miRNA2gene.show # For Report
    }
    
    return(ls.miRNA2gene)
}