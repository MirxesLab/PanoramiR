# =========================================================================== #
# Step 1.1: Preprocess Sample sheet
# --------------------------------------------------------------------------- #
# input:
    # file.samplesheet
# output:
    # df.samplesheet
# =========================================================================== #

# Read in sample sheet
df.samplesheet.ori = read_excel(file.samplesheet, skip = 33)

# Check Sample Sheet
keep.idx <- which(!apply(apply(df.samplesheet.ori[,2:7],
                               2,
                               is.na),
                         1,
                         all))
df.samplesheet = df.samplesheet.ori[keep.idx,]
colnames(df.samplesheet)[1] <- "Samples"

## Check sample name
if (!is.numeric(df.samplesheet$Samples)) {
    stop("Sample name should be number")
}


## Check the comparison 1 ~ 3
n.comp.group <- apply(df.samplesheet[,2:4],
                      2,
                      function(x)
                          !all(is.na(x)))
comparisons <- sapply(1:3,
                      function(x)
                          ifelse(n.comp.group[x],
                                 names(n.comp.group)[x],
                                 "Comparison 0"))
title.compare <- sapply(comparisons,
                        function(x)
                            ifelse(x=="Comparison 0",
                                   NA,
                                   x))

ls.compare.group <- list()

for (i in 1:sum(n.comp.group)) {
    ls.compare.group[[i]] <- df.samplesheet[,c(1,(1+i))]    
}

names(ls.compare.group) <- title.compare[1:sum(n.comp.group)]
