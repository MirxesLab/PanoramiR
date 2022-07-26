# =========================================================================== #
# Step 1: Read in Sample sheet
# --------------------------------------------------------------------------- #
# Input: file.samplesheet

# Output:
    # num.sample
    # comparisons ['Comparison 1', 'Comparison 2', 'Comparison 0'] [character]
    
    # ls.compare.group ['Comparison 1', 'Comparison 2']
        # For data frame named 'Comparison 1' ['Samples', 'Comparison 1']
    # df.sampleID ['Samples', 'SampleID'] [data frame]
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

# Remove rows that do not have 'Unique Sample ID' 
    # Meaning that these samples cannot use miRNA assay panel (due to low concentration of RNA or serious hemolysis)
df.samplesheet = df.samplesheet[!is.na(df.samplesheet[, 5]), ]

colnames(df.samplesheet)[1] <- "Samples"

## Check sample name
if (!is.numeric(df.samplesheet$Samples)) {
    stop("Sorry, sample name should be number")
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

# Remove: 'comparisons' has the same function
# title.compare <- sapply(comparisons,
#                         function(x)
#                             ifelse(x=="Comparison 0",
#                                    NA,
#                                    x))

ls.compare.group <- list()

for (i in 1:sum(n.comp.group)) {
    df.tmp = df.samplesheet[,c(1,(1+i), 7)] 
    df.tmp = df.tmp %>%
        dplyr::filter(!is.na(df.tmp[2]))
    
    # Samples with the same sample type can be compared
    if(length(unlist(unique(df.tmp[,3]))) > 1) {
        stop('Sorry, only samples with the same sample type can be comapred')
    } else {
        df.tmp = df.tmp[, c(1,2)]
    }
    ls.compare.group[[i]] <- df.tmp    
}
names(ls.compare.group) <- comparisons[1:sum(n.comp.group)]


# Data frame for reading in raw data, sampleID : 'projectID sample'
df.sampleID = df.samplesheet[, c(1, 5)]
colnames(df.sampleID) = c('Samples', 'SampleID')

num.sample = nrow(df.samplesheet)
splitID = str_split(df.sampleID$SampleID, ' ')

# For Biofluid, 2 samples in 1 result file, merge sampleID to get the name of result file 
fun.mergeID = function(even) {
    # even: a even number. 
        # num.sample = 2n, even = num.sample
        # num.sample = 2n + 1, even = num.sample - 1
    
    # row in the final data frame
    fun.mergeRow = function(i) {
        row = data.frame(Samples = paste0(i, '_', i + 1),
                         SampleID = paste0(splitID[[1]][1],
                                           ' ',
                                           splitID[[i]][2],
                                           '_',
                                           splitID[[i + 1]][2]))
        return(row)
    }
    
    # Make the sampleID data frame
    for(i in seq.int(1, even, 2)) {
        if (i == 1) {
            df.tmp = fun.mergeRow(i)
        } else {
            row.tmp = fun.mergeRow(i)
            df.tmp = rbind(df.tmp, row.tmp)
        }
    }
    return(df.tmp)
}

if(arg.pipeline == 'Biofluid') {
    if (num.sample %% 2 == 0) {
        df.sampleID = fun.mergeID(num.sample)
    } else {
        df.sampleID.tmp = fun.mergeID(num.sample - 1)
        df.sampleID = rbind(df.sampleID.tmp, df.sampleID[num.sample, ])
    }
}


