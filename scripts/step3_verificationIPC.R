# =========================================================================== #
# Step 3 Verification of IPC [Only for Biofuild & Cancer Panel]
# --------------------------------------------------------------------------- #
# input
    # df.input.data
    # ls.compare.group
    # num.sample

    # name.ipc [Cancer / Biofluid]
    # check.ipc.strict
    

# output:
    # failure.ipc [record the removed sample IDs [1,2,...]]
    # ipc.report [Parameter used in report]
    # IPC.sum 
    # df.input.data.ipc

# =========================================================================== #
# Prepare
index = which(df.input.data$miRNA %in% name.ipc)
IPC = df.input.data[index, ]
IPC[, -c(1,2)] = round(IPC[, -c(1,2)], 3)
group = unique(IPC$Group) # [A, B, C, D] or [A, B]
group = sort(group)

# samples having one IPC Ct > 25 are removed
failure.ipc = which(apply(IPC[, -c(1,2)], 
                          2,
                          function(x) sum(x > threshold.ipc)) > 0)

# IPC summary, save results automatically
for(i in group) {
    ipc = IPC[IPC$Group == i, ]
    ipc.1 = ipc[ipc$miRNA == name.ipc[1], ]
    ipc.2 = ipc[ipc$miRNA == name.ipc[2], ]
    
    sd.1   = round(apply(ipc.1[, -c(1:2)], 2, sd, na.rm = TRUE), 3)
    sd.2   = round(apply(ipc.2[, -c(1:2)], 2, sd, na.rm = TRUE), 3)
    
    if (i == 'A') {
        IPC.sum = data.frame(sd.A.rep1 = sd.1,
                             sd.A.rep2 = sd.2)
    } else {
        df.tmp = data.frame(sd.1, sd.2)
        colnames(df.tmp) = paste0('sd.', i, '.rep',c(1:2))
        IPC.sum = cbind(IPC.sum, df.tmp)
    }
}
ipc.max = apply(IPC[, -c(1,2)], 2, max, na.rm = TRUE)
tmp.colname = colnames(IPC.sum)
IPC.sum = IPC.sum %>%
    dplyr::mutate(max.IPC = ipc.max) %>%
    dplyr::mutate(status = ifelse(ipc.max > threshold.ipc,
                                  'exclude',
                                  'include')) %>%
    dplyr::mutate(sample = rownames(IPC.sum)) %>%
    dplyr::mutate(threshold = threshold.ipc) %>%
    dplyr::select('sample', 'max.IPC', 'threshold', 'status', all_of(tmp.colname))

write.csv(IPC.sum, 
           file.path(dir.out.tbl, 'IPC_Check_Summary.csv'))


# Whether later analysis is applicable, Error conditions below
    # Cond 1: no smaple remains
    # Cond 2: number of remaining samples < 4 [one group in all comparisons have less than 2 samples]
    # Cond 3: All comparisons contain groups with number of remaining samples < 2
stop = FALSE
if(length(failure.ipc) == num.sample) {
    stop = TRUE
    description = paste0('Sorry, all samples contain IPC having Ct value greater than ',
                check.ipc,
                '. Please check raw data or change the threshold of IPC')
    
} else if(num.sample - length(failure.ipc) < 4) {
    stop = TRUE
    description = 'Sorry, after verification of IPC, the remaining samples cannot be used for differential expression analysis as the number of remaining sample is less than 4'
    
    
} else if (length(failure.ipc) != 0) {
    # num.sample - failure.ipc > 4 && failure.ipc != 0
    df.input.data.ipc = df.input.data[, -(failure.ipc + 2)]
    # make sure number of samples in every group in comparisons is bigger than 2
    comp.tmp = c('Comparison 1', 'Comparison 2', 'Comparison 3')
    comp.remove = c()
    
    for (i in comp.tmp) {
        if (i %in% comparisons) {
            comp = ls.compare.group[[i]]
            tmp  = comp[-failure.ipc, ]
            num.A = sum(tmp[[i]] == 'A')
            num.B = sum(tmp[[i]] == 'B')
            if (num.A < 2) { 
                comp.remove = append(comp.remove, i)
            } else if (num.B < 2) {
                comp.remove = append(comp.remove, i)
            }
        }
    }
    idx = which(comparisons %in% comp.remove)
    comparisons[idx] = 'Comparison 0'
    if(sum(comparisons == 'Comparison 0') == 3) {
        stop = TRUE
        description = 'Sorry, after verification of IPC, the remaining samples cannot be used for differential expression analysis. All defined comparisons contain a group of which the number of samples is less than 2'
    } else {
        stop = FALSE
        description = ''
    }
} else {
    df.input.data.ipc = df.input.data
    stop = FALSE
    description = ''
}



    
    

    


