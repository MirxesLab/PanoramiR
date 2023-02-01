# =========================================================================== #
# Step 5 Normalization and Round 2 Filter
# --------------------------------------------------------------------------- #
# input
    # df.input.data.Filt1
    # threshold.sp
    # file.RNAvol
    # cutoff.sp
    

# Output:
    # df.input.data.IPC.norm [data normalized by IPC]      -> Save
    # SP.sum            [summary of Spike-in CT value]     -> Save
    # df.input.data.sp  [data filtered based on Spike-in threshold]
    # df.input.data.SP.norm  [data normalized by SP]       -> Save
    # df.ipc.factor    
    # df.sp.factor

# =========================================================================== #
group = unique(df.input.data.Filt1$Group) %>% sort()
# Cancer: A, B, C, D [4RT group -> 4 * 4 SP/IPC]
# Biofluid: A, B     [2RT group -> 4 * 2 SP/IPC]
# PanoramiR: A1, B1, C1, D1, A2, B2, C2, D2 [8RT group -> 1 * 8 SP/IPC]

fun.norm = function(df, df.factor) {
    # df is the data which need to be normalized [input data]
    # df.factor is the IPC factor or SP factor
    # mode = 'SP' or 'IPC'
    df.norm = df
    if(arg.pipeline == 'PanoramiR') {
        df.factor = df.factor[, -2]
    }
    num.sample = ncol(df) - 2
    for(i in group) {
        input.tmp = df.norm[df.norm$Group == i, ]
        factor.tmp = df.factor[df.factor$Group == i, ]
        
        for (j in 1:num.sample) {
            input.tmp[, j + 2] = input.tmp[, j + 2] - factor.tmp[, j+1]
        }
        df.norm[df.norm$Group == i, ] = input.tmp
    }
    return(df.norm)
}
fun.df.factor = function(df) {
    # df is the data frame of IPC or SP
    
    # calculate factors based on RT group
    fun.factor = function(df.x) {
        # df.x is the sub data frame of df 
        # special situationL when both IPC1 or IPC2 is NA
        mean.IPC1 = apply(df.x[1:2, -c(1,2)], 2, mean)
        mean.IPC2 = apply(df.x[3:4, -c(1,2)], 2, mean)
        
        mean.IPC1.all = mean(mean.IPC1, na.rm = TRUE)
        mean.IPC2.all = mean(mean.IPC2, na.rm = TRUE)
        
        mean.IPC1 = ifelse(is.na(mean.IPC1), mean.IPC1.all, mean.IPC1)
        mean.IPC2 = ifelse(is.na(mean.IPC2), mean.IPC2.all, mean.IPC2)
        
        df.tmp = data.frame(factorIPC1 = mean.IPC1 - mean.IPC1.all,
                            factorIPC2 = mean.IPC2 - mean.IPC2.all)
        
        factor = apply(df.tmp, 1, mean)
        return(factor)
    }
    
    df.factor = data.frame()
    for (i in group) {
        df.tmp = df[df$Group == i, ]
        df.factor.tmp = fun.factor(df.tmp)
        if (i == 'A') {
            df.factor = df.factor.tmp
        } else {
            df.factor = rbind(df.factor, df.factor.tmp)
        }
    }
    tmp.colname = colnames(df.factor)
    
    df.factor = as.data.frame(df.factor) %>%
        dplyr::mutate(Group = group) %>%
        dplyr::select('Group', all_of(tmp.colname))
    rownames(df.factor) = NULL
    
    return(df.factor)
} # Function for calculate normlization factor [SP / IPC]
                    # --------------------------------- #
                    #       IPC Normalization           #
                    # --------------------------------- #
# input: df.input.data.Filter1
# output: df.input.data.IPC.norm
if(arg.pipeline != 'PanoramiR') {
    index.ipc = which(df.input.data.Filt1$miRNA %in% name.ipc)
    
    IPC = df.input.data.Filt1[index.ipc, ]
    df.ipc.factor = fun.df.factor(IPC)
    
    # IPC Normalization
    df.input.data.IPC.norm = df.input.data.Filt1[-index.ipc, ]
    df.input.data.IPC.norm = fun.norm(df.input.data.IPC.norm, df.ipc.factor)
    
    # Save Results
    tmp = as.data.frame(df.input.data.IPC.norm)
    # colnames(tmp)[-c(1,2)] = paste('sample', colnames(tmp)[-c(1,2)])
    
    write.csv(tmp, file.path(dir.out.tbl, 'Data IPCnorm.csv'))
    
} else {df.input.data.IPC.norm = df.input.data.Filt1}



                    # --------------------------------- #
                    #      Check SP (only for RT SP)    #
                    # --------------------------------- #
# input: df.input.data.IPC.norm
# output: df.input.data.sp, SP
# name.sp for PanoramiR is [A1.SP, B1.SP, C1.SP, D1.SP ... D2.SP]
if(arg.pipeline == 'PanoramiR') {
    name.sp = grep('*.SP', 
                   df.input.data.Filt1$miRNA,
                   value = TRUE)}

# Isolation Spike-ins do not need this step
stop = FALSE
if(is.RTsp) {
    index.sp = which(df.input.data.IPC.norm$miRNA %in% name.sp)
    SP = df.input.data.IPC.norm[index.sp, ]
    
    # This function is used for Cancer / Biofluid results
    fun.diff.sp  = function(groupID) {
        # group ID is element in group [A, B, C, D] or [A, B]
        sp = SP[SP$Group == groupID, ]
        sp.1 = sp[sp$miRNA == name.sp[1], -c(1,2)]
        sp.2 = sp[sp$miRNA == name.sp[2], -c(1,2)]
        
        # Plan A
        sd.1 = apply(sp.1, 2, sd, na.rm = TRUE)
        sd.2 = apply(sp.2, 2, sd, na.rm = TRUE)

        df.sd = cbind(sd.1, sd.2) %>% round(3)
        colnames(df.sd) = paste0(paste0('sd.', groupID), c('.rep1', '.rep2'))
        
        # Plan B
        # mean.1 = apply(sp.1, 2, mean, na.rm = TRUE)
        # mean.2 = apply(sp.2, 2, mean, na.rm = TRUE)
        # diff   = mean.1 - mean.2
        # 
        # df.diff = data.frame(diff = diff) %>% round(3)
        # colnames(df.diff) = paste0('diff.', groupID)
        
        return(df.sd)
    }  # Calculate the difference between Spike-ins
    
    # Create the SP summary data frame: [sample, statue, max SP, difference]
    SP.sum = data.frame()
    if (arg.pipeline == 'PanoramiR') {
        mean = apply(SP[, -c(1,2)], 2, mean, na.rm = TRUE)
        sd   = apply(SP[, -c(1,2)], 2, sd, na.rm = TRUE)
        SP.sum = data.frame(mean = mean,
                            sd   = sd)
    } else {
        for (i in group) {
            if (i == 'A') {
                SP.sum = fun.diff.sp(i)
            } else {
                tmp = fun.diff.sp(i)
                SP.sum = cbind(SP.sum, tmp)
            }
        }
    } 
    
    tmp.colname = colnames(SP.sum)
    max = apply(SP[, -c(1,2)], 2, max, na.rm = TRUE)
    SP.sum = SP.sum %>% as.data.frame() %>%
        dplyr::mutate(max.SP = max,
                      status = ifelse(max > threshold.sp,
                                      'exclude',
                                      'include')) %>%
        dplyr::mutate(sample = rownames(SP.sum)) %>%
        dplyr::mutate(threshold = threshold.sp) %>%
        dplyr::select('sample', 'max.SP', 'threshold', 'status', all_of(tmp.colname))
    
    # Save results
    write.csv(SP.sum, file.path(dir.out.tbl, 'SP_Check_Summary.csv'))

    # Remove samples based on SP.sum
    remove.sample = SP.sum$sample[SP.sum$status == 'exclude'] %>% as.numeric()
    num.remove    = length(remove.sample)
    num.sample    = ncol(df.input.data.IPC.norm) - 2
    
    # Check the remaining samples
        # cond 1: num.remove == num.sample
        # cond 2: num.remove > num.sample - 4
        # cond 3: 0 < num.remove <= num.sample - 4
        # cond 4: num.remove = 0
    
    if(num.remove == num.sample) {
        stop = TRUE
        description = paste0('Sorry, all samples are removed because all samples contain at least one Spike-In having Ct value greater than ',
                   threshold.sp, 
                   '. Please check the raw data or change the threshold of Spike-in')
    } else if(num.sample - num.remove < 4) {
        stop = TRUE
        description = 'Sorry, after verification of Spike-in, the remaining samples cannot be used for differential expression analysis as the number of remaining samples is less than 4'
    } else if(num.remove != 0) {
        # make sure number of samples in every group in comparisons is bigger than 2
        comp.tmp = c('Comparison 1', 'Comparison 2', 'Comparison 3')
        comp.remove = c()
        
        for (i in comp.tmp) {
            if (i %in% comparisons) {
                comp = ls.compare.group[[i]]
                index = which(comp$Samples  %in% remove.sample)
                tmp  = comp[-index, ]
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
            description = 'Sorry, after verification of Spike-in, the remaining samples cannot be used for differential expression analysis. All defined comparisons contain a group in which the number of samples is less than 2'
        } else {
            stop = FALSE
            description = ''
            df.input.data.sp = df.input.data.IPC.norm[, -(remove.sample + 2)]
        }
    } else {
        stop = FALSE
        description = ''
        df.input.data.sp = df.input.data.IPC.norm
    }
    
} else {
    stop = FALSE
    description = ''
    df.input.data.sp = df.input.data.IPC.norm
}


                    # --------------------------------- #
                    #          SP Normalization         #
                    # --------------------------------- #
# input: df.input.data.sp, SP
# output: df.input.data.SP.norm
fun.sp.factor = function(sp) {
    # sp is the data frame of spkie-in Ct values [normalized by RNA input volume or not]
    # sp: [Group, miRNA, 1, 2, ....]
    if(arg.pipeline == 'PanoramiR') {
        factor = sp
        factor[, -c(1,2)] = apply(factor[, -c(1,2)], 
                                  1, 
                                  function(x) x - mean(x, na.rm = TRUE)) %>% t()
        factor = as.data.frame(factor)
    } else {
        factor = fun.df.factor(sp) # function in IPC normalization
    }
    return(factor)
} # two modes: PanoramiR / Cancer & Biofluid

# Calculate Spike-in factor: df.sp.factor
index.sp = which(df.input.data.IPC.norm$miRNA %in% name.sp)
SP = df.input.data.IPC.norm[index.sp, ]
df.sp.factor = fun.sp.factor(SP) # RNA volume is not taken into consideration
# ----- History ------
#if(is.RTsp) {
#    # Calculate the Spike-in Factor
#    df.sp.factor = fun.sp.factor(SP)
    
#} else {
#    # Read in RNA input volume and check the file: all samples should exist in RNA input volume
#     RNAvol = read_excel(file.RNAvol)
#     sample = colnames(df.input.data.sp)[-c(1,2)]
#     index  = which(!(sample %in% RNAvol$sample))
#     if(length(index > 0)) {
#         stop(paste('Sorry, cannot find sample', 
#                    sample[index], 
#                    'in RNA volume excel file. Please Check RNA volume excel file'))
#     } else {
#         RNAvol = RNAvol[RNAvol$sample %in% sample, ]
#     }
#     
#     # Calculate Input Scaling Factor
#     RNAvol = RNAvol %>% dplyr::mutate(log2vol = log2(RNAvol$volume))
#     mean.log2vol = mean(RNAvol$log2vol)
#     RNAvol = RNAvol %>%
#         dplyr::mutate(InputSF = RNAvol$log2vol - mean.log2vol)
#     
#     # Calculate RNA volume normalized SP
#     SP.volNorm = SP
#     num.sample = ncol(df.input.data.sp) - 2
#     for (i in 1:num.sample) {
#         SP.volNorm[, i+2] = SP.volNorm[, i+2] - RNAvol$InputSF[i]
#     }
#     
#     # Calculate Spike-in factor: Cancer&Biofluid / PanoramiR
#     df.sp.factor = fun.sp.factor(SP.volNorm)
# } # two modes: normalized by RNA input volume / no normalization

# Normalization -----
index.sp = which(df.input.data.sp$miRNA %in% name.sp)
df.input.data.SP.norm = df.input.data.sp[-index.sp, ]
df.input.data.SP.norm = fun.norm(df.input.data.SP.norm, df.sp.factor)



    