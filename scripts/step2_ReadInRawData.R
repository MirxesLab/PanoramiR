# =========================================================================== #
# Step 2 Read in raw Data
# --------------------------------------------------------------------------- #
# input
    # df.sampleID [step1_ReadInSamplesheet: 'Samples', 'SampleID']
    # file.ref.miRList [group, position, miRNA ID]
    # result.skip = 46 [by default]
    # num.sample

# output:
    # df.input.data [colnames: 'Group', 'miRNA', '1', '2', ...  or '1_2', '3_4', ...]
    # num.file [only PanoramiR and Cancer: num.file = num.sample]
# =========================================================================== #

# Relationship between well position and miRNA ID
# --------------------------------------------------------------------------- #
# Three important columns ['Group', 'Position', 'miRNA ID']
miRList = read_excel(file.ref.miRList) %>%  
    dplyr::select('Group', 'Position', 'miRNA ID')

colnames(miRList) = c('Group','Position', 'miRNA') 
# Read in Results files
# --------------------------------------------------------------------------- #
# Read in files based on data frame: df.sampleID
num.file = nrow(df.sampleID)
# for PanoramiR and Cancer : num.file == num.sample
# for Biofluid: num.file == num.sample %/% 2 + 1

# Define function to read in one result file
fun.readRes = function(i) {
    # i is the i_th row of df.sampleID
    
    # Record sample and unique sample id
    no.sample = df.sampleID$Samples[i]
    id.sample = df.sampleID$SampleID[i]
    pat       = paste0(id.sample, '_')
    file      = grep(pat, list.files(dir.input, all.files = FALSE, pattern = '.xlsx'), value = TRUE)
    
    # Read in raw data
    if (length(file) != 0 ){
        res.ori = read_excel(file.path(dir.input, file), 
                             'Results', 
                             skip = result.skip)
        res.ext = res.ori[, c('Well Position', 'CT')]
        colnames(res.ext) = c('Position', no.sample)
        str_sub(res.ext$Position[which(str_length(res.ext$Position) == 2)], 
                2, 
                1) = 0 # A1 -> A01
    } else {
        res.ext = data.frame(Position = miRList$Position,
                             CT = NA)
        colnames(res.ext) = c('Position', no.sample)
    }
    
    res.ext = res.ext[1:nrow(miRList), ]
    return(res.ext)
}

# Read in files and generate df.input.data
### Note: row is miRNAs, columns: [ 'Position', df.sampleID$Samples ]
for(i in 1:num.file) {
    if ( i == 1 ) {
        df.input.data = fun.readRes(i) 
    } else {
        df.tmp = fun.readRes(i)
        df.input.data = df.input.data %>%
            dplyr::inner_join(df.tmp, by = 'Position')
    }
}
# Some elements in data.frame is 'undetermined', these elements are converted into NA
df.input.data[, -1] = apply(df.input.data[, -1], 2, as.numeric)

# Link CT value with miRNA ID by the 'Position'
df.input.data = df.input.data %>%
    dplyr::inner_join(miRList, by = 'Position') %>%
    dplyr::select('Group', 'miRNA', all_of(colnames(df.input.data)[-1]))

# Split input data [Only for Biofluid]
if(arg.pipeline == 'Biofluid') {
    # Create New input.data only having 192 rows
    fun.newInput = function(filenumber) {
        # if num.sample is an even number, filenumber = num.file
        # if num.sample is an even number, filenumber = num.file - 1
        
        df.input.data.tmp = data.frame()
        
        fun.splitInput = function(df) {
            # df is df.tmp
            df.tmp.1 = df[df$Group %in% c('A1', 'B1'), ]
            df.tmp.2 = df[df$Group %in% c('A2', 'B2'), ]
            
            df.tmp = data.frame(df.tmp.1[, 3], df.tmp.2[, 3])
            colnames(df.tmp) = no.sample
            
            return(df.tmp)
        }
        
        for (i in 1:filenumber) {
            sample = df.sampleID$Samples[i]
            no.sample = unlist(str_split(sample, '_'))
            df.tmp = df.input.data[, c('Group', 'miRNA', sample)]
            df.tmp = fun.splitInput(df.tmp)
            if(i == 1) {
                df.input.data.tmp = df.tmp
            } else {
                df.input.data.tmp = cbind(df.input.data.tmp, df.tmp)
            }
        }
        df.tmp = df.input.data[df.input.data$Group %in% c('A1', 'B1'), c(1,2)]
        df.tmp$Group = factor(df.tmp$Group)
        levels(df.tmp$Group) = c('A', 'B')
        
        df.input.data.tmp = cbind(df.tmp, df.input.data.tmp)
        
        return(df.input.data.tmp)
    }
    
    if (num.sample %% 2 == 0) {
        df.input.data = fun.newInput(num.file)
    } else {
        # For the last column, Group.A1 and Group.B1 is the data for the last sample
        df.input.data.tmp = fun.newInput(num.file - 1) # ncol is an even number
        # Extract results of the last sample
        df.tmp = df.input.data[df.input.data$Group %in% c('A1', 'B1'),
                               num.file + 2]
        df.input.data = cbind(df.input.data.tmp, df.tmp)
    }
}


# ---------------------------------------------------------------------------- #
# Save Results
write.csv(df.input.data, file.path(dir.out.tbl, 'Data Raw.csv'))
