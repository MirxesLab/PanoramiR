# =========================================================================== #
# Step 1: Pre-processing: Extract Ct values and Filter based on Recommend Cut-off
# --------------------------------------------------------------------------- #
# input
    # file.ref.miRList
    # file.input.data

# Params:
    # cutoff.max
    # cutoff.min

# output:
# 1. df.input.data: [miRNA, 1, 2, 3, ...]
# 2. df.input.data.Filt1 [miRNA, 1, 2, 3, ...]
# 3. miRList [Group, Position, miRNA]
# 4. num.sample

# =========================================================================== #

# Read in file recording relationship between well position and miRNA ID
miRList = read_excel(file.ref.miRList) %>%  # [Group, Position, miRBASE v22 Accession, miRNA ID]
    dplyr::select(Position, 'miRNA ID')
colnames(miRList) = c('Position', 'miRNA')  # [Position, miRNA]

# --------------------------------------------------------------------------- #
# Read in Results files
num.sample = length(file.input.data)

for (i in 1:num.sample) {
    file.tmp = file.input.data[i]
    
    # Record No.sample
    no.sample = as.numeric(unlist(strsplit(unlist(strsplit(file.tmp, ' '))[3], '_'))[1])
    
    # Read in result file
    res.ori = read_excel(file.tmp, skip = result.skip)
    res.ext = res.ori[, c('Well Position', 'CT')]
    colnames(res.ext) = c('Position', 'CT')
    str_sub(res.ext$Position[which(str_length(res.ext$Position) == 2)], 2, 1) = 0 # A1 -> A01
    
    # Link CT values with miRN ID
    res.ext = inner_join(res.ext, miRList, by = 'Position') %>%
        dplyr::select('miRNA', 'CT')
    colnames(res.ext) = c('miRNA', no.sample)
    
    # Merge CT values of Samples
    if(i == 1) {
        df.input.data = res.ext
    } else {
        res.ext = res.ext[res.ext$miRNA==df.input.data$miRNA, ]
        df.input.data = cbind(df.input.data, res.ext[, 2])
    }
}

# Change the order of columns
df.input.data = df.input.data[, c('miRNA', as.character(1:num.sample))]

# --------------------------------------------------------------------------- #
# First round filter: based on cutoff.max, cutoff.min
df.input.data.Filt1 = df.input.data

df.input.data.Filt1.tmp = apply(df.input.data[, -1], 
                                2,
                                function(x) ifelse(x < cutoff.min, NA, x)) %>%
    apply(2,
          function(x) ifelse(x > cutoff.max, NA, x))
df.input.data.Filt1[, -1] = df.input.data.Filt1.tmp
# --------------------------------------------------------------------------- #
# Save files
df.input.data.tmp = df.input.data
colnames(df.input.data.tmp)[-1] = paste0('sample_', colnames(df.input.data.tmp)[-1]) 
write.table(df.input.data.tmp, 
          file = file.path(dir.out.tbl, 'rawData.tsv'),
          sep = '\t',
          row.names = FALSE)


df.input.data.Filt1.tmp = df.input.data.Filt1
colnames(df.input.data.Filt1.tmp)[-1] = paste0('sample_', colnames(df.input.data.Filt1.tmp)[-1])
write.table(df.input.data.Filt1.tmp, 
            file = file.path(dir.out.tbl, 'filterData_round1.tsv'),
            sep = '\t',
            row.names = FALSE)
