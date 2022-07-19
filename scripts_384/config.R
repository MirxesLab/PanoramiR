# =========================================================================== #
# Configure file for panoramiR.Rmd
# =========================================================================== #

# Load Packages
library(dplyr)
library(readxl)
library(stringi)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(knitr)
library(kableExtra)
library(topGO)
library(Gmisc)
library(glue)
library(grid)

# Set Path
dir.top     = '/Users/plateau/Documents/GitHub/PanoramiR'
dir.out.fig = 'output_384/figures'
dir.out.tbl = 'output_384/tables'
dir.input   = 'input_384'
dir.resource = 'resources'


# Set working directory
setwd(dir.top)

# Create directory for output
dir.create(dir.out.fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.out.tbl, recursive = TRUE, showWarnings = FALSE)

# Set file path

file.input.data  = file.path(dir.input, 
                             grep('*_Results.xlsx', 
                                  list.files(dir.input, all.files = FALSE), 
                                  value = TRUE))    #[*.Results.xlsx]
file.samplesheet = 'input_384/samplesheet.xlsx'
png.workflow     = file.path(dir.resource, 'knowledge_384.png')
file.ref.miRList = file.path(dir.resource, 'miRList_384.xlsx')        # All miRNA: Group, Position, miRBASE v22 Accession, miRNAID
file.ref.RNAvol  = file.path(dir.resource, 'RNAvolume.xlsx')
file.ref.target  = file.path(dir.resource, 'Predicted_Targets_Human_Filtered_Rearranged_GeneID.txt') # Will be changed into miTarBase data
name.sp          = c('Spike-in RNA Ctr 1', 'Spike-in RNA Ctr 2')
name.ipc         = c('Inter-plate Calibrator 1', 'Inter-plate Calibrator 2')

# Set parameters
is.basic       = TRUE        # Basic worflow execute Global Normalization. Non-basic workflow also execute normalization based on stable miRNA
is.RTsp        = FALSE       # If the spike-in is Reverse Transcript Spike-in

result.skip    = 46          # The number depends on the machine
cutoff.max     = 33          # The first round filter, before normalization
cutoff.min     = 9
cutoff.sp      = 32          # The second round filter, for spike-in normalization

threshold.impute    = 0.1    # No more than 10% missing value in miRNA
threshold.DE.pvalue = 0.05   # The p value of T test
threshold.DE.dCt    = 1      # The log2 fold change of geometric mean
threshold.GO.pvalue = 0.01
GO.nodesize         = 10


# Set Color
col.compare.1 = brewer.pal(6, "Dark2")[1:2] # Green - Orange
col.compare.2 = brewer.pal(6, "Dark2")[3:4] # Purple - Pink
col.compare.3 = brewer.pal(6, "Dark2")[5:6] # Green - Yellow
col.others = brewer.pal(11,'RdYlBu')[c(1, 11)] # Red - Blue
col.grey   = brewer.pal(8, "Set2")[8] # grey
col.heatmap = brewer.pal(11, 'PRGn')[c(1, 11)] # Purple - Green
col.heatmap = c(col.heatmap[1], "#FFFFFF", col.heatmap[2])



