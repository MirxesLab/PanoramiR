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
dir.top      = '/Users/plateau/Documents/GitHub/PanoramiR'
dir.input    = 'input_ncRNA'
dir.resource = 'resources'
dir.out.fig  = 'output/figures'
dir.out.tbl  = 'output/tables'

# Set working directory
setwd(dir.top)

# Create directory for output
dir.create(dir.out.fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.out.tbl, recursive = TRUE, showWarnings = FALSE)

# Set file path
file.input.data  = file.path(dir.input, 
                             grep('*_Results.xls', 
                                  list.files(dir.input, all.files = FALSE), 
                                  value = TRUE))    #[*.Results.xlsx]
file.samplesheet = file.path(dir.input, 'Sample Manifest Form.xlsx')
file.ref.miRList = file.path(dir.resource, 'miRList.xlsx')      # All miRNA: Group, Position, miRBASE v22 Accession, miRNAID
file.ref.miRsp   = file.path(dir.resource, 'miRSpike.xlsx' )       # Spike-in miRNA: Position, miRNA ID, SP, Group
file.ref.miRthld = file.path(dir.resource, 'miRThreshold.xlsx')    # Only for PanoramiR: miRNA, Threshold
file.ref.RNAvol  = file.path(dir.resource, 'RNAvolume.xlsx')
file.ref.target  = file.path(dir.resource, 'Predicted_Targets_Human_Filtered_Rearranged_GeneID.txt') # Will be changed into miTarBase data

png.workflow     = file.path(dir.resource, 'workflow.png')
png.mSMRT.qPCR   = file.path(dir.resource, 'mSMRT_qPCR.png')

# Set parameters
is.basic       = TRUE        # Basic worflow execute Global Normalization. Non-basic workflow also execute normalization based on stable miRNA
is.RTsp        = FALSE       # If the spike-in is Reverse Transcript Spike-in

result.skip    = 46          # The number depends on the machine
cutoff.max     = 33          # The first round filter, before normalization
cutoff.min     = 5
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



