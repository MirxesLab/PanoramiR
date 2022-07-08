# =========================================================================== #
# Configure file for panoramiR.Rmd
# =========================================================================== #

# Load Packages
library(dplyr)
library(readxl)
library(stringi)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(plotly)
library(knitr)
library(kableExtra)

# Set Path
dir.top     = '/Users/plateau/Documents/GitHub/PanoramiR'
dir.out.fig = 'output/figures'
dir.out.tbl = 'output/tables'

setwd(dir.top)

file.input.data  = file.path('input', 
                             grep('*_Results.xlsx', 
                                  list.files('input', all.files = FALSE), 
                                  value = TRUE))    #[*.Results.xlsx]
file.samplesheet = 'input/Sample Manifest Form.xlsx'
file.ref.miRList = 'resources/miRList.xlsx'         # All miRNA: Group, Position, miRBASE v22 Accession, miRNAID
file.ref.miRsp   = 'resources/miRSpike.xlsx'        # Spike-in miRNA: Position, miRNA ID, SP, Group
file.ref.miRthld = 'resources/miRThreshold.xlsx'    # Only for PanoramiR: miRNA, Threshold
file.ref.RNAvol  = 'resources/RNAvolume.xlsx'

# Create directory for output
dir.create(dir.out.fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir.out.tbl, recursive = TRUE, showWarnings = FALSE)

# Set parameters
is.basic       = TRUE        # Basic worflow execute Global Normalization. Non-basic workflow also execute normalization based on stable miRNA
is.RTsp        = FALSE       # If the spike-in is Reverse Transcript Spike-in

result.skip    = 46          # The number depends on the machine
cutoff.max     = 33          # The first round filter, before normalization
cutoff.min     = 5
cutoff.sp      = 32          # The second round filter, for spike-in normalization

threshold.DE.pvalue = 0.05   # The p value of T test
threshold.DE.log2fc = 1      # The log2 fold change of geometric mean

# Set Color
col.compare.1 = brewer.pal(6, "Dark2")[1:2] # Green - Orange
col.compare.2 = brewer.pal(6, "Dark2")[3:4] # Purple - Pink
col.compare.3 = brewer.pal(6, "Dark2")[5:6] # Green - Yellow
col.others = brewer.pal(11,'RdYlBu')[c(1, 11)] # Red - Blue
col.grey   = brewer.pal(8, "Set2")[8] # grey





