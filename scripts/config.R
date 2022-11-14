# =========================================================================== #
# Configure file for panoramiR.Rmd
# =========================================================================== #

# Load Packages
# --------------------------- #
library(dplyr)
library(readxl)
library(stringi)
library(stringr)

library(org.Hs.eg.db)
library(topGO)

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(plotly) # For interactive
library(Gmisc)
library(glue)
library(grid)
library(cowplot)

library(knitr)
library(kableExtra)



# Set Directory
# --------------------------- #
dir.out.fig  = file.path(dir.out, 'figure')
dir.out.tbl  = file.path(dir.out, 'table')

dir.resource.fig = 'resources/figure'
dir.resource.tbl = 'resources/table'

# Set working directory
# setwd(dir.top)

# Create directory for output
dir.create(dir.out.tbl, recursive = TRUE, showWarnings = FALSE)

# Set file path
# --------------------------- #
# Pattern of input file names: 2020-11-16_093607 2007 001_Results.xlsx [Cancer]
# Pattern of input file names: 2020-11-16_093607 2007 001_002_Results.xlsx [Biofluid]
    # Experiment Time: 2020-11-16_093607
    # Unique sample ID: 2007 001 or 2007 001_002 [Recorded in sample manifest form]

##### Common Resources
file.ref.target  = file.path(dir.resource.tbl, 'Predicted_Targets_Human_Filtered_Rearranged_GeneID.txt') 
png.mSMRT.qPCR   = file.path(dir.resource.fig, 'mSMRT_qPCR.jpeg') # Haven't been inserted into report

##### Pipeline based resources
if(arg.pipeline == 'PanoramiR') { # PanoramiR Pipeline
    file.ref.miRList = file.path(dir.resource.tbl, 'miRList_PanoramiR.xlsx') # [Group, Position, miRBASE v22 Accession, miRNAID]
    file.ref.miRsp   = file.path(dir.resource.tbl, 'miRSpike.xlsx' )         # [Spike-in miRNA: Position, miRNA ID, SP, Group]
    # file.ref.miRthld = file.path(dir.resource.qc, 'miRThreshold.xlsx')      # [miRNA, Threshold]
    # Use Fixed threshold now 9 ~ 33
    png.workflow     = file.path(dir.resource.fig, 'workflow_PanoramiR.png')
    
} else if (arg.pipeline == 'Cancer') { # Cancer Pipeline
    file.ref.miRList = file.path(dir.resource.tbl, 'miRList_Cancer.xlsx')
    png.workflow     = file.path(dir.resource.fig, 'workflow_Cancer.png')
    name.sp          = c('Spike-in RNA Ctr 1', 'Spike-in RNA Ctr 2')
    name.ipc         = c('Inter-plate Calibrator 1', 'Inter-plate Calibrator 2')
    
} else if (arg.pipeline == 'Biofluid') { # Biofluid Pipeline
    file.ref.miRList = file.path(dir.resource.tbl, 'miRList_Biofluid.xlsx')
    png.workflow     = file.path(dir.resource.fig, 'workflow_Biofluid.png')
    name.sp          = c('Spike-in RNA Ctr 1', 'Spike-in RNA Ctr 2')
    name.ipc         = c('Inter-plate Calibrator 1', 'Inter-plate Calibrator 2')
}

# Set parameters
# --------------------------- #
if(!is.basic) {
    dir.create(dir.out.fig, recursive = TRUE, showWarnings = FALSE)
}

sd.ipc         = 0.5         # Standard deviation greater than 0.5 indicates pipetting error and caution should be taken when interpreting results
diff.sp        = 0.5         # only when is.RTsp = TRUE. larger difference is usually a sign of pipetting error.

result.skip    = 46          # The number depends on the machine
cutoff.max     = 33          # The first round filter, before normalization
cutoff.min     = 9           # Can be changed based on customer' requirement
cutoff.sp      = 32          # The second round filter, for spike-in normalization

threshold.ipc  = 25          # If IPC has Ct values greater than 25, it indicates qPCR reaction has failed. -> remove those samples
threshold.sp   = 30          # all Spike-Ins from all samples have Ct value lower than this value -> remove those samples

threshold.impute    = 0.1    # No more than 10% missing value in miRNA
threshold.DE.pvalue = 0.05   # The p value of T test
threshold.DE.dCt    = 1      
threshold.GO.pvalue = 0.01
GO.nodesize         = 6


# Set Color
col.compare.1 = brewer.pal(6, "Dark2")[1:2] # Green - Orange
col.compare.2 = brewer.pal(6, "Dark2")[3:4] # Purple - Pink
col.compare.3 = brewer.pal(6, "Dark2")[5:6] # Green - Yellow
col.others = c(brewer.pal(11,'RdYlBu')[1], '#001F36') # Red - Blue
col.grey   = brewer.pal(8, "Set2")[8] # grey
col.heatmap = brewer.pal(11, 'PRGn')[c(1, 11)] # Purple - Green
col.heatmap = c(col.heatmap[1], "#FFFFFF", col.heatmap[2])



