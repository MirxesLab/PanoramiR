# ----- Path -----
dir.input        = 'P2203_Raw_data'
dir.out          = 'P2203_results'
file.samplesheet = file.path(dir.input, 'samplesheet.xlsx')

# ----- Parameters -----
skip.samplesheet = 0           # the sample information start from which line
arg.pipeline     = 'PanoramiR'    # ['Cancer', 'Biofluid', 'PanoramiR']
is.RTsp          = FALSE       # Whether the filter samples by spike-in Ct values
is.basic         = FALSE       # For basic tier, no report generated
