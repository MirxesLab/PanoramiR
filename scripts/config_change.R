# ----- Path -----
dir.input        = '/path/to/input'        # absolute path
dir.out          = '/output/results/path'    # absolute path
file.samplesheet = '/input/samplesheet/file' # absolute path

# ----- Parameters -----
skip.samplesheet = 0           # the sample information start from which lines
arg.pipeline     = 'Cancer'    # ['Cancer', 'PanoramiR'] # 'Biofluid' haven't been tested by real data yet
is.RTsp          = FALSE       # Whether the filter samples by spike-in Ct values
is.basic         = FALSE       # For basic tier, no report generated
threshold.DE.dCt = 1
