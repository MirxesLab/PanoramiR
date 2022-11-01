# ----- Path -----
dir.input        = 'input_ncRNA'
dir.out          = 'panoramiR_output'
file.samplesheet = file.path(dir.input, 'Sample Manifest Form.xlsx')

# ----- Parameters -----
arg.pipeline   = 'PanoramiR'    # ['Cancer', 'Biofluid', 'PanoramiR']
is.RTsp        = FALSE       # Whether the filter samples by spike-in Ct values
is.basic       = FALSE       # For basic tier, no report generated
