# ----- Path -----
dir.input        = 'input_sample/input_cancer'
dir.out          = 'cancer_output'
file.samplesheet = file.path(dir.input, 'Sample Manifest Form.xlsx')

# ----- Parameters -----
arg.pipeline   = 'Cancer'    # ['Cancer', 'Biofluid', 'PanoramiR']
is.RTsp        = FALSE       # Whether the filter samples by spike-in Ct values
is.basic       = FALSE       # For basic tier, no report generated
