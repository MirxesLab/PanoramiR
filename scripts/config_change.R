# ----- Path -----
dir.input        = '/Users/plateau/OneDrive - MiRXES/1. Ongoing task/Preferred Medicine_2022/Results/Final report/qPCR data/'
dir.out          = '/Users/plateau/OneDrive - MiRXES/1. Ongoing task/Preferred Medicine_2022/QualityCheck/Res_Norm'
#file.samplesheet = file.path(dir.input, 'samplesheet_Tall.xlsx')
file.samplesheet = '/Users/plateau/OneDrive - MiRXES/1. Ongoing task/Preferred Medicine_2022/QualityCheck/samplesheet.xlsx'

# ----- Parameters -----
skip.samplesheet = 0           # the sample information start from which line
arg.pipeline     = 'Cancer'    # ['Cancer', 'Biofluid', 'PanoramiR']
is.RTsp          = FALSE       # Whether the filter samples by spike-in Ct values
is.basic         = T       # For basic tier, no report generated
