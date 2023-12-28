library(rmarkdown)

dir.script="/path/to/RMDfiles"                    # absolute path
panoramiR = file.path(dir.script, "panoramiR.Rmd") 
cancer = file.path(dir.script, "CancerPanel.Rmd")

## to process results of cancer panel 
rmarkdown::render(input = cancer, 
                  output_file = 'CancerPanel.html', 
                  output_dir = '/path/to/store/HTML',  # absolute path
                  params = list(config = '/path/to/config_change.R'))

## to process results of panoramiR
# rmarkdown::render(input = panoramiR, 
#                   output_file = 'panoramiR.html', 
#                   output_dir = '/path/to/store/HTML',  # absolute path
#                   params = list(config = '/path/to/config_change.R'))



