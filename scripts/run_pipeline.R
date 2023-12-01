library(rmarkdown)

dir.script="/path/to/RMD/files"                    # absolute path
panoramiR = file.path(dir.script, "panoramiR.Rmd") 
cancer = file.path(dir.script, "CancerPanel.Rmd")

rmarkdown::render(input = cancer, 
                  output_file = 'CancerPanel.html', 
                  output_dir = '/path/to/store/HTML',  # absolute path
                  params = list(config = '/path/to/config_change.R'))

