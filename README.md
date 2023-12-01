# Pipeline for Cancer Panel & PanoramiR Panel
## Prepare environment
The pipeline requires certain R packages. Please refer to `resources/installPackage.R` for installation instructions.

## Prepare sample infomation sheet
1. Sample information sheet is in `xlsx` format.
2. The sheet should have 7 columns with the following specifications:

    1. **First column** : labeled as `S/N`. This column records sample names defined by the customers. These sample names will be used by figures.
    2. **Second ~ Forth columns**: These must be labeled as `Comparison 1`, `Comparison 2`, and `Comparison 3`. They represent pairwise comparisons with "A" as the experiment group and "B" as the control group. ∆CT is calculated as $$avg(CT)_B - avg(CT)_A$$, where positive values indicate upregulation in the experiment group, and negative values indicate downregulation.
    3. **Fifth column**: IDs assigned to the qPCR results in `xlsx` format. Each ID should be unique to only one qPCR Excel file. Examples of good IDs are 001, 010, 100, while examples of bad IDs are 1, 10, 100.
    4. **Sixth column**: sample description, which can be left empty.
    5. **Seventh column**: This column must be labeled as `Sample Type`. The pipeline checks the sample type for each comparison, and within each comparison, the sample types should be the same.

## Prepare the configration files
Modify the paths and parameters in `config_change.R` to customize the pipeline for different datasets and situations.

## Run the pipeline `run_pipeline.R`
1. Change the diretory to the ones where you stored the scripts. 
2. Lines 7 to 10 contain the code to execute the pipeline. Please update the parameters used by the `rmarkdown::render` function:

    1. `input`: the `RMD` file used to compile the report
    2. `output_file`: the file name of the report
    3. `output_dir`: the directory where the report will be stored.
    4. `params`: the path of configration file prepared earlier.
```
rmarkdown::render(input = cancer, 
                  output_file = 'CancerPanel.html', 
                  output_dir = '/path/to/store/HTML',  # absolute path
                  params = list(config = '/path/to/config_change.R'))
```
## Trouble shooting
1. Please ensure that all required R packages are installed. Note that some packages may not be listed in `installPackage.R`. Please let me know if you encounter any missing packages.
2. Ensure you use absolute paths to guarantee the pipeline can locate all necessary files.
3. The pipeline relies on the 5th column in the sample information sheet to read qPCR results. Please double-check for **typos** and confirm that the names in the 5th column are **uniquely matched with the names of qPCR results**.

## Licenses of packages
| package      | License            | Version | URL                                         |
|--------------|--------------------|---------|---------------------------------------------|
| BiocManager  | Artistic-2.0       | 1.30.22 | https://bioconductor.github.io/BiocManager/ |
| rmarkdown    | GPL-3              | 2.25    | https://github.com/rstudio/rmarkdown        |
| tidyverse    | MIT + file LICENSE | 2.0.0   | https://tidyverse.tidyverse.org             |
| readxl       | MIT + file LICENSE | 1.4.3   | https://readxl.tidyverse.org                |
| stringi      | file LICENSE       | 1.7.12  | https://stringi.gagolewski.com/             |
| stringr      | MIT + file LICENSE | 1.5.0   | https://stringr.tidyverse.org               |
| plotly       | MIT + file LICENSE | 4.10.2  | https://plotly-r.com                        |
| Gmisc        | GPL (>= 3)         | 3.0.3   | https://gforge.se                           |
| glue         | MIT + file LICENSE | 1.6.2   | https://github.com/tidyverse/glue           |
| grid         | Part of R 4.3.2    | 4.3.2   | https://cran.r-project.org/web/packages/grid/index.html                                          |
| knitr        | GPL                | 1.44    | https://yihui.org/knitr/                    |
| kableExtra   | MIT + file LICENSE | 1.3.4   | http://haozhu233.github.io/kableExtra/      |
| rmdformats   | GPL (>= 2)         | 1.0.4   | https://github.com/juba/rmdformats          |
| reshape2     | MIT + file LICENSE | 1.4.4   | https://github.com/hadley/reshape           |
| org.Hs.eg.db | Artistic-2.0       | 3.17.0  | https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html                                          |
| topGO        | LGPL               | 2.52.0  | https://bioconductor.org/packages/release/bioc/html/topGO.html                                          |
| IRanges      | Artistic-2.0       | 2.34.1  | https://github.com/Bioconductor/IRanges   |
| pheatmap     | GPL-2              | 1.0.12  | https://cran.r-project.org/web/packages/pheatmap/index.html                                          |
| RColorBrewer | Apache License 2.0 | 1.1-3   | https://cran.r-project.org/web/packages/RColorBrewer/index.html                                          |
| cowplot      | GPL-2              | 1.1.1   | https://wilkelab.org/cowplot/               |