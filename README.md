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
