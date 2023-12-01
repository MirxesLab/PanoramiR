# Pipeline for Cancer Panel & PanoramiR Panel
## Prepare Environment
The pipeline requires some R packages. Please refer to `resources/installPackage.R`.

## Prepare sample infomation sheet
1. Sample information sheet is in `xlsx` format.
2. There should be 7 columns in the sample information sheet.

    1. First column: sample names defined by the customers. These sample names will be used by figures.
    2. Second ~ forth columns must be `Comparison 1`, `Comparison 2`, `Comparison 3`. They are pairwise comparison, using "A" as the experiment group and "B" as the control group. ∆CT is calculated by 
        $$
        avg(CT)_B - avg(CT)_A
        $$
        Therefore, possitive values of ∆CT represents miRNAs are upregulated in experiment groups and negative values of ∆CT represents miRNAs are downregulated in experiment groups. 
    
    3. Fifth column should be the IDs given to the qPCR results. Pipeline uses these ID as "patterns" to read in qPRC results. **Each pattern should be uniqe to only 1 qPRC excel file**. Good IDs: 001, 010, 100. Bad IDs: 1, 10, 100. 
    4. Sixth column is sample description which can be left empty.
    5. Seventh column is `Sample Type` (cannot change the name of this column). The pipeline will check the sample type of each comparison. In each comparison, the sample type should be the same.

## Prepare the configration files
Change the paths and parameters in `config_change.R` to apply the pipeline to different datasets and different situations. 

## Run the pipeline `run_pipeline.R`
1. Change the diretory to the ones in which you stored the scripts. 
2. For line 7 to line 10 is the code to run the pipeline. Please change the parameters used by the `rmarkdown::render` function.

    1. `input`: the `RMD` file used to compile the report
    2. `output_file`: the file name of the report
    3. `output_dir`: the directory stored the report
    4. `params`: the path of configration file prepared above.
```
rmarkdown::render(input = cancer, 
                  output_file = 'CancerPanel.html', 
                  output_dir = '/path/to/store/HTML',  # absolute path
                  params = list(config = '/path/to/config_change.R'))
```
## Trouble shooting
1. Please use absolute path to make sure that the pipeline can find all necessary files.
2. The pipeline uses the 5th column in sample information sheet to read in qPCR results. Please make sure that **there is no typo** and the names in 5th column is **the unique part of the names of qPCR results**. 
