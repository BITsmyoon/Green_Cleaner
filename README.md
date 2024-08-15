
# Green Cleaner: A Methodology for Removing Contaminants from Urine 16S rRNA Sequencing Data

## Introduction

**Green Cleaner** is a methodology designed to remove contaminants from urine 16S rRNA sequencing data. All scripts in this repository have been thoroughly tested.

### Repository Structure

- **DATA**: Contains all the necessary data for analysis.
- **SCRIPT**: Includes all analysis code, organized in the order of execution, consisting of R and shell scripts.
- **OUTPUT**: Stores the results of the analysis.

### Setting Up the Analysis Environment

The analysis environment is defined in `data/Green_Cleaner_conda_env.yml`. To set up the environment, run the following commands:

    conda env create -f ${path_variable}/Green_Cleaner_conda_env.yml
    conda activate Green_Cleaner_conda

Here, `${path_variable}` should be set to the location where you downloaded the repository and the data directory. For example:

    path_variable=${download_path}/data

### Running the Analysis
To perform the analysis, navigate to the script directory and use Rscript or bash to execute the scripts, passing ${path_variable} as a variable. For example:

    cd ${path_variable}/../script
    Rscript 01.Green_Cleaner_decontamination_in_batch.R ${path_variable}

This example demonstrates how to execute the first script in the analysis pipeline. Replace the script name as needed for subsequent steps.
