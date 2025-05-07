<p align="center">
  <img width="717" alt="Image" src="https://github.com/user-attachments/assets/c97fe22d-f5f3-4f2d-a6fc-a93669f5ccf4" />
</p>


## Overview
VESPER is a variant interpretation pipeline that profiles each read of input SAM files to identify single nucleotide polymorphisms (SNPs) using the alignment's CS tag. It classifies detected variants on a per-read basis as single-nucleotide variants (SNVs), multi-nucleotide variants (MNVs), or a combination of both. These variants are then reformatted to follow HGVS nomenclature before utilising Ensembl's Variant Effect Predictor (VEP) to assess the functional consequences and potential pathogenicity of each variant, with particular emphasis on their implications in Dilated Cardiomyopathy. VESPER also includes downstream analyses to compare variant frequencies across different cell populations, elucidating population-specific variant distributions. High level architecture of the workflow is as shown in the schematic below.

## Vesper's Workflow
<p align="center">
  <img width="1469" alt="Image" src="https://github.com/user-attachments/assets/f0c57602-489c-4867-a343-2ee19b4c3b5c" />
</p>

#### Pre-Requisites
* Python 3.9
    * Ensure that Python 3.9 is installed on your system. Install it otherwise before proceeding.
* Singularity image of [Ensembl VEP](https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#singularity)
    * Refer to link for installation instructions

#### Installation of Vesper
1. Navigate into directory where vesper is to be  installed in
   ```
   wget vesper link (to be finalised)
   tar -xvzf vesper.tar.gz
#### Set Up Vespers's Environment
2. Navigate into vesper directory
   ```
   cd vesper
3. Create and activate vesper_env
   ```
   python3.9 -m venv vesper_env         # create new python virtual env
   source vesper_env/bin/activate       # activate vesper_env
4. Install packages
   ```
   python3 setup_vesper_env.py 
#### Configure vesper_config.yaml file
5. Edit minigene_config.yaml file to suit your data

    Option 1: use vim directly in terminal to edit
    ```
    vim minigene_config.yaml
   
   # press 'i' to edit and make changes; when done, press 'esc', ':wq' and hit 'enter' to save and exit
   ```
    
    Option 2: edit locally on your preferred text editor

#### Execution
```
python3 vesper.py -c vesper_config.yaml
