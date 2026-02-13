# Mesoscopic Imaging Analysis Pipeline

This repository provides a workflow for preprocessing MiniMe mesoscopic imaging data with HPC.  
The pipeline includes ROI selection, blue/green frame separation, hemodynamics correction and alignment with the Allen Brain Atlas.

---

## Installation and Setup

### Step 1: Install Anaconda

Download and install **Anaconda for Python** from:  
https://www.anaconda.com/download

Anaconda is used to manage all dependencies and environments.

---

### Step 2: Create the Conda Environment

After installing Anaconda, open a terminal (or Anaconda Prompt) and navigate to this project folder.

Create a new environment named **analysis** using the provided `analysis.yml` file:

```bash
conda env create -f analysis.yml
````

Activate the environment:

```bash
conda activate analysis
```

---

## Step 3: ROI Selection and Alignment Preparation

Run the following script to perform ROI selection, separate blue/green frames, and generate the transformation matrix (`tform`) for alignment with the Allen Brain Atlas:

```bash
python MesoROI_square_rotation.py
```

This step performs:

* ROI selection
* Blue/green channel separation
* Generation of the alignment transformation (`tform`)

---

## Step 4: Upload the Pipeline to HPC and Configure Data Paths

After completing ROI selection locally, upload the meso_code folder to your HPC working directory.

#### 4.1 Update Data Paths in `defineIODirs.m`
Open `meso_code/defineIODirs.m`. Modify the following lines to match your HPC data directory:

```matlab
fixedInputDir = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';
fixedOutputDir = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training';
```

#### 4.2 Update Data Path in `separate_channel.sh`
Open `meso_code/cluster_code/separate_channel.sh`. Locate the `DATA_FOLDER` variable and update it:
```bash
DATA_FOLDER=/path/to/your/data_folder
```

## Step 5: Run the Processing Pipeline

Navigate to the `meso_code/cluster_code` folder:

```bash
cd meso_code/cluster_code
```

Run the main processing pipeline for a specific mouse dataset:

```bash
bash pipeline_MiniMe.sh <mouse_folder_id> <normalization_method>
```

Example:

```bash
bash pipeline_MiniMe.sh HD_Mouse5_0928_Stg3 dff
```

---

## Project Structure

```
├── analysis.yml                # Conda environment setup file
├── mouseMesoInfo.xlsx          # Example Mouse Info file for ROI selection and frame separation
├── tformMouse11.mat            # Example tform file for alignment with Allen Map
├── parcels_updated12522.mat    # Allen Map for the parcellation of Meso Imaging
├── MesoROI_square_rotation.py  # ROI selection and alignment preparation script
├── meso_code                   # Folder with auxiliary scripts
|   ├── MesoProcessing-master       # Helper function and older code
│   ├── cluster_code                # Main processing pipeline for HPC
│   |   ├── pipeline_MiniMe.sh          # Pipeline Code
└── README.md                   # Project documentation
```

---

## Notes

* Make sure the `analysis` environment is activated before running any scripts.
* All dependencies are installed automatically from `environment.yml`.
* Tested with Python 3.9+ and Anaconda 2024.

---

## License

This project is released for research use.
Please cite or acknowledge this repository if it contributes to your work.


