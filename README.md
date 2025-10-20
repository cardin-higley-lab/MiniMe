# Mesoscopic Imaging Analysis Pipeline

This repository provides a workflow for preprocessing and analyzing widefield mesoscopic imaging data.  
The pipeline includes ROI selection, blue/green frame separation, and alignment with the Allen Brain Atlas.

---

## Installation and Setup

### Step 1: Install Anaconda

Download and install **Anaconda for Python** from:  
https://www.anaconda.com/download

Anaconda is used to manage all dependencies and environments.

---

### Step 2: Create the Conda Environment

After installing Anaconda, open a terminal (or Anaconda Prompt) and navigate to this project folder.

Create a new environment named **analysis** using the provided `environment.yml` file:

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

## Step 4: Run the Processing Pipeline

Navigate to the `meso-aux-scripts` folder:

```bash
cd meso-aux-scripts
```

Run the main processing pipeline for a specific mouse dataset:

```bash
bash pipeline_id.sh <mouse_folder_id>
```

Example:

```bash
bash pipeline_id.sh HD_Mouse5_0928_Stg3
```

---

## Project Structure

```
├── analysis.yml                # Conda environment setup file
├── mouseMesoInfo.xlsx          # Example Mouse Info file for ROI selection and frame separation
├── tformMouse11.mat            # Example tform file for alignment with Allen Map
├──parcels_updated12522.mat     # Allen Map for the parcellation of Meso Imaging
├── MesoROI_square_rotation.py  # ROI selection and alignment preparation script
├── meso-aux-scripts/           # Folder with auxiliary scripts
│   ├── pipeline_id.sh          # Main processing pipeline
│   └── ...
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


