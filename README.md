Mesoscopic Imaging Analysis Pipeline

This repository provides a streamlined workflow for preprocessing and analyzing widefield mesoscopic imaging data, including ROI selection, channel separation, and alignment with the Allen Brain Atlas.

🧰 Installation and Setup
1. Install Anaconda

First, download and install Anaconda for Python
.
Anaconda provides an isolated environment and all necessary dependencies for the analysis.

2. Create the Conda Environment

After installing Anaconda, open a terminal (or Anaconda Prompt) and navigate to the directory containing this repository.
Then create the environment using the provided .yml file:

conda env create -f analysis.yml


Once the environment is created, activate it:

conda activate analysis

🧪 Step 1: ROI Selection and Alignment Preparation

Run the following script to perform ROI selection, separate blue and green frames, and generate the transformation matrix (tform) for alignment with the Allen Brain Atlas:

python MesoROI_square_rotation.py


This script will:

Allow you to manually select ROIs

Separate blue/green channel frames

Generate the alignment transformation (tform) file for later use

⚙️ Step 2: Run the Processing Pipeline

Navigate to the meso-aux-scripts folder:

cd meso-aux-scripts


Run the main processing pipeline for a specific mouse dataset using:

bash pipeline_id.sh <mouse_folder_id>


Replace <mouse_folder_id> with the ID or path of your mouse data folder.

Example:

bash pipeline_id.sh 230415_mouse11

📁 Project Structure
├── environment.yml             # Conda environment configuration
├── MesoROI_square_rotation.py  # ROI selection and alignment preparation script
├── meso-aux-scripts/           # Auxiliary scripts and pipelines
│   ├── pipeline_id.sh          # Main processing pipeline
│   └── ...
└── README.md                   # Project documentation

🧩 Requirements

All required dependencies are automatically installed via the environment.yml file.
You just need Anaconda and Python ≥ 3.9.
