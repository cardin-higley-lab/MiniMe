# This code is for analysis of the regular folder behaviors in hpc
# need more modification
import deeplabcut
import matplotlib
import time
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import cv2
import csv
import os





root_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training'
path_config_file='/gpfs/gibbs/project/higley/hd362/Head_Direction_Vis_Sti_On-Hao-2024-11-16/config.yaml'

def find_main_folder(root_folder):
# List to store only directories
    main_folders = []

    # Loop through all items in the directory
    for item in os.listdir(root_folder):
        # Construct the full path of the item
        main_folders.append(os.path.join(root_folder, item))

    return main_folders


# Function to extract the clips
def extract_and_save_clips(video_path, start_indices, output_folder):
    # Open the video file
    cap = cv2.VideoCapture(video_path)
    
    if not cap.isOpened():
        print("Error: Could not open video.")
        return

    fps = cap.get(cv2.CAP_PROP_FPS)
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))  # Original width
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))  # Original height

    # Loop through each start index
    for start_index in start_indices:
        start_index = start_index - 30*4
        
        # Set the video capture to the start index
        cap.set(cv2.CAP_PROP_POS_FRAMES, start_index)

        # Specify the output video path with .mp4 extension
        output_video_path = f'{output_folder}/clip_{start_index}.mp4'
        
        # Use the 'mp4v' codec (H.264) for mp4 format
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Codec for .mp4 format
        
        # Create a VideoWriter object with the specified codec, FPS, and frame size
        out = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

        # Write the frames for this clip
        for _ in range(600):  # Write 210 frames (adjust if needed)
            ret, frame = cap.read()
            if not ret:
                print(f"End of video reached before completing clip starting at frame {start_index}.")
                break
            out.write(frame)

        # Release the VideoWriter for this clip
        out.release()
        print(f"Saved clip: {output_video_path}")

    # Release the VideoCapture object
    cap.release()
    

# find all the start_frame file in the main folder
def start_frame_csv(date_folder):
    csv_files = os.path.join(date_folder, 'process', 'behavior_start.csv')
    return csv_files

'''
def create_locs_folder(main_folder):
    # Create folders locs1, locs2, locs3, locs4
    for i in range(1, 5):  # This will loop from 1 to 4
        folder_name = os.path.join(main_folder, f"locs{i}")
        
        # Create the main folder
        os.makedirs(folder_name, exist_ok=True)

        # Create the two subfolders inside each main folder
        subfolders = ["Behavior Clips", "Angles"]
        for subfolder in subfolders:
            os.makedirs(os.path.join(folder_name, subfolder), exist_ok=True)
'''
## Previous code did not find the correct start frames
# So I use MATLAB code to find the start frames
# Load them here and create all the videos
def save_clips_all_locs(video_path, date_folder, csv_files): 
    # Iterate over locs1, locs2, locs3, and locs4
    
    # Find the output folder
    output_folder = os.path.join(date_folder, 'behavior', 'Behavior Clips')
    
    os.makedirs(output_folder, exist_ok=True)
    
    start_frame_filename = csv_files  # Use the single element directly

    # Load the CSV file into a list of integers (start frames)
    start_frames = []
    # Open the CSV file
    with open(start_frame_filename, newline='', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile)
        start_frames = [int(item) for row in csvreader for item in row]
    
    # Function to extract and save clips
    extract_and_save_clips(video_path, start_frames, output_folder)
        
def find_mid_behavior(date_folder):
    # List all files and directories in the specified directory
    files = os.listdir(date_folder)

    # Filter files that start with 'mid'
    mid_behavior = [file for file in files if file.startswith('mid')]
    mid_behavior_fullpath = os.path.join(date_folder, mid_behavior[0])

    return mid_behavior_fullpath

def find_behavior_clips_folder(root_dir):
    # List to store the paths of all "Behavior Clips Long" folders found
    found_folders = []
    
    # Walk through all directories and subdirectories
    for root, dirs, files in os.walk(root_dir):
        # If "Behavior Clips Long" is in the current directory's list of subdirectories
        if "Behavior Clips" in dirs:
            # Append the full path to the "Behavior Clips Long" folder to the list
            found_folders.append(os.path.join(root, 'Behavior Clips'))
    
    # Return the list of all found folders
    return found_folders

def DLC(videoFolder, path_config_file):

    mp4_files = []
    for root, dirs, files in os.walk(videoFolder):
        for file in files:
            if file.endswith('.mp4'):
                mp4_files.append(os.path.join(root, file))

    path_video = mp4_files

    deeplabcut.evaluate_network(path_config_file)
    deeplabcut.analyze_videos(path_config_file, path_video, save_as_csv=True)

main_folder = find_main_folder(root_folder)
for date_folder in main_folder:
    
    # create_locs_folder(date_folder)
    video_path = find_mid_behavior(date_folder)
    if len(video_path)>0: 
        print(video_path)
        csv_files = start_frame_csv(date_folder)
        save_clips_all_locs(video_path, date_folder, csv_files)


for date_folder in main_folder:
    #print(date_folder)
    clipFolder = find_behavior_clips_folder(date_folder)
    for item in clipFolder:
        DLC(item, path_config_file)

'''
date_folder = main_folder[13]
video_path = find_mid_behavior(date_folder)
csv_files = start_frame_csv(date_folder)
save_clips_all_locs(video_path, date_folder, csv_files)

print(csv_files[-1:])

for file in csv_files[-1:]:
    print(file)
''' 



# extract_and_save_clips(video_path, start_frames)


# extract the hemo corrected video
# start_frame_hemo = [frames - 2000 for frames in start_frames]
# hemo_video_path = 'W:/Hao/Globus/2024_12_12_14_52/output_video_range3_mask_600000.avi'
# output_folder = 'W:/Hao/Globus/2024_12_12_14_52/Neuro Clips'
# extract_and_save_clips(hemo_video_path, start_frame_hemo, output_folder)
'''
# Plot one of the column to check the shape and start
if df.shape[1] >= 5:
    # Plot the fifth column
    plt.plot(df.iloc[:, 2])  # iloc[:, 4] accesses the fifth column (0-indexed)
    plt.title('Fifth Column Data')
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.xlim(0, 80000)
    plt.show()
else:
    print("The CSV file does not have a fifth column.")

# Plot the one column with x axis of the first column
if df.shape[1] >= 5:
    # Extract the first and fifth columns
    x = pd.to_numeric(df.iloc[:, 0], errors='coerce')  # Convert to numeric, coerce errors
    y = pd.to_numeric(df.iloc[:, 2], errors='coerce')  # Convert to numeric, coerce errors

    # Drop NaN values that may have resulted from coercion
    valid_indices = x.notna() & y.notna()
    x = x[valid_indices]
    y = y[valid_indices]

    # Create the plot
    plt.plot(x, y)
    plt.title('Fifth Column vs First Column')
    plt.xlabel('First Column')
    plt.ylabel('Fifth Column')
    plt.grid()
    plt.show()
else:
    print("The CSV file does not have enough columns.")

# Find all the start index of the Screen Vis Stimuli Start
if df.shape[1] >= 5:
    # Extract the fifth column
    fifth_column = df.iloc[:, 2]

    # Define a threshold to identify the start of the signal
    threshold = 0.5  # Adjust this based on your signal characteristics

    # Find the start of the signal
    start_index = fifth_column[fifth_column < threshold].index.tolist()

    if start_index is not None:
        print("Start indices of the signal:", start_index)
    else:
        print("No signal found above the threshold.")
else:
    print("The CSV file does not have a fifth column.")
    
# print(start_index[2])

'''




