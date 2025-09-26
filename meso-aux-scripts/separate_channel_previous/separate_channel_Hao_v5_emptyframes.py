# -*- coding: utf-8 -*-
"""
Project Name: Hemodynamics correction pipline

Created on Thu Aug  1 14:46:06 2024
Modified on Monday Mar 31 11:10 2025

@author: Clayton
@update: Hao
"""

import sys
import glob
import cv2
import re
import numpy as np
import json
from scipy.signal import butter
from scipy.signal import sosfiltfilt
from scipy.io import savemat, loadmat
import keyboard
from natsort import natsorted
import os
from tqdm import tqdm
import skimage
import pandas as pd
from collections import defaultdict

# Main folder and mice folder
movement_artifact_thresh = 3
main_folder = '/vast/palmer/scratch/higley/hd362/contrast_adjusting/contrast_adjusting'
# main_folder = r'W:\Hao\Vis Stimuli Behavior\contrast_adjusting'
mouse_folder = glob.glob(os.path.join(main_folder,'mouse*'))

# mice info input
mouse_input = sys.argv[1]
date_input = sys.argv[2]

'''
# For temporary testing
main_folder = 'W:/Hao/Vis Stimuli Behavior/contrast_adjusting'
mouse_folder = glob.glob(os.path.join(main_folder,'mouse*'))
mouse_input = 'mouse2'
date_input = '0221'

'''
# Seclect the mouse and date file
mouse_select = [folder for folder in mouse_folder if mouse_input in folder]
date_folder = glob.glob(os.path.join(mouse_select[0], '*/'))
date_select = [folder for folder in date_folder if date_input in folder]

'''
for mfolder in mouse_folder:
    
    subfolder = glob.glob(os.path.join(mfolder, '*/'))
    date_folders.append(subfolder)
'''


# find the meso data folder
meso_folder = glob.glob(os.path.join(date_select[0], '*', 'My_V4_Miniscope'), recursive=True)
print(meso_folder)

# funtion for finding the frame drop location
def findFrameDrop(meso_folder):

    # find the frame drop locs and num of frames
    timeStampsFile = glob.glob(os.path.join(meso_folder[0], 'timeStamps*'))
    timeStampsDf = pd.read_csv(timeStampsFile[0])
    timeStamps = timeStampsDf.values
    timeDiff = np.diff(timeStamps[:, 1])
    meanDiff = np.mean(timeDiff)


    # Initialize variables
    frameCount = 0
    frameSave = []
    dropFrame = []

    # Iterate through time differences
    for n, diff in enumerate(timeDiff):
        numDiff = round((diff - meanDiff) / meanDiff)
        if numDiff == 0:  # Changed from ~numDiff to numDiff == 0
            frameCount += 1
            frameSave.append(frameCount-1)
        else:
            for _ in range(numDiff):  # Changed m to _
                frameCount += 1
                dropFrame.append(frameCount-1)
            frameCount += 1
            frameSave.append(frameCount-1)

    # Create the frameCheck DataFrame
    frameCheck = pd.DataFrame(np.ones((frameCount, 2), dtype=int))  # Initialize with ones
    frameCheck.iloc[dropFrame, :] = 0  # Set drop frames to 0
    frameCheck.iloc[frameSave, 1] = np.arange(1, len(frameSave) + 1)  # Assign frame numbers

    fDropData = []  # list of dicts to mimic struct array

    for i in range(len(dropFrame)):
        drop_idx = dropFrame[i]
        pre_idx = drop_idx - 1
        
        preFrame = frameCheck.iloc[pre_idx, 1]  # get value before dropped frame

        if preFrame:
            # start a new drop record
            fDropData.append([preFrame,  1])
        else:
            # increase the drop count of the last record
            if fDropData:
                fDropData[-1][1] += 1

    # Remove entries where 'drop_locs' is empty or evaluates to False (e.g., 0, None)
    fDropData = [entry for entry in fDropData if entry[0]]
    
    # Save to .mat file
    fDropFolder = os.path.join(os.path.dirname(meso_folder[0]), 'frameDrop.mat')
    savemat(fDropFolder, {'frameDrop': fDropData})
    
    return fDropData

'''
fDropFolder = os.path.join(os.path.dirname(meso_folder[0]), 'frameDrop.mat')
fDrop = loadmat(fDropFolder)
fDropData = fDrop['frameDrop'][0]
'''

# Find the frame drop location and numbers
fDropData = findFrameDrop(meso_folder)
fDropList = []
fList = []
for item in fDropData:

    frame_locs = item[0] # e.g., '6.avi'# e.g., 538
    count = item[1]      # e.g., 1
    fDropList.append([frame_locs, count])
    fList.append(frame_locs)
    print(f"Frame: {frame_locs}, Count: {count}")


'''
# temporary num for testing purpose
num_folder = 7
# Find mouse ID from the filename
'''

# Find mice ID for cropping
parts = os.path.normpath(meso_folder[0]).split(os.sep)
mouseID = int(parts[8][5:])
print(mouseID) 

# In meso data folder, find all the .avi video path   
data_dir = meso_folder[0]
data_files = natsorted(glob.glob(os.path.join(data_dir,'*.avi')))
print(data_dir)

# extract the frame_rate, roi, and filenumber information
with open(os.path.join(data_dir,'metaData.json')) as f:
    meta_data = json.load(f)

frame_rate = meta_data['frameRate']
roi = meta_data['ROI']
frames_per_file = meta_data['framesPerFile']

output_shape = (256,256)

mouseInfofile = glob.glob(os.path.join(main_folder, mouse_input, date_input, 'mouseMesoInfo*'))
print(mouseInfofile)
#mouseInfofile = glob.glob(os.path.join(main_folder, '*', 'mouseMesoInfo*'))
#mouseInfofile = 'W:/Hao/Vis Stimuli Behavior/contrast_adjusting/process/mouseMesoInfo.xlsx'
mouseInfo = pd.read_excel(mouseInfofile[0])
row_num = mouseID-2
print(mouseInfo.loc[row_num].x_start)

# cropping params
y_start = mouseInfo.loc[row_num].y_start
y_end = mouseInfo.loc[row_num].y_end

x_start = mouseInfo.loc[row_num].x_start
x_end = mouseInfo.loc[row_num].x_end

# seperate G/R channel params
x1 = mouseInfo.loc[row_num].x1
x2 = mouseInfo.loc[row_num].x2
y1 = mouseInfo.loc[row_num].y1
y2 = mouseInfo.loc[row_num].y2
print(x1, x2, y1, y2)


def is_bright_enough(frame, x1, x2, y1, y2, brightness_threshold=230):
    # Extract the region of interest (ROI)
    roi = frame[y1:y2, x1:x2]
    
    # Calculate the average brightness of the ROI
    avg_brightness = np.mean(roi)
    
    # Return True if the average brightness of the ROI is less than the threshold
    return avg_brightness > brightness_threshold

def seperateGR(data_files, data_dir, seperateIndB, seperateIndG, frame_count, globalB, globalG, global_count):
    # Iterate over each AVI file in the folder and process them
    B_frames = []
    G_frames = []
    for i, filename in enumerate(data_files):
        video_path = os.path.join(data_dir, filename)
            
        # Open the video
        cap = cv2.VideoCapture(video_path)

        # Check if video opened successfully
        if not cap.isOpened():
            print(f"Error: Could not open video {filename}.")
            continue

        # Skip first 4 frames only for the first video
        if frame_count == 0:
            for _ in range(4):
                cap.read()  # Skip the first 4 frames
            frame_count = 4
            global_count = 4
            

        # Read frames and process
        while True:
            ret, frame = cap.read()
            if not ret:
                break  # Exit if no more frames are available

            frame_count += 1
            global_count += 1

            # Check if the specified region meets the brightness condition
            if is_bright_enough(frame, x1, x2, y1, y2):
                # Crop the frame to the specified final region (x1, x2, y1, y2)
                cropped_frame = frame[y_start:y_end, x_start:x_end]

                # Convert the frame to grayscale
                grayscale_frame = cv2.cvtColor(cropped_frame, cv2.COLOR_BGR2GRAY).astype(np.float64)

                # Resize the grayscale frame to 256x256
                resized_frame = cv2.resize(grayscale_frame, output_shape)

                # Append the resized grayscale frame to the selected frames list
                G_frames.append(resized_frame)
                seperateIndG.append(frame_count)
                globalG.append(global_count)
            elif np.mean(frame[y1:y2, x1:x2])>20:
                cropped_frame = frame[y_start:y_end, x_start:x_end]

                # Convert the frame to grayscale
                grayscale_frame = cv2.cvtColor(cropped_frame, cv2.COLOR_BGR2GRAY).astype(np.float64)

                # Resize the grayscale frame to 256x256
                resized_frame = cv2.resize(grayscale_frame, (256, 256))

                # Append the resized grayscale frame to the selected frames list
                B_frames.append(resized_frame)
                seperateIndB.append(frame_count)
                globalB.append(global_count)
                
            else:
                if frame_count>1:
                    # If there is a blank frame, check the previous frame
                    # if the previous is G, then blank is blue
                    if seperateIndG[-1] > seperateIndB[-1]:
                        B_frames.append(B_frames[-1])
                        seperateIndB.append(frame_count)
                        globalB.append(global_count)
                    else:
                        G_frames.append(G_frames[-1])
                        seperateIndG.append(frame_count)
                        globalG.append(global_count)
                        
            # Check the frame drop before looping to next frame
            if frame_count in fList:
                for row in fDropList:
                    if row[0] == frame_count:
                        # get the number of dropped frames
                        num_drop = row[1]
                        
                        # According to the previous frame, copy the number of previous frames
                        if seperateIndG[-1] > seperateIndB[-1]:
                            for m in range(num_drop):
                                global_count += 1
                                if m%2 == 0:                                        
                                    B_frames.append(B_frames[-1])
                                    globalB.append(global_count)
                                    
                                else:
                                    G_frames.append(G_frames[-1])
                                    globalG.append(global_count)
                                    
                        else:
                            for m in range(num_drop):
                                global_count += 1
                                if m%2 == 0:                                        
                                    G_frames.append(G_frames[-1])
                                    globalG.append(global_count)
                                    
                                else:
                                    B_frames.append(B_frames[-1])
                                    globalB.append(global_count)
                                      
        # Release the video object
        cap.release()
    return B_frames, seperateIndB, G_frames, seperateIndG, frame_count, globalB, globalG, global_count


def insertEmpty(data_files, data_dir, frame_count, global_count):
    # Iterate over each AVI file in the folder and process them
    full_frames = []
    for i, filename in enumerate(data_files):
        video_path = os.path.join(data_dir, filename)
            
        # Open the video
        cap = cv2.VideoCapture(video_path)
        frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
        frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
        
        empty_frame = np.zeros(output_shape, dtype=np.uint8)

        # Check if video opened successfully
        if not cap.isOpened():
            print(f"Error: Could not open video {filename}.")
            continue

        # Skip first 4 frames only for the first video
        if frame_count == 0:
            for _ in range(4):
                cap.read()  # Skip the first 4 frames
            frame_count = 4
            global_count = 4
            

        # Read frames and process
        while True:
            ret, frame = cap.read()
            if not ret:
                break  # Exit if no more frames are available

            frame_count += 1
            global_count += 1

            # Crop the frame to the specified final region (x1, x2, y1, y2)
            cropped_frame = frame[y_start:y_end, x_start:x_end]

            # Convert the frame to grayscale
            grayscale_frame = cv2.cvtColor(cropped_frame, cv2.COLOR_BGR2GRAY).astype(np.uint8)

            # Resize the grayscale frame to 256x256
            resized_frame = cv2.resize(grayscale_frame, output_shape)
            
            full_frames.append(resized_frame)

                        
            # Check the frame drop before looping to next frame
            if frame_count in fList:
                for row in fDropList:
                    if row[0] == frame_count:
                        # get the number of dropped frames
                        num_drop = row[1]
                           
                        for m in range(num_drop):
                            global_count += 1
                            full_frames.append(empty_frame)                                 
                                      
        # Release the video object
        cap.release()
    return np.array(full_frames), frame_count, global_count


def savecolumnmat(output_file, frames_full, colInd, chInd):
    # Prepare the slice data for saving (transpose if necessary)
    slice_data = frames_full[:,:,colInd].T

    ch1_filename = f"{output_file}/RawDemixedCH{chInd}Col{colInd}.mat"
    if os.path.exists(ch1_filename):
        # Load the existing data
        existing_data_ch1 = loadmat(ch1_filename)
        existing_column_ch1 = existing_data_ch1.get("column", np.array([]))
        
        # Append the new data to the existing data
        updated_column_ch1 = np.concatenate([existing_column_ch1, slice_data], axis=1)
    else:
        # If file doesn't exist, set the new data as the column
        updated_column_ch1 = slice_data

    # Save (or overwrite) the updated data for channel 1
    savemat(ch1_filename, {"column": updated_column_ch1})
    #print(f"Saved or updated {ch1_filename}")

# Function to process files in a chunk
def process_files(files_chunk, data_dir, output_file, B_frames_full, G_frames_full, seperateIndB, seperateIndG, frame_count, globalB, globalG, global_count):
    # Your processing code here (e.g., reading, analyzing files, etc.)
    [B_frames, seperateIndB, G_frames, seperateIndG, frame_count, globalB, globalG, global_count] = seperateGR(files_chunk, data_dir, seperateIndB, seperateIndG, frame_count, globalB, globalG, global_count)
    
    # Calculate the difference between two channels and add frames to the less channel
    frameDiff = abs(len(globalB)-len(globalG))
    print(f"The difference between two channels is: {frameDiff}")
    if len(globalB) > len(globalG):
        for n in range(frameDiff):
            global_count += 1
            G_frames.append(G_frames[-1])
            globalG.append(global_count)
    else:
        for n in range(frameDiff):
            global_count += 1
            B_frames.append(B_frames[-1])
            globalB.append(global_count)
    
    # After seperation, save the frames into the full data
    if len(B_frames_full) == 0 or len(G_frames_full) == 0:
        B_frames_full = np.array(B_frames)
        G_frames_full = np.array(G_frames)
    else:
        B_frames_full = np.concatenate((B_frames_full, np.array(B_frames)), axis=0)
        G_frames_full = np.concatenate((G_frames_full, np.array(G_frames)), axis=0)



    # Convert processed data to a numpy array (or any format you need)
    # Ensure the output directory exists
    os.makedirs(output_file, exist_ok=True)

    # Save the full data into .mat format
    # Loop through the 256 slices
    for colInd in range(256):
        # Prepare the slice data for saving (transpose if necessary)
        savecolumnmat(output_file, B_frames_full, colInd, '1')
        savecolumnmat(output_file, G_frames_full, colInd, '2')
    return frame_count, seperateIndB, seperateIndG, globalB, globalG, global_count

def saveseqindex(sep_ind_file, seqIndex, chInd):
    os.makedirs(sep_ind_file, exist_ok=True)
    seq_filename = f"{sep_ind_file}/sepIndex{chInd}.mat"
    savemat(seq_filename, {"seqIndex": seqIndex})
    
def saveglobalindex(global_ind_file, globalIndex, chInd):
    os.makedirs(global_ind_file, exist_ok=True)
    seq_filename = f"{global_ind_file}/globalIndex{chInd}.mat"
    savemat(seq_filename, {"globalIndex": globalIndex})
    
globalB = []
globalG = []
global_count = 0
seperateIndB = []
seperateIndG = []
B_frames_full = []
G_frames_full = []
frame_count = 0

# Set up the process folder
num_files = len(data_files)
output_dir = os.path.dirname(data_dir)
output_file = os.path.join(output_dir, 'RawDemixed')
sep_ind_file = os.path.join(output_dir, 'SeqIndex')
global_ind_file = os.path.join(output_dir, 'globalIndex')
output_empty = os.path.join(output_dir, 'Empty')


[full_frames, frame_count, global_count] = insertEmpty(data_files, data_dir, frame_count, global_count)

os.makedirs(output_empty, exist_ok=True)

# Save the full data into .mat format
print(type(full_frames))
print(len(full_frames))
print(type(full_frames[0]))
print(np.array(full_frames).shape)
indicatorname = os.path.join(output_empty, 'indicator.mat')
scaleRatio = 256/(x_end-x_start)
x1scale = round((x1-x_start)*scaleRatio)
x2scale = round((x2-x_start)*scaleRatio)
y1scale = round((y1-y_start)*scaleRatio)
y2scale = round((y2-y_start)*scaleRatio)
savemat(indicatorname, {'indicator': full_frames[:, y1scale:y2scale, x1scale:x2scale]})


    

'''
# Process the chunk of files and save the result
[frame_count, seperateIndB, seperateIndG, globalB, globalG, global_count] = process_files(data_files, data_dir, output_file, B_frames_full, G_frames_full, seperateIndB, seperateIndG, frame_count, globalB, globalG, global_count)

saveseqindex(sep_ind_file, seperateIndB, '1')
saveseqindex(sep_ind_file, seperateIndG, '2')

saveglobalindex(global_ind_file, globalB, '1')
saveglobalindex(global_ind_file, globalG, '2')

print(f"Processed and saved {output_file}")
'''



'''
if selected_frames:
    # Convert selected frames to a numpy array of type float32
    selected_frames_np = np.array(selected_frames, dtype=np.float32)

    # Reorder the array to have shape (256, 256, num_frames)
    # selected_frames_np.shape will be (num_frames, 256, 256)
    # We need to reshape it to (256, 256, num_frames)
    frames_final = np.transpose(selected_frames_np, (1, 2, 0))

    # Create the output filename for .mat file
    output_mat_filename = os.path.join(output_folder_path, f"processed_{filename.split('.')[0]}.mat")

    # Save the frames to a .mat file (256, 256, num_frames)
    sio.savemat(output_mat_filename, {'frames': frames_final})

    print(f"Processed frames saved as {output_mat_filename}")
else:
    print(f"No frames selected for {filename}. Skipping video creation.")


#data = np.empty((output_shape[0],output_shape[1] ,frames_per_file*(end_file-start_file)),np.float32)*np.nan
data = np.empty((output_shape[0],output_shape[1] ,frames_per_file*len(data_files)),np.float32)*np.nan
frame_it = 0
for df in tqdm(data_files[start_file:end_file]):
    cap = cv2.VideoCapture(df)
    for it in tqdm( range(int(cap.get(cv2.CAP_PROP_FRAME_COUNT))),leave=False):
        ret, frame = cap.read()
        data[:,:,frame_it] = skimage.transform.resize(frame[y_start:y_end,x_start:x_end,0].astype(np.float32),output_shape)
        frame_it+=1


bs = data[:,:,::2] #backscatter
blue = data[:,:,1::2] #fluorescence
for i in range(256):
    savemat(f"/vast/palmer/scratch/higley/hd362/2025_02_03_16_47_18/RawDemixed/RawDemixedCH1Col{i}.mat",{"column": blue[:,i,:]})
    savemat(f"/vast/palmer/scratch/higley/hd362/2025_02_03_16_47_18/RawDemixed/RawDemixedCH2Col{i}.mat",{"column": bs[:,i,:]})


#savemat(r'W:\Hao\Miniscope Data\Miniscope\2024_08_01\12_28_29\mobile_data.mat',{'mobile_data' : data})

'''
