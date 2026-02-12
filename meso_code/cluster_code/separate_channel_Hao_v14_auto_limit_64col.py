# -*- coding: utf-8 -*-
"""
Project: Hemodynamics correction pipeline

Created: Thu Aug  1 14:46:06 2024
Modified: Mon May 18 11:10 2025

Authors:
    Hao Dong

Description:
    Processes mesoscope imaging data by:
    - Detecting frame drops from timestamps and inserting empty frames
    - Cropping/resizing frames and computing channel indicators
    - Inserting/correcting frame assignments based on indicator thresholds
    - Separating G/B channels and balancing frame counts
    - Saving column-wise .mat outputs for downstream processing
"""

# ======================================================================
# Imports
# ======================================================================

# Standard library
import os
import sys
import glob
import json
from collections import defaultdict

# Third-party
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from natsort import natsorted
from scipy.io import savemat, loadmat
from scipy.signal import butter, sosfiltfilt  # (kept for parity with original)
import skimage  # (kept for parity with original)
import keyboard  # (kept for parity with original)

# ======================================================================
# Global parameters (defaults)
# ======================================================================

startIndex = 0
endIndex = 0

# NOTE: These globals are used inside functions; they are set in main before use.
main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training'
mouse_folder = glob.glob(os.path.join(main_folder, '*Mouse*'))

# I/O + processing globals (populated later in main)
output_shape = (256, 256)
fDropList = []
fList = []

# Cropping/ROI globals (populated in main from Excel)
x_start = x_end = y_start = y_end = 0
x1 = x2 = y1 = y2 = 0
uppbdy = lowbdy = 0

# ======================================================================
# Functions
# ======================================================================

def findFrameDrop(meso_folder):
    """
    Identify frame drops using timeStamps*.csv within the meso folder.
    Saves results to frameDrop.mat (list of [preFrame, count] pairs).
    Returns:
        fDropData : list of [preFrame, count]
    """
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
        if numDiff == 0:
            frameCount += 1
            frameSave.append(frameCount - 1)
        else:
            for _ in range(numDiff):
                frameCount += 1
                dropFrame.append(frameCount - 1)
            frameCount += 1
            frameSave.append(frameCount - 1)

    # Create the frameCheck DataFrame
    frameCheck = pd.DataFrame(np.ones((frameCount, 2), dtype=int))
    frameCheck.iloc[dropFrame, :] = 0
    frameCheck.iloc[frameSave, 1] = np.arange(1, len(frameSave) + 1)

    fDropData = []
    for i in range(len(dropFrame)):
        drop_idx = dropFrame[i]
        pre_idx = drop_idx - 1
        preFrame = frameCheck.iloc[pre_idx, 1]
        if preFrame:
            fDropData.append([preFrame, 1])
        else:
            if fDropData:
                fDropData[-1][1] += 1

    # Filter empty/falsey entries
    fDropData = [entry for entry in fDropData if entry[0]]

    # Save to .mat
    fDropFolder = os.path.join(os.path.dirname(meso_folder[0]), 'frameDrop.mat')
    savemat(fDropFolder, {'frameDrop': fDropData})

    return fDropData


def insertEmpty(data_files, data_dir, frame_count, global_count):
    """
    Insert empty frames at detected drop locations while reading & preprocessing video frames.
    Uses globals:
        output_shape, x_start, x_end, y_start, y_end, fDropList, fList
    Returns:
        full_frames : np.array of shape (N, 256, 256)
        frame_count, global_count : updated counters
    """
    full_frames = []
    for i, filename in enumerate(data_files):
        video_path = os.path.join(data_dir, filename)

        # Open the video
        cap = cv2.VideoCapture(video_path)
        frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
        frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

        empty_frame = np.zeros(output_shape, dtype=np.uint8)

        if not cap.isOpened():
            print(f"Error: Could not open video {filename}.")
            continue

        # Skip first 4 frames only for the first video
        if frame_count == 0:
            for _ in range(4):
                cap.read()
            frame_count = 4
            global_count = 4

        # Read frames and process
        while True:
            ret, frame = cap.read()
            if not ret:
                break

            frame_count += 1
            global_count += 1

            # Crop to ROI (x_start:x_end, y_start:y_end)
            cropped_frame = frame[y_start:y_end, x_start:x_end]

            # Grayscale and resize to 256x256
            grayscale_frame = cv2.cvtColor(cropped_frame, cv2.COLOR_BGR2GRAY).astype(np.uint8)
            resized_frame = cv2.resize(grayscale_frame, output_shape)

            full_frames.append(resized_frame)

            # Insert empty frames at drop locations
            if frame_count in fList:
                for row in fDropList:
                    if row[0] == frame_count:
                        num_drop = row[1]
                        for _ in range(num_drop):
                            global_count += 1
                            full_frames.append(empty_frame)

        cap.release()

    return np.array(full_frames), frame_count, global_count


def insertCorrect(full_frames, meanvalue, uppbdy, lowbdy):
    """
    Insert/correct frames in-place guided by meanvalue thresholds.
    """
    for n in range(len(meanvalue) - 1):
        # neigh_diff is calculated but unused (preserved to keep logic identical)
        neigh_diff = meanvalue[n + 1] - meanvalue[n]
        if abs(meanvalue[n + 1]) < lowbdy:
            if meanvalue[n] >= uppbdy:
                [full_frames, meanvalue] = find_near_b(full_frames, meanvalue, n, uppbdy, lowbdy)
            elif meanvalue[n] >= lowbdy:
                [full_frames, meanvalue] = find_near_g(full_frames, meanvalue, n, uppbdy, lowbdy)
    return full_frames, meanvalue


def find_near_b(full_frames, meanvalue, dindex, uppbdy, lowbdy):
    """
    Backtrack to find nearest 'blue' (between lowbdy and uppbdy) to copy forward.
    """
    for m in range(1, dindex + 1):
        if lowbdy <= meanvalue[dindex - m] < uppbdy:
            meanvalue[dindex + 1] = meanvalue[dindex - m]
            full_frames[dindex + 1] = full_frames[dindex - m]
            break
    return full_frames, meanvalue


def find_near_g(full_frames, meanvalue, dindex, uppbdy, lowbdy):
    """
    Backtrack to find nearest 'green' (>= uppbdy) to copy forward.
    """
    for m in range(1, dindex + 1):
        if meanvalue[dindex - m] >= uppbdy:
            meanvalue[dindex + 1] = meanvalue[dindex - m]
            full_frames[dindex + 1] = full_frames[dindex - m]
            break
    return full_frames, meanvalue


def correctpair(full_frames, meanvalue, uppbdy, lowbdy):
    """
    Pairwise correction across frames (even/odd pairs).
    """
    for n in range(int(len(meanvalue) / 2)):
        frame1value = meanvalue[n * 2]
        frame2value = meanvalue[n * 2 + 1]
        if frame1value > uppbdy and frame2value > uppbdy:
            [full_frames, meanvalue] = find_near_b(full_frames, meanvalue, n * 2, uppbdy, lowbdy)
        elif lowbdy < frame1value < uppbdy and lowbdy < frame2value < uppbdy:
            [full_frames, meanvalue] = find_near_g(full_frames, meanvalue, n * 2, uppbdy, lowbdy)
    return full_frames, meanvalue


def chnlDiff(meanvalue, uppbdy, lowbdy):
    """
    Count difference between frames above uppbdy (G) and between (lowbdy, uppbdy) (B).
    """
    g_chan = [i for i, val in enumerate(meanvalue) if val > uppbdy]
    b_chan = [i for i, val in enumerate(meanvalue) if lowbdy < val < uppbdy]
    gb_diff = len(g_chan) - len(b_chan)
    return gb_diff


def sepGBChnl(full_frames, meanvalue, uppbdy):
    """
    Separate frames into G (meanvalue > uppbdy) and B (lowbdy < meanvalue < uppbdy).
    Uses global lowbdy to preserve original behavior.
    """
    G_frames = []
    B_frames = []
    for n in range(len(meanvalue)):
        if meanvalue[n] > uppbdy:
            G_frames.append(full_frames[n, :, :])
        elif lowbdy < meanvalue[n] < uppbdy:
            B_frames.append(full_frames[n, :, :])
    return G_frames, B_frames


#def savecolumnmat(output_file, frames_full, colInd, chInd):
#    """
#    Save column-wise slices to .mat files, appending if file exists.
#    """
#    slice_data = frames_full[:, :, colInd].T
#    ch1_filename = f"{output_file}/RawDemixedCH{chInd}Col{colInd}.mat"
#
#    if os.path.exists(ch1_filename):
#        existing_data_ch1 = loadmat(ch1_filename)
#        existing_column_ch1 = existing_data_ch1.get("column", np.array([]))
#        updated_column_ch1 = np.concatenate([existing_column_ch1, slice_data], axis=1)
#    else:
#        updated_column_ch1 = slice_data
#
#    savemat(ch1_filename, {"column": updated_column_ch1})
    
def savecolumnmat(output_file, frames_full, colstart, colInd, chInd):
    # Prepare the slice data for saving (transpose if necessary)
    slice_data = frames_full[:,:,colInd].T
    colname = colstart/4+1

    channel = 'Blue' if chInd == '1' else 'UV'
    ch1_filename = f"{output_file}/Data{channel}Col{int(colname)}.mat"


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


def saveseqindex(sep_ind_file, seqIndex, chInd):
    """
    Save sequence index to .mat.
    """
    os.makedirs(sep_ind_file, exist_ok=True)
    seq_filename = f"{sep_ind_file}/sepIndex{chInd}.mat"
    savemat(seq_filename, {"seqIndex": seqIndex})


def saveglobalindex(global_ind_file, globalIndex, chInd):
    """
    Save global index to .mat.
    """
    os.makedirs(global_ind_file, exist_ok=True)
    seq_filename = f"{global_ind_file}/globalIndex{chInd}.mat"
    savemat(seq_filename, {"globalIndex": globalIndex})


# ======================================================================
# Main workflow
# ======================================================================

if __name__ == "__main__":
    # --------------------------------------------------------------
    # Inputs & discovery
    # --------------------------------------------------------------
    mouse_input = sys.argv[1]
    normalization_method = sys.argv[2]
    mouse_folder = glob.glob(os.path.join(main_folder, '*Mouse*'))
    mouse_select = [folder for folder in mouse_folder if mouse_input in folder]

    # Find the meso data folder
    meso_folder = glob.glob(os.path.join(mouse_select[0], '*', 'My_V4_Miniscope'), recursive=True)
    print(meso_folder)

    # --------------------------------------------------------------
    # Frame drop detection
    # --------------------------------------------------------------
    fDropData = findFrameDrop(meso_folder)
    fDropList = []
    fList = []
    for item in fDropData:
        frame_locs = item[0]  # e.g., 538
        count = item[1]       # e.g., 1
        fDropList.append([frame_locs, count])
        fList.append(frame_locs)
        print(f"Frame: {frame_locs}, Count: {count}")

    # --------------------------------------------------------------
    # Mouse ID (from path) and data files
    # --------------------------------------------------------------
    parts = os.path.normpath(meso_folder[0]).split(os.sep)
    try:
        mouseID = int(parts[7][8:9])
    except ValueError:
        mouseID = parts[7][8:9]
    print(mouseID)

    data_dir = meso_folder[0]
    data_files = natsorted(glob.glob(os.path.join(data_dir, '*.avi')))
    print(data_dir)

    # --------------------------------------------------------------
    # Meta info (frameRate, ROI, framesPerFile)
    # --------------------------------------------------------------
    with open(os.path.join(data_dir, 'metaData.json')) as f:
        meta_data = json.load(f)

    frame_rate = meta_data['frameRate']
    roi = meta_data['ROI']
    frames_per_file = meta_data['framesPerFile']

    output_shape = (256, 256)

    # --------------------------------------------------------------
    # Load per-mouse ROI/thresholds from Excel
    # --------------------------------------------------------------
    mouseInfofile = glob.glob(os.path.join(mouse_select[0], 'mouseMesoInfo*'))
    print(mouseInfofile)
    mouseInfo = pd.read_excel(mouseInfofile[0])
    # Get the row matching the mouse_input
    row = mouseInfo.loc[mouseInfo["MouseID"] == mouse_input]

    if row.empty:
        raise ValueError(f"MouseID {mouse_input} not found in mouseInfo file.")

    # If only one match is expected, take the first row
    row = row.iloc[0]

    uppbdy  = int(row['uplimit'])
    lowbdy  = int(row['lowlimit'])

    y_start = int(row['y_start'])
    y_end   = int(row['y_end'])
    x_start = int(row['x_start'])
    x_end   = int(row['x_end'])

    x1 = int(row['x1'])
    x2 = int(row['x2'])
    y1 = int(row['y1'])
    y2 = int(row['y2'])

    print(x1, x2, y1, y2)
    '''
    row_num = mouseInfo.index[mouseInfo["MouseID"] == mouse_input].tolist()

    uppbdy = int(mouseInfo.loc[row_num].uplimit)
    lowbdy = int(mouseInfo.loc[row_num].lowlimit)

    y_start = int(mouseInfo.loc[row_num].y_start)
    y_end = int(mouseInfo.loc[row_num].y_end)
    x_start = int(mouseInfo.loc[row_num].x_start)
    x_end = int(mouseInfo.loc[row_num].x_end)

    x1 = int(mouseInfo.loc[row_num].x1)
    x2 = int(mouseInfo.loc[row_num].x2)
    y1 = int(mouseInfo.loc[row_num].y1)
    y2 = int(mouseInfo.loc[row_num].y2)
    print(x1, x2, y1, y2)
    '''
    # --------------------------------------------------------------
    # Output folders and counters
    # --------------------------------------------------------------
    globalB = []
    globalG = []
    global_count = startIndex * 1000
    seperateIndB = []
    seperateIndG = []
    B_frames_full = []
    G_frames_full = []
    frame_count = startIndex * 1000

    num_files = len(data_files)
    output_folder = os.path.dirname(data_dir)
    output_dir = os.path.join(output_folder, normalization_method)
    output_file = os.path.join(output_dir, 'aligned')
    sep_ind_file = os.path.join(output_dir, 'SeqIndex')
    global_ind_file = os.path.join(output_dir, 'globalIndex')
    output_empty = os.path.join(output_dir, 'Empty')
    os.makedirs(output_empty, exist_ok=True)

    # --------------------------------------------------------------
    # Read frames + insert empties at drop locations
    # --------------------------------------------------------------
    if endIndex == 0:
        endIndex = len(data_files)

    full_frames, frame_count, global_count = insertEmpty(
        data_files[startIndex:endIndex], data_dir, frame_count, global_count
    )

    # --------------------------------------------------------------
    # Indicator window (scaled to 256x256) for channel separation
    # --------------------------------------------------------------
    scaleRatio = 256 / (x_end - x_start)
    x1scale = round((x1 - x_start) * scaleRatio)
    x2scale = round((x2 - x_start) * scaleRatio)
    y1scale = round((y1 - y_start) * scaleRatio)
    y2scale = round((y2 - y_start) * scaleRatio)

    indicator = full_frames[:, y1scale:y2scale, x1scale:x2scale]
    meanvalue = np.mean(indicator, axis=(1, 2))

    savemat(os.path.join(output_empty, 'indicatorraw.mat'), {'indicator': indicator})
    savemat(os.path.join(output_empty, 'meanvalueraw.mat'), {'meanvalue': meanvalue})

    # --------------------------------------------------------------
    # Insertion (drop) correction
    # --------------------------------------------------------------
    gb_diff0 = chnlDiff(meanvalue, uppbdy, lowbdy)
    print(f"The difference between two channels(before insertion) is: {gb_diff0}")

    full_frames, meanvalue = insertCorrect(full_frames, meanvalue, uppbdy, lowbdy)
    gb_diff1 = chnlDiff(meanvalue, uppbdy, lowbdy)
    print(f"The difference between two channels(after drop insertion) is: {gb_diff1}")

    indicator = full_frames[:, y1scale:y2scale, x1scale:x2scale]
    meanvalue = np.mean(indicator, axis=(1, 2))
    savemat(os.path.join(output_empty, 'indicatorinsert.mat'), {'indicator': indicator})
    savemat(os.path.join(output_empty, 'meanvalueinsert.mat'), {'meanvalue': meanvalue})

    # --------------------------------------------------------------
    # Pair correction
    # --------------------------------------------------------------
    full_frames, meanvalue = correctpair(full_frames, meanvalue, uppbdy, lowbdy)
    gb_diff2 = chnlDiff(meanvalue, uppbdy, lowbdy)
    print(f"The difference between two channels(after pair correction) is: {gb_diff2}")

    indicator = full_frames[:, y1scale:y2scale, x1scale:x2scale]
    meanvalue = np.mean(indicator, axis=(1, 2))
    savemat(os.path.join(output_empty, 'indicatorpair.mat'), {'indicator': indicator})
    savemat(os.path.join(output_empty, 'meanvaluepair.mat'), {'meanvalue': meanvalue})

    # --------------------------------------------------------------
    # Separate channels and balance counts
    # --------------------------------------------------------------
    G_frames, B_frames = sepGBChnl(full_frames, meanvalue, uppbdy)

    frameDiff = abs(len(B_frames) - len(G_frames))
    print(f"The length of G_frames is: {len(G_frames)}, the length of B_frames is: {len(B_frames)}")
    print(f"The difference between two channels is: {frameDiff}")

    if len(B_frames) > len(G_frames):
        for _ in range(frameDiff):
            G_frames.append(G_frames[-1])
    else:
        for _ in range(frameDiff):
            B_frames.append(B_frames[-1])

    frameDiff = abs(len(B_frames) - len(G_frames))
    print(f"The difference after final insertion is: {frameDiff}")

    G_frames = np.array(G_frames).astype(np.float64)
    B_frames = np.array(B_frames).astype(np.float64)

    os.makedirs(output_file, exist_ok=True)

    # --------------------------------------------------------------
    # Save per-column .mat files
    # --------------------------------------------------------------
    #for colInd in range(256):
    #    savecolumnmat(output_file, B_frames, colInd, '1')
    #    savecolumnmat(output_file, G_frames, colInd, '2')
    
    for colstart in range(0, 256, 4):
        colInd = list(range(colstart, colstart + 4)) 
        # Prepare the slice data for saving (transpose if necessary)
        savecolumnmat(output_file, B_frames, colstart, colInd, '1')
        savecolumnmat(output_file, G_frames, colstart, colInd, '2')

    # --------------------------------------------------------------
    # Save final indicator and meanvalue snapshots
    # --------------------------------------------------------------
    print(type(full_frames))
    print(len(full_frames))
    print(type(full_frames[0]))
    print(np.array(full_frames).shape)

    savemat(os.path.join(output_empty, 'indicator.mat'),
            {'indicator': full_frames[:, y1scale:y2scale, x1scale:x2scale]})
    savemat(os.path.join(output_empty, 'meanvalue.mat'),
            {'meanvalue': meanvalue})

    print("\n[INFO] Pipeline finished successfully.")
