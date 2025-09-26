# -*- coding: utf-8 -*-
"""
Project: Hemodynamics correction pipeline (pre-separated channels)

Description:
    Processes pre-separated G/B mesoscope AVI videos:
    - Crop frames to ROI (x_start:x_end, y_start:y_end)
    - Resize frames to 256x256
    - Save per-column .mat outputs
"""

import os
import sys
import glob
import cv2
import numpy as np
from scipy.io import savemat
import pandas as pd

# ======================================================================
# Functions
# ======================================================================

def savecolumnmat(output_file, frames_full, colInd, chInd):
    slice_data = frames_full[:, :, colInd].T
    ch_filename = os.path.join(output_file, f"RawDemixedCH{chInd}Col{colInd}.mat")

    if os.path.exists(ch_filename):
        from scipy.io import loadmat
        existing_data = loadmat(ch_filename)
        existing_column = existing_data.get("column", np.array([]))
        updated_column = np.concatenate([existing_column, slice_data], axis=1)
    else:
        updated_column = slice_data

    savemat(ch_filename, {"column": updated_column})


def read_and_process_channel(video_path, output_file, chInd, x_start, x_end, y_start, y_end, output_shape=(256, 256)):
    """
    Read AVI, crop to ROI, resize to output_shape, and save per-column .mat files.
    """
    os.makedirs(output_file, exist_ok=True)
    cap = cv2.VideoCapture(video_path)
    frames = []

    while True:
        ret, frame = cap.read()
        if not ret:
            break
        if len(frame.shape) == 3:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        # Crop to ROI
        frame = frame[y_start:y_end, x_start:x_end]
        # Resize to 256x256
        frame_resized = cv2.resize(frame, output_shape)
        frames.append(frame_resized)

    cap.release()
    frames = np.array(frames).astype(np.float64)

    print(f"[INFO] Finished saving channel {chInd}, {len(frames)} frames.")
    return frames


# ======================================================================
# Main workflow
# ======================================================================

if __name__ == "__main__":
    mouse_input = sys.argv[1]  # e.g., "Mouse1"
    main_folder = '/vast/palmer/scratch/higley/hd362/HD_Mouse_Training'

    mouse_folder = glob.glob(os.path.join(main_folder, '*Mouse*'))
    mouse_select = [folder for folder in mouse_folder if mouse_input in folder]

    if not mouse_select:
        raise ValueError(f"No mouse folder found for {mouse_input}")

    data_dir = glob.glob(os.path.join(mouse_select[0], '*', 'My_V4_Miniscope'))[0]

    # Pre-separated AVI files
    G_video_path = os.path.join(data_dir, 'All_GrnGrn2.avi')
    B_video_path = os.path.join(data_dir, 'All_RedGrn1.avi')

    output_dir = os.path.dirname(data_dir)
    output_file = os.path.join(output_dir, 'RawDemixed')
    os.makedirs(output_file, exist_ok=True)

    # Load ROI from Excel
    mouseInfofile = glob.glob(os.path.join(mouse_select[0], 'mouseMesoInfo*'))[0]
    mouseInfo = pd.read_excel(mouseInfofile)
    row_num = mouseInfo.index[mouseInfo["MouseID"] == mouse_input].tolist()[0]

    # Use .at to avoid FutureWarning
    x_start = int(mouseInfo.at[row_num, 'x_start'])
    x_end = int(mouseInfo.at[row_num, 'x_end'])
    y_start = int(mouseInfo.at[row_num, 'y_start'])
    y_end = int(mouseInfo.at[row_num, 'y_end'])

    output_shape = (256, 256)

    # Read and process channels
    G_frames = read_and_process_channel(G_video_path, output_file, chInd='2',
                                        x_start=x_start, x_end=x_end,
                                        y_start=y_start, y_end=y_end,
                                        output_shape=output_shape)

    B_frames = read_and_process_channel(B_video_path, output_file, chInd='1',
                                        x_start=x_start, x_end=x_end,
                                        y_start=y_start, y_end=y_end,
                                        output_shape=output_shape)

    # --------------------------------------------------------------
    # Truncate longer video immediately after reading
    # --------------------------------------------------------------
    n_frames = min(len(G_frames), len(B_frames))
    if len(G_frames) != len(B_frames):
        print(f"[INFO] Truncating frames to {n_frames} to match both channels.")
        G_frames = G_frames[:n_frames]
        B_frames = B_frames[:n_frames]

    # --------------------------------------------------------------
    # Save per-column .mat files
    # --------------------------------------------------------------
    for colInd in range(output_shape[1]):
        savecolumnmat(output_file, B_frames, colInd, '1')
        savecolumnmat(output_file, G_frames, colInd, '2')

    print("\n[INFO] Pipeline finished successfully.")

