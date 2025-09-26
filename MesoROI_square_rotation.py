import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from skimage.transform import SimilarityTransform, warp

video_path = r"C:\Users\hao\Downloads\2.avi"
parcellated_path = r'C:\Users\hao\Downloads\parcels_updated12522.mat'
xlsx_path = r"C:\Users\hao\Downloads\mouseMesoInfo.xlsx"

# -------------------- 1. USER INPUT --------------------
mouse_id = input("Enter Mouse ID: ").strip()

# --- Extract middle part (e.g. "MouseF") ---
# Assumes ID looks like "SP_MouseF_0820_pre"
parts = mouse_id.split("_")
if len(parts) >= 2:
    core_id = parts[1]   # "MouseF"
else:
    raise ValueError("Mouse ID format not recognized (expected at least 2 parts separated by '_').")
tform_save_path = f'C:/Users/hao/Downloads/tform{core_id}.mat'

# -------------------- 2. LOAD VIDEO --------------------
cap = cv2.VideoCapture(video_path)
if not cap.isOpened():
    raise IOError("❌ Cannot open video")

ret, first_frame = cap.read()
ret2, second_frame = cap.read()
if not ret or not ret2:
    raise IOError("❌ Cannot read frames")

added_frame = cv2.add(first_frame, second_frame)
fixed_size = 420
roi1_pos = [50, 50]
dragging = False
clone = added_frame.copy()

# -------------------- 3. INTERACTIVE FIXED ROI --------------------
def draw_roi(event, x, y, flags, param):
    global roi1_pos, dragging
    if event == cv2.EVENT_LBUTTONDOWN:
        dragging = True
    elif event == cv2.EVENT_MOUSEMOVE and dragging:
        roi1_pos = [x, y]
    elif event == cv2.EVENT_LBUTTONUP:
        dragging = False

cv2.namedWindow("Move Fixed ROI (Press Enter)")
cv2.setMouseCallback("Move Fixed ROI (Press Enter)", draw_roi)

while True:
    frame = clone.copy()
    x, y = np.clip(roi1_pos, [0, 0], [added_frame.shape[1]-fixed_size, added_frame.shape[0]-fixed_size])
    cv2.rectangle(frame, (x, y), (x+fixed_size, y+fixed_size), (0, 255, 0), 2)
    cv2.imshow("Move Fixed ROI (Press Enter)", frame)
    if cv2.waitKey(20) == 13:  # Enter
        break

cv2.destroyWindow("Move Fixed ROI (Press Enter)")
roi1_pos = [x, y]
roi1_corners = [(x, y), (x+fixed_size, y), (x+fixed_size, y+fixed_size), (x, y+fixed_size)]
print("ROI 1 Corners:", roi1_corners)

# -------------------- 4. SELECT ROI 2 --------------------
def select_square_roi(window_name, frame):
    while True:
        x, y, w, h = map(int, cv2.selectROI(window_name, frame, fromCenter=False))
        s = min(w, h)
        cv2.destroyWindow(window_name)
        return (x, y, s, s)

x2, y2, w2, h2 = select_square_roi("Select ROI 2", added_frame)
roi2_corners = [(x2, y2), (x2+w2, y2), (x2+w2, y2+h2), (x2, y2+h2)]
print("ROI 2 Corners:", roi2_corners)

# -------------------- 5. EXTRACT ROI2 SIGNAL --------------------
roi2_means = []
cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
while True:
    ret, frame = cap.read()
    if not ret: break
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    roi_crop = gray[y2:y2+h2, x2:x2+w2]
    roi2_means.append(np.mean(roi_crop))
cap.release()

plt.plot(roi2_means)
plt.title(f"Mouse {mouse_id} - ROI 2 Average Intensity")
plt.xlabel("Frame")
plt.ylabel("Mean Intensity")
plt.grid(True)
plt.show()

# -------------------- 6. ADD UP and LOW LIMIT --------------------

uplimit = float(input("Enter uplimit value: ").strip())
lowlimit = float(input("Enter lowlimit value: ").strip())

# -------------------- 7. SAVE ROI COORDINATES --------------------
row = {
    "MouseID": mouse_id,
    "x_start": roi1_pos[0],
    "x_end": roi1_pos[0] + fixed_size,
    "y_start": roi1_pos[1],
    "y_end": roi1_pos[1] + fixed_size,
    "x1": x2,
    "x2": x2 + w2,
    "y1": y2,
    "y2": y2 + h2,
    "uplimit": uplimit,
    "lowlimit": lowlimit,
}

if os.path.exists(xlsx_path):
    df = pd.read_excel(xlsx_path)
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
else:
    df = pd.DataFrame([row])
df.to_excel(xlsx_path, index=False)
print(f"✅ Mouse info saved to: {xlsx_path}")

# -------------------- 8. GET FIRST CROPPED + RESIZED FRAME --------------------
cap = cv2.VideoCapture(video_path)
if not cap.isOpened():
    raise IOError("❌ Cannot reopen video")

ret, frame = cap.read()
if not ret:
    raise IOError("❌ Cannot read first frame for cropping")

x, y = roi1_pos
x = min(max(0, x), frame.shape[1] - fixed_size)
y = min(max(0, y), frame.shape[0] - fixed_size)

cropped = frame[y:y+fixed_size, x:x+fixed_size]
image_to_align = cv2.resize(cropped, (256, 256))

cap.release()
print("✅ Got first cropped and resized frame for alignment.")

# -------------------- 9. LOAD PARCELLATION --------------------
data = loadmat(parcellated_path)
combined_parcels = data['allen_parcels'][0, 0]['CombinedParcells']
mask = (combined_parcels > 0).astype(np.uint8) * 255

# -------------------- 10. ALIGN GUI (Interactive Points) --------------------

image_gray_raw = cv2.cvtColor(image_to_align, cv2.COLOR_BGR2GRAY)
kernel_sharpen = np.array([
    [0, -1, 0],
    [-1, 5, -1],
    [0, -1, 0]
])
# image_gray = cv2.filter2D(image_gray_raw, -1, kernel_sharpen)
image_gray = cv2.equalizeHist(image_gray_raw)
fixed_points = np.array([
    [70.8576, 62.6474],
    [187.8377, 62.6474],
    [129.0651, 184.9211]
], dtype=np.float32)

moving_points = np.array([
    [64, 64],
    [192, 64],
    [128, 192]
], dtype=np.float32)

drag_index = -1
scale_factor = 2

def draw_points(img, points, color=(0, 0, 255)):
    for pt in points:
        pt_scaled = (pt * scale_factor).astype(int)
        cv2.circle(img, tuple(pt_scaled), 5, color, -1)

def apply_transform(moving_pts, fixed_pts, image):
    t = SimilarityTransform()
    t.estimate(moving_pts, fixed_pts)
    warped = warp(image, inverse_map=t.inverse, output_shape=(256, 256))
    return (warped * 255).astype(np.uint8), t

def mouse_callback(event, x, y, flags, param):
    global drag_index, moving_points
    pos = np.array([x, y]) / scale_factor
    if event == cv2.EVENT_LBUTTONDOWN:
        for i, pt in enumerate(moving_points):
            if np.linalg.norm(pt - pos) < 10:
                drag_index = i
                break
    elif event == cv2.EVENT_MOUSEMOVE and drag_index != -1:
        moving_points[drag_index] = pos
    elif event == cv2.EVENT_LBUTTONUP:
        drag_index = -1

cv2.namedWindow("Align Points GUI")
cv2.setMouseCallback("Align Points GUI", mouse_callback)

while True:
    aligned_img, tform = apply_transform(moving_points, fixed_points, image_gray)
    overlay = cv2.addWeighted(aligned_img, 0.9, mask, 0.1, 0)
    display = cv2.cvtColor(overlay, cv2.COLOR_GRAY2BGR)
    display_large = cv2.resize(display, (256 * scale_factor, 256 * scale_factor))
    draw_points(display_large, moving_points)
    cv2.imshow("Align Points GUI", display_large)
    if cv2.waitKey(20) in [27, ord('q')]:
        break

cv2.destroyAllWindows()

# -------------------- 11. SAVE TRANSFORM TO .mat --------------------
tform_struct = {
    'Dimensionality': 2,
    'Scale': tform.scale,
    'RotationAngle': np.rad2deg(tform.rotation),
    'Translation': np.array(tform.translation).reshape(1, 2),
    'R': tform.params[:2, :2],
    'A': tform.params
}
savemat(tform_save_path, {'tform': tform_struct})
print(f"✅ Saved transformation to: {tform_save_path}")
