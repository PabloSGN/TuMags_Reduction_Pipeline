# import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# from datetime import datetime
from scipy.ndimage import affine_transform

# # === Step 1: Read CSV ===
# df = pd.read_csv("pinholes.csv")  # <-- make sure the file path is correct

# # Convert day, hour, min to a single timestamp value (minutes since start)
# df['timestamp'] = df['day'] * 1440 + df['hour'] * 60 + df['min']

# # === Step 2: Interpolation Function ===
# def interpolate_filter(timestamp_query):
#     timestamps = df['timestamp'].values.astype(float)
    
#     results = {}
#     for filter_name in ['filter1', 'filter2', 'filter3']:
#         for axis in ['rot', 'x', 'y']:
#             col = f"{filter_name}_{axis}"
#             y = df[col].values.astype(float)

#             # Use cubic interpolation for smooth curve (can mimic sine shape)
#             interp_func = interp1d(timestamps, y, kind='cubic', fill_value="extrapolate")
#             results[f"{filter_name}_{axis}"] = interp_func(timestamp_query)
    
#     return results

# # === Step 3: Plotting ===
# def plot_with_interpolation(timestamp_query):
#     timestamps = df['timestamp'].values
#     plt.figure(figsize=(15, 8))
    
#     for i, filter_name in enumerate(['filter1', 'filter2', 'filter3'], 1):
#         for axis in ['rot', 'x', 'y']:
#             col = f"{filter_name}_{axis}"
#             y = df[col].values.astype(float)
            
#             # Interpolation
#             interp_func = interp1d(timestamps, y, kind='cubic', fill_value="extrapolate")
#             y_interp = interp_func(timestamp_query)
            
#             # Plot
#             plt.subplot(3, 3, (i - 1) * 3 + ['rot', 'x', 'y'].index(axis) + 1)
#             plt.plot(timestamps, y, 'o', label='Original Data')
#             plt.plot(timestamp_query, y_interp, 'rx', label='Interpolated', markersize=10)
#             plt.title(f'{filter_name} - {axis}')
#             plt.xlabel('Timestamp')
#             plt.ylabel(axis)
#             plt.legend()
    
#     plt.tight_layout()
#     plt.show()

# # # === Step 4: Example Usage ===
# # # Example input (day 13, hour 9, min 30)
# # day, hour, minute = 13, 9, 30
# # query_timestamp = day * 1440 + hour * 60 + minute

# # result = interpolate_filter(query_timestamp)
# # print(f"Interpolated values for Day {day}, Hour {hour}, Min {minute}:\n")
# # for k, v in result.items():
# #     print(f"{k}: {v:.6f}")

# # plot_with_interpolation(query_timestamp)


def interpolate_filter(df, timestamp_query):
    timestamps = df['timestamp'].values.astype(float)
    results = {}

    # Diccionario con valores por defecto por parámetro
    default_values = {
        "rotation_angle_deg": 0.0,
        "translation_y": 0.0,
        "translation_x": 0.0,
        "center_y": 0.0,
        "center_x": 0.0,
        "scale_x": 1.0,
        "scale_y": 1.0,
        "shear_x": 0.0,
        "shear_y": 0.0
    }

    for key in default_values:
        if key in df.columns:
            y = df[key].values.astype(float)
            interp_func = interp1d(timestamps, y, kind='cubic', fill_value="extrapolate")
            results[key] = float(interp_func(timestamp_query))
        else:
            results[key] = default_values[key]

    return results

def apply_transform(image, angle_deg, t, center,
                    scale_x=1.0, scale_y=1.0, shear_x=0.0, shear_y=0.0):
    angle_rad = np.radians(angle_deg)
    cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    Distortion = np.array([[scale_x, shear_x], [shear_y, scale_y]])
    A = Distortion @ R
    A_inv = np.linalg.inv(A)
    offset = center - A_inv @ (center + t)

    return affine_transform(
        image, matrix=A_inv, offset=offset, order=3, mode='nearest', cval=0.0
    )