import numpy as np
import glob
import matplotlib.pyplot as plt


import config as cf


import field_stop_finder as fsf
from utils import read_Tumag


#iron = np.load("iron.npy")
mag = np.load("magnesium.npy")

"""fig, axs = plt.subplots(2, 8, figsize = (15, 5))


for i in range(8):


    axs[0, i].set_title(f"{cf.om_config['2.06']['lambda_array'][i]} mA")
    axs[0, i].imshow(iron[0, i, 0, 200:1750, 200:1750], cmap = "gray", vmin = 0, vmax = 1.5)
    axs[1, i].imshow(iron[1, i, 0, 200:1750, 200:1750], cmap = "gray", vmin = 0, vmax = 1.5)

    axs[0, i].set_xticks([])
    axs[0, i].set_yticks([])

    axs[1, i].set_xticks([])    
    axs[1, i].set_yticks([])    


axs[0, 0].set_ylabel("cam 1")
axs[1, 0].set_ylabel("cam 2")

plt.tight_layout()
plt.show()"""


"""
fig, axs = plt.subplots(2, 10, figsize = (15, 5))
for i in range(10):
    axs[0, i].set_title(f"{cf.om_config['1']['lambda_array'][i]} mA")
    axs[0, i].imshow(mag[0, i, 0, 200:1750, 200:1750], cmap = "gray", vmin = 0, vmax = 0.8)
    axs[1, i].imshow(mag[1, i, 0, 200:1750, 200:1750], cmap = "gray", vmin = 0, vmax = 0.8)

    axs[0, i].set_xticks([])
    axs[0, i].set_yticks([])

    axs[1, i].set_xticks([])    
    axs[1, i].set_yticks([])    
axs[0, 0].set_ylabel("cam 1")
axs[1, 0].set_ylabel("cam 2")

plt.tight_layout()
plt.show()
"""

pinholes1 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/36"
pinholes2 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/37"
pinholes3 = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/Spectral_calibration/38"
ims_c1 = sorted(glob.glob(f"{pinholes1}/*_0_*"))
ims_c2 = sorted(glob.glob(f"{pinholes1}/*_1_*"))

flats = "/home/pablo/Desktop/TuMAGDATA/flare_set_ordered/2.02/39"


flats_c1 = sorted(glob.glob(f"{flats}/*_0_*"))
flats_c2 = sorted(glob.glob(f"{flats}/*_1_*"))


first_flat, _ = read_Tumag(flats_c1[0])

ff = np.zeros((2, np.shape(first_flat)[0],np.shape(first_flat)[0]))
for img_ind, img_path in enumerate(flats_c1):
    I, _ = read_Tumag(img_path)
    ff[0] += I 

for img_ind, img_path in enumerate(flats_c2):
    I, _ = read_Tumag(img_path)
    I = np.flip(I, axis = -1)
    ff[1] += I 

ff[0] /= len(flats_c1)
ff[1] /= len(flats_c2)

f1, _ = read_Tumag(flats_c1[0])
f2, _ = read_Tumag(flats_c2[1])

f1 = f1 / np.max(f1)
f2 = f2 / np.max(f2)

fsc1, fsc2 = fsf.compute_alignment(flat_cam1=ff[0], flat_cam2=ff[1], pinhole_c1_path=ims_c1[0], pinhole_c2_path=ims_c2[0], plot_flag=True, method = "flats")

