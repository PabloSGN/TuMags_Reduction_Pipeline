
# Built-in 
import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Own Libs
sys.path.append("/home/users/dss/orozco/Tumag/PabloTests")
import config as cf
from utils import read_Tumag
import image_handler as ih

# CONFIG
paths_517 = ih.get_images_paths("D10-5735-5768")
paths_502 = ih.get_images_paths("D10-5769-5802")
paths_506 = ih.get_images_paths("D10-5803-5836")


int_517 = []
volts_517 = []

int_502 = []
volts_502 = []

int_506 = []
volts_506 = []

for img in paths_517:
        
    I, H = ih.read(img)

    if H['cam'] == 0:
        int_517.append(np.mean(I[500:1500, 500:1500]))
        volts_517.append(H["hvps_read_volts"])
            
    print('Images size : ', np.shape(I))

for img in paths_502:
        
    I, H = ih.read(img)

    if H['cam'] == 0:
        int_502.append(np.mean(I[500:1500, 500:1500]))
        volts_502.append(H["hvps_read_volts"])
            
    print('Images size : ', np.shape(I))


for img in paths_506:
        
    I, H = ih.read(img)

    if H['cam'] == 0:
        int_506.append(np.mean(I[500:1500, 500:1500]))
        volts_506.append(H["hvps_read_volts"])
            
    print('Images size : ', np.shape(I))

int_517  = np.array(int_517)
volts_517 = np.array(volts_517)

int_502  = np.array(int_502)
volts_502 = np.array(volts_502)

int_506  = np.array(int_506)
volts_506 = np.array(volts_506)
    
# Compute mean value

fig, axs = plt.subplots(3, 1, figsize = (9, 9))
fig.suptitle("Voltage scan - Disk Center")
#axs[0].set_title('Voltage Scan - Disk Center - 517 nm')
axs[0].plot(volts_517, int_517 / H["nAcc"], color = 'crimson', lw =3, label = "517 nm Pre-filter")
axs[0].scatter(volts_517, int_517 / H["nAcc"], marker = 'x', c = 'k', s = 60)
axs[0].set_xticks(volts_517)

axs[1].plot(volts_502, int_502 / H["nAcc"], color = 'crimson', lw =3, label = "525.02 nm Pre-filter")
axs[1].scatter(volts_502, int_502 / H["nAcc"], marker = 'x', c = 'k', s = 60)
axs[1].set_xticks(volts_502)

axs[2].plot(volts_506, int_506 / H["nAcc"], color = 'crimson', lw =3, label = "525.06 nm Pre-filter")
axs[2].scatter(volts_506, int_506 / H["nAcc"], marker = 'x', c = 'k', s = 60)
axs[2].set_xticks(volts_506)

for i in range(3):
    axs[i].set_ylabel('Intensity [per accumulation]')
    axs[i].set_xlabel('Voltage [V]')
    axs[i].legend(edgecolor = 'k')
    axs[i].tick_params(axis='x', labelrotation=60)
    axs[i].grid(True, color = 'k', alpha = 0.2)

plt.tight_layout()
plt.show()
