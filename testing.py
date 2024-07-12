import matplotlib.pyplot as plt
import numpy as np 

from field_stop_finder import find_fieldstop, fieldstopping_and_shifting

from utils import read_Tumag

cam1 = "/home/pablo/Desktop/TuMAGDATA/FS_modes_2_02/2024_07_11_22_17_00_484_0_542375.img"
cam2 = "/home/pablo/Desktop/TuMAGDATA/FS_modes_2_02/2024_07_11_22_17_00_484_1_541848.img"

c1, _ = read_Tumag(cam1)
c2, _ = read_Tumag(cam2)

c2 = np.flip(c2, axis = -1)

fsc1, fsc2 = find_fieldstop("auto", cam1 = c1, cam2 = c2, verbose = True, plot_flag=False)

fsc2[0][0] -= 2


c1, c2 = fieldstopping_and_shifting(c1, c2, fsc1, fsc2)


fig, axs = plt.subplots(1, 3, figsize = (15, 5))

c1 = c1 / np.max(c1)
c2 = c2 / np.max(c2)

axs[0].imshow(c1, cmap = "gray")
axs[1].imshow(c2, cmap = "gray")
axs[2].imshow(c1 - c2, cmap = "gray")

plt.show()