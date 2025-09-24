# ============================= IMPORTS
===================================== #
import os
#os.environ["OMP_NUM_THREADS"] = "1"
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import time
inicio=time.time()
# Own Libs
sys.path.append("/scratch/sunriseIII/TuMags_Reduction_Pipeline")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate, demodulate_quadrants
import alignment as al

# ============================= CONFIG ===================================== #

dc_paths = ih.get_images_paths("D10-5620-5719")
### ff_paths = ih.get_images_paths("D11-46203-47161") ### Otro flat
### Meter aqui el flat del modo correspondiente
ff_paths = ih.get_images_paths("D10-4340-5619")
Obs_paths = ih.get_images_paths("D10-304-2739") ### Todo el scan
#Obs_paths = ih.get_images_paths("D11-30968-32605") ### Todo el scan
#Obs_paths = ih.get_images_paths("D11-30968-32000") ### Todo el scan
###Obs_paths = ih.get_images_paths("D11-30968-31605")

# ======================= Darks and flats
================================== #

#dc = compute_master_darks(dc_paths[:-4], verbose = True)
#np.save('dc_01_QSUN_mode_1.npy',dc)
#dc = np.load('dc_01_QSUN_mode_1.npy')

#ff_data, ff_info = compute_master_flat_field(ff_paths,dc,  verbose =
True)
#np.save('ff_01_QSUN_mode_202.npy',ff_data)
### Los he guardado par no calcularlos cada vez.
#ff_data=np.load('ff_quietsun_mode202.npy')
#ff_data=np.load('ff_quietsun_mode1.npy')
# ======================= Selecting Observation Mode
======================= #

Ocs = ih.separate_ocs(Obs_paths, verbose = True, flat_fieldmode= False)
#sys.exit()
### occs para hacer la reduccion de todas las immagenes a la vez ###
### Modo 2.02 son los numeros impares en occs. Mode 1 los numeros pares
###
occs=np.arange(81,114,1)
occs1=np.zeros(int(len(occs)/2.))
occs202=np.zeros(int(len(occs)/2.))
for i in range(int(len(occs)/2.)):
     occs1[i]=occs[i*2]
     occs202[i]=occs[i*2+1]

### PARALELO ###
import concurrent.futures
def reduce_image(i):
     dc = np.load('dc_01_QSUN_mode_1.npy')
     ff_data = np.load('ff_01_QSUN_mode_202.npy')
     om = ih.nominal_observation("2.02", Ocs[int(i)]["ims"], dc)
     om_data = om.get_data()
     om_info = om.get_info()

     ff_cropped = ff_data[:, :, :, 300:-300, 300:-300]
     om_data_cropped = om_data[:, :, :, 300:-300, 300:-300]

     del ff_data, om_data

     om_corr = np.zeros(np.shape(om_data_cropped))
     for mod in range(om_info["Nmods"]):
         for lamb in range(om_info["Nlambda"]):
             om_corr[0, lamb, mod] = om_data_cropped[0, lamb, mod] /
ff_cropped[0, lamb, mod]
             om_corr[1, lamb, mod] = om_data_cropped[1, lamb, mod] /
ff_cropped[1, lamb, mod]

     quadrants = al.reshape_into_16_quadrants(om_corr,
om_info["Nlambda"], om_info["Nmods"])
     aligned, shifts = al.align_quadrants(quadrants)
     dual, demodulated = demodulate_quadrants(aligned,
om_info["Nlambda"], om_info["Nmods"], "517")

     # Guardar resultados (ajusta la ruta)

np.save(f'/scratch/sunriseIII/level1/20240710_speed_test/01_QSUN_Mode202_{int(i)}',
dual)

     return int(i)

# Ejecutar en paralelo
if __name__ == "__main__":
     with concurrent.futures.ProcessPoolExecutor(max_workers=12) as
executor:
         futures = executor.map(reduce_image, occs202)
         for result in futures:
             print(f"Finalizado: {result}")
fin=time.time()
print("Tiempo empleado (s):")
print(fin-inicio)