
#os.environ["OMP_NUM_THREADS"] = "1"

# tewst: python3 ../../TuMags_Reduction_Pipeline/process_data_main.py -f ../../TuMags_Reduction_Pipeline/config_tumag.yaml 
# 
# conda (with miniconda) enviroment creation:
# conda create --name TuMag
# conda activate TuMag
#
# ============================= IMPORTS ===================================== #
# import sys, os, argparse, logging
from sys import exit
from sys import path as syspath
from os import path, makedirs, listdir
from argparse import ArgumentParser
from pathlib import Path

#location of the tumag software:
syspath.append("/Users/orozco/IdAdA Dropbox/David orozco suárez/Python/TuMAG_codes/TuMags_Reduction_Pipeline")

import logging
import numpy as np
from tqdm import tqdm

from functools import partial
# import concurrent.futures
from concurrent.futures import ProcessPoolExecutor #, as_completed
import multiprocessing
# import matplotlib.pyplot as plt

#loading of TuMag software needed programs
import config as cf
import image_handler as ih 
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
from image_filtering import filter_frecuencies
from fits_files_handling import generate_fits, update_header
from astropy.io import fits
from pandas import read_csv
from alignment import align_obsmode
from advanced_alignment import apply_transform, update_alignment_csv, advance_alignment
from demodulation import demodulate
# from image_alignment_suit import image_alignment_affine
from xtalk_jaeggli import fit_mueller_matrix, fit_mueller_matrix_2d,fit_mueller_matrix_2d_interference
import pd_functions_v22 as phased
from destretch import destretch

from process_data_utils import (ConfigLoader,parse_range, 
                              print_shifts_by_cam, plt_darks, 
                              plt_flats,plt_level, _timestamp_from_filename, 
                              balance, parse_header_time, minutes_from_dt, 
                              interpolate_filter, correct_image_fft) #, format_dict_two_rows
from process_data_timelines import obs_dict # this brings into memory timeline_names and obs_dict

# Global variables
try:
    workspacePath = path.dirname(path.abspath(__file__))
except:
    workspacePath = './'
errLogFilename = f'{workspacePath}/errors.log'
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

# ======================= processing programs =======================
# =======================  reduce_image_0_5   =======================

def reduce_image_0_5(ocs, OCs, cfg, dc_real, ff_data, obs_ID, ff_paths, dc_paths, process_line_index):

     process_name = multiprocessing.current_process().name
     logging.basicConfig(
          level=logging.INFO,
          format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
          force=True  # resets logging config for each subprocess
     )

     logging.info(f' processing ocs: {ocs} ........... ')
     om_value = OCs[ocs]['OM']
     if cfg['process_line'] == om_value:
          logging.info(f' processing ocs: {ocs} which corresponds to {om_value} with {len(OCs[ocs]["ims"])} total images')
          # print(OCs[ocs]['ims'])
          if len(OCs[ocs]['ims']) != obs_dict[obs_ID]['obs_size'][process_line_index]:
               logging.error(f"  >> Error. The ocs: {ocs} number of images {len(OCs[ocs]['ims'])} does not coincide with the timeline info: {obs_dict[obs_ID]['obs_size'][process_line_index]}")
               return
               # sys.exit()

          obs_data = ih.nominal_observation(cfg['process_line'], OCs[ocs]["ims"], dc_real,modify_linearity=([1539,1540],[1.0,1.0]))
          data = obs_data.get_data()
          om_info = obs_data.get_info()  # Get observation mode info

          date_str = om_info["Images_headers"]["wv_0"]["M0"]["Date"].strftime("%d%m%YT%H%M%S")
          dataid = obs_ID+"_TM_"+cf.om_config[cfg['process_line']]["name"]+'_'+str(cf.om_config[cfg['process_line']]["Nlambda"])+"_"
          filename = dataid + date_str
          extended_filename = f"{filename}_LV_0.5_v{cfg['proc_version']}.fits"
          logging.info(f' Output filename: {extended_filename}')

          logging.info(f'  >> FF correction..........')
          data = np.where(np.isfinite(ff_data), data/ff_data, 0) 

          logging.info(f'  >> data cropping..........')
          data = data[:, :, :, cfg['centro'][0] - cfg['corte']:cfg['centro'][0] + cfg['corte'], cfg['centro'][1] - cfg['corte']:cfg['centro'][1] + cfg['corte']]

          if cfg['level_05']['filtering']:
               data = filter_frecuencies(data,band='fixed',verbose=True,pad=500,N=350,cam=1)

          # logging.info(f' Saving filename: {cfg['output_folder']+obs_ID}'/'{extended_filename}')

          generate_fits(data.astype(np.float32), 
               cfg['output_folder']+obs_ID+'/', extended_filename, '0.5', cfg['proc_version'], om_info = om_info, 
               zkes = None,  shifts = None, fitted_muller = None, datatype='SCIENCE',
               DARK_ID = dc_paths, FLAT_ID = ff_paths[process_line_index])

          if cfg['plots']['plot_level0_5']:
               plt_level(data, 
                         cfg['plots']['roi_plots'], 
                         cfg['output_folder']+obs_ID, 
                         f"{filename}_LV_0.5_v{cfg['proc_version']}",
                         '0.5','flat_corrected')

     
# =======================  reduce_image_0_7   =======================

def reduce_image_0_7(input_data_filename, cfg, df, line, from_label = "LV_0.5", to_label ="LV_0.7"):#ocs, OCs, cfg, dc_real, ff_data, obs_ID, ff_paths, dc_paths, process_line_index):

     try:
          align_sequence = cfg['level_07']['align_sequence']
     except:
          align_sequence = 0

     process_name = multiprocessing.current_process().name
     logging.basicConfig(
          level=logging.INFO,
          format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
          force=True  # resets logging config for each subprocess
     )

     logging.info(f' processing file: {input_data_filename} ')
     filename = path.splitext(path.basename(input_data_filename.replace(from_label, to_label)))[0]
     obs_ID = cfg['obs_ID']

     #read data
     with fits.open(input_data_filename) as hdul:
          data = hdul[0].data
          header = hdul[0].header

     cn, wn, pn, xs, ys = data.shape
     roi = cfg['level_07']['align_roi']
     if roi[1] == -1 or roi[-1] == -1:
          roi[1] = ys
          roi[3] = xs
          
     scale, gamma = balance(data[0, 0], data[1, 0], roi=roi)#, clip_percentiles=(40,60))
     data[1] = data[1] * scale
     logging.info(f"Balance (4x): scale={scale:.6f}, gamma={gamma:.6f}")

     timestamp = _timestamp_from_filename(input_data_filename)
     logging.info(f'  timestamp: {timestamp} ')

     try:
          if cfg['level_07']['advanced_alignment']:
               logging.info(f'  ADVANCE ROT MODE')
               advanced = cfg['level_07']['advanced_alignment']
          else:
               advanced = False
          if cfg['level_07']['advanced_wave']:
               advanced_wave = cfg['level_07']['advanced_wave']
          else:
               advanced_wave = 0
          if cfg['level_07']['advance_save']:
               advance_save = cfg['level_07']['advance_save']
          else:
               advance_save = False
          if cfg['level_07']['advance_overwrite']:
               advance_overwrite = cfg['level_07']['advance_overwrite']
          else:
               advance_overwrite = False
     except:
          advanced = False
          advanced_wave = 0
          advance_save = False
          advance_overwrite = False

     if cfg['debug']['filter_test']:
          cm,wv,pl,_,_ = data.shape
          data_ = np.copy(data)
          if line['line'] == '525.02':
               logging.info(f'  FILTERING NEW FREQ IN 525.02......')
               for cam in range(cm):
                    for wl in range(wv):
                         for pn in range(pl):
                              if wl == 0:
                                   kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24)
                                   ky = (-3, -4, -5, -2, -3, -4, -5, -2, -3, -4)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
                              if wl >= 0 and wl <= 6:
                                   kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24)
                                   ky = (3 , 4 , 5 , 2 , 3 , 4 , 5 , 2 , 3 , 4)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
                              if wl == 4:
                                   kx = (16, 17, 18, 16, 17, 18)
                                   ky = (-3, -3, -3, -4, -4, -4)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
                              if wl == 7:
                                   kx = (20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24)
                                   ky = (-4, -5, -6, -7, -4, -5, -6, -7, -3, -4, -5, -6, -7, -3, -4, -5, -6, -7, -3, -4, -5, -6, -7)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
               data = np.copy(data_)
               del data_
          if line['line'] == '525.06':
               pass
          if line['line'] == '517':
               logging.info(f'  FILTERING NEW FREQ IN 517......')
               for cam in range(cm):
                    for wl in range(wv):
                         for pn in range(pl):
                              if wl == 0:
                                   kx = (22, 22, 22, 23, 23, 23, 24, 24, 24)
                                   ky = (-2, -3, -4, -2, -3, -4, -2, -3, -4)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
                              if wl >= 7 and wl <= 8:
                                   kx = (22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24)
                                   ky = (-2, -3, -4, -2, -3, -4, -5, -2, -3, -4, -5)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
                              if wl == 9:
                                   kx = (20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24)
                                   ky = (-2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5, -2, -3, -4, -5)
                                   data_[cam,wl,pn,:,:] = \
                                        correct_image_fft(data[cam,wl,pn,:,:],kx=kx,ky=ky, 
                                                          h=1,mode='intensity')
               data = np.copy(data_)
               del data_               # kx = (-24, -24)
               # ky = (4, 4)
               # data = correct_data_fft(data,kx=kx,ky=ky, h=4,mode='intensity')
               
     if advanced:

          data, result = advance_alignment(data, line['line'], advanced_wave)

          # time stamp es el que corresponda con la longitud de onda. Así mejoro la base de datos.

          dt0 = parse_header_time(header[f'WV_{int(advanced_wave)}_M0'])
          dt1 = parse_header_time(header[f'WV_{int(advanced_wave)}_M1'])
          dt2 = parse_header_time(header[f'WV_{int(advanced_wave)}_M2'])
          dt3 = parse_header_time(header[f'WV_{int(advanced_wave)}_M3'])

          timestamp = int((minutes_from_dt(dt0) + minutes_from_dt(dt1) + minutes_from_dt(dt2) + minutes_from_dt(dt3)) / 4.0 )

          if advance_save:

               # Guardar resultado en CSV con UPDATE
               MODULE_DIR = Path(__file__).resolve().parent
               alignment_csv_filename = MODULE_DIR / "CALDATA" / "alignment_results.csv"
               update_alignment_csv(
                    alignment_csv_filename,
                    timestamp,
                    result.x,
                    result.fun,
                    overwrite = advance_overwrite,
               )

     else:

          result = interpolate_filter(df, timestamp)

          def _print_result_pretty(result, title="Resultado interpolado"):
               print("\n" + "="*70)
               print(f" {title}")
               print("="*70)
               for key, value in result.items():
                    print(f"{key:<12} : {value}")
               print("="*70 + "\n")

          _print_result_pretty(result)

          # result['rotation_angle_deg'] = 0.049591
          # result['translation_y'] =  -1.251800
          # result['translation_x'] =  -3.444528
          # result["center_y"] = 792.282173
          # result["center_x"] = 1108.208125
          # result["scale_x"] = 0.999956
          # result["scale_y"] = 0.999962
          # result["shear_x"] = -0.000269
          # result["shear_y"] = 0.000398

          for i in range(wn):
               for j in range(pn):
                    data[1, i, j] = apply_transform(
                              data[1, i, j], 
                              result['angle_deg'], #0.0675  #0.05192....
                              np.array([result['t_y'], result['t_x']]),
                              np.array([result["center_y"],result["center_x"]]),
                              scale_x=result["scale_x"], scale_y=result["scale_y"], 
                              shear_x=result["shear_x"], shear_y=result["shear_y"]
                              )

     if cfg['level_07']['align_mode'] == 'fourier':

          # dd = data.copy()
          try:
               verb = cfg['level_07']['align_verbose']
          except:
               verb = False

          try:
               print(cfg['debug'])
               debug = cfg['debug']
          except:
               debug = False

          try:
               print(cfg['level_07']['align_modulations'])
               align_modulations = cfg['level_07']['align_modulations']
          except:
               align_modulations = True

          data,shifts,_  = align_obsmode(data, acc=cfg['level_07']['align_accuracy'],
                                               verbose=verb,
                                               filter=line['line'],
                                               returnshifts=True,
                                               roi=roi,
                                               quadrants = cfg['level_07']['align_quadrants'],
                                               align_sequence = align_sequence,
                                               debug = debug,
                                               align_modulations = align_modulations)

          if cfg['level_07']['align_quadrants'] == 0:
               update_header(header, 'REALIGN', 1)
               update_header(header, 'ALIGMETH', 'sicairos')
               update_header(header, 'ALIGN_ID', f"roi = {roi[0]}:{roi[1]},{roi[2]}{roi[3]}")
               for i in range(wn):
                    for j in range(pn):
                         key = f'row0_{i}{j}'
                         val = float(shifts[i, 0, 0, j])
                         comment = f'cam 0 row shift wn/pn [{i},{j}]'
                         header[key] = (val, comment)
                         key = f'col0_{i}{j}'
                         val = float(shifts[i, 0, 1, j])
                         comment = f'cam 0 col shift wn/pn [{i},{j}]'
                         header[key] = (val, comment)
                         key = f'row1_{i}{j}'
                         val = float(shifts[i, 1, 0, j])
                         comment = f'cam 1 row shift wn/pn [{i},{j}]'
                         header[key] = (val, comment)
                         key = f'col1_{i}{j}'
                         val = float(shifts[i, 1, 1, j])
                         comment = f'cam 1 col shift wn/pn [{i},{j}]'
                         header[key] = (val, comment)

               print_shifts_by_cam(shifts)

     if cfg['level_07']['align_mode'] == 'destretch':
          # import 
          ## destrecching to be implementes
          # We will get the shift and rot from the destrieching components and use them in 
          data, _ = destretch(data,
                              n_iterations = cfg['level_07']['destretch_n_iterations'],
                              aling_cam=cfg['level_07']['destretch_aling_cam'],
                              ngrid=cfg['level_07']['destretch_ngrid'],
                              lr=cfg['level_07']['destretch_lr'],
                              lambda_tt=cfg['level_07']['destretch_lambda_tt'],
                              filter = line['line'])#aling_cam='partial')

     logging.info(f'  demodulation: ')


     # data[0], mmatrix = fit_mueller_matrix(data[0],  
     #                          method = cfg['level_10']['crosst_mode'],
     #                          norm = False,
     #                          verbose = cfg['level_10']['crosst_verbose'],
     #                          plots=cfg['level_10']['plot_crosst_method'], 
     #                          last_wvl=cfg['level_10']['crosst_last_wave'],
     #                          ctmethod='linfit',
     #                          pthresh=cfg['level_10']['crosst_threshold'],
     #                          region = cfg['level_10']['crosst_region'])
     # data[1], mmatrix = fit_mueller_matrix(data[1],  
     #                          method = cfg['level_10']['crosst_mode'],
     #                          norm = False,
     #                          verbose = cfg['level_10']['crosst_verbose'],
     #                          plots=cfg['level_10']['plot_crosst_method'], 
     #                          last_wvl=cfg['level_10']['crosst_last_wave'],
     #                          ctmethod='linfit',
     #                          pthresh=cfg['level_10']['crosst_threshold'],
     #                          region = cfg['level_10']['crosst_region'])

     # _, data_both = demodulate(data, line['line'],dmod_matrices = cfg['level_07']['demod_matrix'], BothCams=True)

     # try:
     # #     from astropy.io import fits

     #      hdu = fits.PrimaryHDU(data)
     #      hdu.header['COMMENT'] = 'I_cam1 - I_cam2 at last wavelength'
     #      out_name = "anted_dual_bean.fits"
     #      hdu.writeto(out_name, overwrite=True)
     #      logging.info(f"Wrote difference FITS: {out_name}")
     # except Exception as exc:
     #      logging.warning(f"Could not write FITS (cam1-cam2): {exc}")
     # plt_level(data_both, 
     #                cfg['plots']['roi_plots'], 
     #                cfg['output_folder']+obs_ID, 
     #                filename,
     #                '0.5','after_demod')

     # exit()

     data = demodulate(data, line['line'],dmod_matrices = cfg['level_07']['demod_matrix'],mode=cfg['level_07']['demod_mode'])

     # data = dual_align(data_both)

     # print(cfg['plots']['roi_plots'],cfg['output_folder']+obs_ID,filename)
     if cfg['plots']['plot_level0_7']:
          plt_level(data, 
                    cfg['plots']['roi_plots'], 
                    cfg['output_folder']+obs_ID, 
                    filename,
                    '0.7',
                    label = '_'+cfg['level_07']['align_mode']+'demodulated'+cfg['level_07']['add_level_07_label'])

     out_file = input_data_filename.replace(from_label, to_label+cfg['level_07']['add_level_07_label'])
     logging.info(f' Saving filename: {out_file}')
     with fits.open(input_data_filename) as hdu_list:
          hdu_list[0].data = data
          hdu_list[0].header = header
          hdu_list.writeto(out_file, overwrite=True)

# =======================  reduce_image_1_0   =======================

def reduce_image_1_0(input_data_filename, cfg, from_label = "LV_0.7", to_label = "LV_1.0"):

     process_name = multiprocessing.current_process().name
     logging.basicConfig(
          level=logging.INFO,
          format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
          force=True  # resets logging config for each subprocess
     )

     logging.info(f'  >> processing file: {input_data_filename} ')
     filename = path.splitext(path.basename(input_data_filename.replace(from_label, to_label)))[0]
     obs_ID = cfg['obs_ID']


     with fits.open(input_data_filename) as hdul:
          data = hdul[0].data
          header = hdul[0].header

     roi = cfg['level_10']['crosst_roi']

     # if cfg['plots']['plot_level1_0']:
     #      plt_level(data, 
     #                     cfg['plots']['roi_plots'], 
     #                     cfg['output_folder']+obs_ID, 
     #                     filename,
     #                     '1.0','demod_before_norma')

     #read data
     # 
     if cfg['level_10']['normalization'] == 0:
          norm_factor=np.median(data[-1,0,roi[0]:roi[1],roi[2]:roi[3]])
          logging.info(f'  >> Normalization factor: {norm_factor} ')
     else:
          norm_factor = cfg['level_10']['normalization']
     data = data / norm_factor

     # if cfg['plots']['plot_level1_0']:
     #      plt_level(data, 
     #                     cfg['plots']['roi_plots'], 
     #                     cfg['output_folder']+obs_ID, 
     #                     filename,
     #                     '1.0','demod_before')

     if cfg['level_10']['crosst_quadrants'] == 0:
          data, mmatrix = fit_mueller_matrix(data,  
                                   method = cfg['level_10']['crosst_mode'],
                                   norm = False,
                                   verbose = cfg['level_10']['crosst_verbose'],
                                   plots=cfg['level_10']['plot_crosst_method'], 
                                   last_wvl=cfg['level_10']['crosst_last_wave'],
                                   ctmethod='linfit',
                                   pthresh=cfg['level_10']['crosst_threshold'],
                                   region = cfg['level_10']['crosst_region'])
     else:
          data, _, _ = fit_mueller_matrix_2d(data,  
                                   method = cfg['level_10']['crosst_mode'],
                                   verbose = cfg['level_10']['crosst_verbose'],
                                   plots=cfg['level_10']['plot_crosst_method'], 
                                   last_wvl=cfg['level_10']['crosst_last_wave'],
                                   divisions= cfg['level_10']['crosst_quadrants'],
                                   pthresh = cfg['level_10']['crosst_threshold'],
                                   region = cfg['level_10']['crosst_region'])

     # Allow crosst_interference to be 0 (disabled), an integer (number of divisions)
     # or a numpy file path (.npy / .npz) containing precomputed slope/intercept (and optionally corrected data)
     ci = cfg['level_10']['crosst_interference']
     slope = None
     intercept = None

     # If it's a string, try to interpret as a filepath or as an int string
     if isinstance(ci, str):
          # try file first
          if path.isfile(ci):
               arr = np.load(ci, allow_pickle=True)
               # .npz file (NpzFile) supports keys
               divisions = arr['divisions']
               slope = arr['slope']
               intercept = arr['intercept']

               data = fit_mueller_matrix_2d_interference(
                    data,
                    verbose = cfg['level_10']['crosst_verbose'],
                    divisions= divisions,
                    slope = slope,
                    intercept = intercept)
               
          else:
               # try to parse string as integer
               try:
                    ci = int(ci)
               except Exception as e:
                    raise ValueError(f"crosst_interference is a string but not a file nor an int: {ci}") from e

     # If ci is numeric and non-zero call the computation routine
     elif isinstance(ci, (int, np.integer)) and ci != 0:
          data, slope, intercept = fit_mueller_matrix_2d_interference(
                                   data,
                                   verbose = cfg['level_10']['crosst_verbose'],
                                   divisions= int(ci),
                                   pthresh = cfg['level_10']['crosst_threshold'],
                                   pthresh_intensity = cfg['level_10']['crosst_intensity_threshold'])
          
          update_header(header, 'CROSTALK', 1, after = 'ALIGMETH', comment='Was crosstalk correction applied?')
          header['INTERC'] = (cfg['level_10']['crosst_interference'], 'cuadrants for crosstalk interference')
          for i in range(slope.shape[0]):
               for j in range(slope.shape[1]):
                    key = f'IS_{i}{j}'
                    val = float(slope[i, j])
                    comment = f'Slope interference element wave, q [{i},{j}]'
                    header[key] = (val, comment)
          for i in range(intercept.shape[0]):
               for j in range(intercept.shape[1]):
                    key = f'II_{i}{j}'
                    val = float(intercept[i, j])
                    comment = f'intercept interference element wave, q [{i},{j}]'
                    header[key] = (val, comment)

     else:
          pass

     if cfg['plots']['plot_level1_0']:
          plt_level(data, 
                         cfg['plots']['roi_plots'], 
                         cfg['output_folder']+obs_ID, 
                         filename,
                         '1.0','demod'+cfg['level_10']['add_level_10_label'])

     if cfg['level_10']['crosst_mode'] == 'jaeggli' and cfg['level_10']['crosst_quadrants'] == 0:
          update_header(header, 'CROSTALK', 1, after = 'ALIGMETH', comment='Was crosstalk correction applied?')
          for i in range(4):
               for j in range(4):
                    key = f'MMAT_{i}{j}'
                    val = float(mmatrix[i, j])
                    comment = f'Mueller matrix element [{i},{j}]'
                    header[key] = (val, comment)

     if cfg['level_10']['crosst_mode'] == 'standard' and cfg['level_10']['crosst_quadrants'] and not cfg['level_10']['crosst_last_wave'] == 0:

          update_header(header, 'CROSTALK', 1, after = 'ALIGMETH', comment='Was crosstalk correction applied?')
          key = f'SQ'
          val = float(mmatrix[3])
          comment = f'slope Q'
          header[key] = (val, comment)
          key = f'SU'
          val = float(mmatrix[4])
          comment = f'slope U'
          header[key] = (val, comment)
          key = f'SV'
          val = float(mmatrix[5])
          comment = f'slope V'
          header[key] = (val, comment)
          key = f'IQ'
          val = float(mmatrix[0])
          comment = f'offset Q'
          header[key] = (val, comment)
          key = f'IU'
          val = float(mmatrix[1])
          comment = f'offset U'
          header[key] = (val, comment)
          key = f'IV'
          val = float(mmatrix[2])
          comment = f'offset V'
          header[key] = (val, comment)

     if cfg['level_10']['crosst_mode'] == 'standard' and cfg['level_10']['crosst_quadrants'] != 0:
          pass 

     out_file = input_data_filename.replace(from_label, to_label+cfg['level_10']['add_level_10_label'])
     logging.info(f' Saving filename: {out_file}')
     with fits.open(input_data_filename) as hdu_list:
          hdu_list[0].data = data
          hdu_list[0].header = header
          hdu_list.writeto(out_file, overwrite=True)


# =======================  reduce_image_1_1   =======================

def reduce_image_1_1(input_data_filename, cfg, zk = None):

     process_name = multiprocessing.current_process().name
     logging.basicConfig(
          level=logging.INFO,
          format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
          force=True  # resets logging config for each subprocess
     )

     logging.info(f'  >> processing file: {input_data_filename} ')
     filename = path.splitext(path.basename(input_data_filename.replace("LV_1.0", "LV_1.1")))[0]
     obs_ID = cfg['obs_ID']

     #read data
     with fits.open(input_data_filename) as hdul:
          data = hdul[0].data
          header = hdul[0].header

     wn, pn, xs, ys = data.shape

     with tqdm(total=wn*pn) as pbar:
          for wl in range(wn):
               for pl in range(pn):
                    # test,_ = pd.restore_ima(data[wl,pl],
                    #     zk,pd=0,low_f=0.2,noise='default',reg1=0.05,reg2=1,cobs=32.4)
                    if pl == 0:
                         data[wl,pl],noise_filter = phased.restore_ima(data[wl,pl],
                              zk,pd=0,low_f=0.2,reg1=0.05,reg2=1,cobs=32.4, epsilon=0.02, sigma= 5000, stray='moffat')
                    else:
                         data[wl,pl],_ = phased.restore_ima(data[wl,pl],
                              zk,pd=0,low_f=0.2,noise=noise_filter,reg1=0.05,reg2=1,cobs=32.4, epsilon=0.02, sigma= 5000, stray='moffat')

               pbar.update(1)


     if cfg['plots']['plot_level1_1']:
          plt_level(data, 
                         cfg['plots']['roi_plots'], 
                         cfg['output_folder']+obs_ID, 
                         filename,
                         '1.0','pd')

     with fits.open(input_data_filename) as hdu_list:
          hdu_list[0].data = data
          hdu_list[0].header = header
          hdu_list.writeto(input_data_filename.replace("LV_1.0", "LV_1.1"), overwrite=True)

# =======================  reduce_image_0_7pd   =======================

def reduce_image_0_6(input_data_filename, cfg, zk = None, from_label = "LV_0.5", to_label ="LV_0.6"):

     process_name = multiprocessing.current_process().name
     logging.basicConfig(
          level=logging.INFO,
          format=f'%(asctime)s [{process_name}] %(levelname)s: %(message)s',
          force=True  # resets logging config for each subprocess
     )

     logging.info(f'  >> processing file: {input_data_filename} ')
     filename = path.splitext(path.basename(input_data_filename.replace(from_label, to_label)))[0]
     obs_ID = cfg['obs_ID']

     #read data
     with fits.open(input_data_filename) as hdul:
          data = hdul[0].data
          header = hdul[0].header

     cam, wn, pn, xs, ys = data.shape

     with tqdm(total=wn*pn*cam) as pbar:
          for cm in range(cam):
               for wl in range(wn):
                    for pl in range(pn):
                         # test,_ = pd.restore_ima(data[wl,pl],
                         #     zk,pd=0,low_f=0.2,noise='default',reg1=0.05,reg2=1,cobs=32.4)
                         if pl == 0:
                              data[cm,wl,pl],noise_filter = phased.restore_ima(data[cm,wl,pl],
                                   zk,pd=0,low_f=0.2,reg1=0.05,reg2=1,cobs=32.4, epsilon=0.02, sigma= 5000, stray='moffat')
                         else:
                              data[cm,wl,pl],_ = phased.restore_ima(data[cm,wl,pl],
                                   zk,pd=0,low_f=0.2,noise=noise_filter,reg1=0.05,reg2=1,cobs=32.4, epsilon=0.02, sigma= 5000, stray='moffat')

                         pbar.update(1)

     if cfg['plots']['plot_level0_5']:
          plt_level(data, 
                    cfg['plots']['roi_plots'], 
                    cfg['output_folder']+obs_ID, 
                    f"{filename}_{from_label}_v{cfg['proc_version']}",
                    '0.5','pd')

     with fits.open(input_data_filename) as hdu_list:
          hdu_list[0].data = data
          hdu_list[0].header = header
          hdu_list.writeto(input_data_filename.replace(from_label, to_label), overwrite=True)

# ===================================================================
# =======================     MAIN PROGRAM    =======================
# ===================================================================

if __name__ == "__main__":

     logging.info('-----------------------------------')
     logging.info('  >> Running process_data_main ')

     parser = ArgumentParser(description="Process TuMags data.")
     parser.add_argument("-f", "--config", type=str, default="config_tumag.yaml", help="Path to config file")
     args = parser.parse_args()

     # load information from config_tumag.yaml or provided file
     cfg = ConfigLoader(args.config)
     cfg_dict = cfg.config  # get its underlying dict

     obs_ID = cfg.obs_ID
     logging.info(f"  >> Input obs_ID: {obs_ID}")
     logging.info(f"  >> Folder data location: {cfg.tumag_data_location}")
     logging.info(f"  >> Output folder location: {cfg.output_folder}")
     

     if cfg.force_redo['level0_5'] or cfg.force_redo['redo_flat'] or cfg.force_redo['redo_dark']:

          # first thing is to get the paths.
          # darks
          dc_paths = obs_dict[obs_ID]['darks'][0] # usually one ID for dark but may be a problem
          # darks
          ff_paths = obs_dict[obs_ID]['flats']  # usually more than one ID for flats but may be a problem
          # obs
          obs_paths = obs_dict[obs_ID]['obsdata'][0] # usually one ID for obs but may be a problem

          logging.info(f"  >> Darks: {dc_paths}, flats: {ff_paths}, obs: {obs_paths} ")
          logging.info('-----------------------------------')
          logging.info(f'  >> checking if {cfg.Organized_files_local_folder_name} exist and continuing')
          if not path.exists(cfg.tumag_data_location + cfg.Organized_files_local_folder_name):
               logging.error(f' >> the Organized_files_local_folder_name: {cfg.Organized_files_local_folder_name} not found')
               exit()

          ih.Organization_folder_files = path.join(cfg.tumag_data_location, cfg.Organized_files_local_folder_name) 

          logging.info(f'  >> checking if out folder {cfg.output_folder+obs_ID} exist and creating it if it does not')
          if not path.exists(cfg.output_folder+obs_ID):
               makedirs(cfg.output_folder+obs_ID)
               makedirs(cfg.output_folder+obs_ID+'/pngs')
          logging.info('-----------------------------------')

          # DARKS PROCESSING

          logging.info('  >> Running dark calculation ')
          dark_output_file = cfg.output_folder+obs_ID+'/'+dc_paths+'.npz'
          logging.info(f'  >> dark output file will be: {dark_output_file}')

          if path.exists(dark_output_file):
               logging.info(f'  >> dark output file exist.')
               if not cfg.force_redo["redo_dark"]:
                    logging.info(f'  >> loading dark')
                    dark = np.load(dark_output_file, allow_pickle=True)
                    dc_real = dark['dc_real']
               else:
                    logging.info(f'  >> but force_redo is True.')
                    dark_paths = ih.get_images_paths(dc_paths)
                    dc_real, head, rms_dark = compute_master_darks(dark_paths[cfg.darks_indexes["dark_from"]:cfg.darks_indexes["dark_to"]], verbose = True)
                    np.savez(dark_output_file,dc_real=dc_real.astype(np.float32))
          else:
               logging.info(f' >> dark output file does not exist.')
               dark_paths = ih.get_images_paths(dc_paths)
               dc_real, head, rms_dark = compute_master_darks(dark_paths[cfg.darks_indexes["dark_from"]:cfg.darks_indexes["dark_to"]], verbose = True)
               np.savez(dark_output_file,dc_real=dc_real.astype(np.float32))

          if cfg.plots["plot_darks"]:
               plt_darks(dc_real)

          # FLATS PROCESSING
          logging.info('-----------------------------------')

          logging.info('  >> Running FLAT calculation ')


          # Check if cfg.process_line exists in obs_dict[obs_ID]['obs_is'] and get its position
          try:
               obs_is_list = obs_dict[obs_ID]['obs_is']
               if cfg.process_line in obs_is_list:
                    process_line_index = obs_is_list.index(cfg.process_line)
                    logging.info(f"  >> process_line '{cfg.process_line}' found at position {process_line_index} in obs_is list.")
               else:
                    logging.error(f"  >> process_line '{cfg.process_line}' not found in obs_is list: {obs_is_list}")
                    exit()
          except Exception as e:
               logging.error(f"  >> Error checking process_line in obs_is list: {e}")
               exit()
          # Check if cfg.process_line exists in obs_dict[obs_ID]['obs_is'] and get its position
          if cfg.flat_file:
               logging.info(f'  >> reading flat: {cfg.flat_file}')
               flat = np.load(cfg.flat_file, allow_pickle=True)
               ff_data = flat['ff_data']
               ff_info = flat['ff_info']

          else:

               flat_output_file = cfg.output_folder+obs_ID+'/'+ff_paths[process_line_index]+'.npz'
               logging.info(f'  >> flat output file will be: {flat_output_file}')

               if path.exists(flat_output_file):
                    logging.info(f'  >> flat output file exist.')
                    if not cfg.force_redo["redo_flat"]:
                         logging.info(f'  >> loading flat')
                         flat = np.load(flat_output_file, allow_pickle=True)
                         ff_data = flat['ff_data']
                         ff_info = flat['ff_info']
                    else:
                         logging.info(f'  >> but force_redo is True.')
                         flat_paths = ih.get_images_paths(ff_paths[process_line_index])
                         ff_data, ff_info = compute_master_flat_field(flat_paths, dc = dc_real, verbose = True,
                                                            modify_linearity=([1539,1540],[1.0,1.0]),
                                                            norm_roi = cfg.flat_norm_roi,
                                                            norm_method = cfg.norm_method, 
                                                            remove_prefilter = cfg.remove_prefilter,
                                                            pref_model = cfg.pref_model, 
                                                            import_blueshift_guess=cfg.import_blueshift_guess)
                                                            # remove_prefilter=True,
                                                            # pref_model="prefilter_model_517.pkl")
                         np.savez(flat_output_file,ff_data=ff_data.astype(np.float32),ff_info=ff_info)
               else:
                    logging.info(f' >> flat output file does not exist.')
                    flat_paths = ih.get_images_paths(ff_paths[process_line_index])
                    ff_data, ff_info = compute_master_flat_field(flat_paths, dc = dc_real, verbose = True,
                                                            modify_linearity=([1539,1540],[1.0,1.0]),
                                                            norm_roi = cfg.flat_norm_roi,
                                                            norm_method = cfg.norm_method, 
                                                            remove_prefilter = cfg.remove_prefilter,
                                                            pref_model = cfg.pref_model, 
                                                            import_blueshift_guess=cfg.import_blueshift_guess)#,
                                                            # remove_prefilter=True,
                                                            # pref_model="prefilter_model_517.pkl")
                    np.savez(flat_output_file,ff_data=ff_data.astype(np.float32),ff_info=ff_info)

               if cfg.plots["plot_flats"]:
                    plt_flats(ff_data)

          
          logging.info('-----------------------------------')

          # DATA PROCESSING 1) First is needed to check the OCs.
          if cfg.force_redo['level0_5']:
               ocs_output_file = cfg.output_folder+obs_ID+'/'+obs_ID+'_ocs.npz'
               logging.info(f'  >> ocs output file will be: {ocs_output_file}')

               if path.exists(ocs_output_file):
                    logging.info(f'  >> ocs output file exist.')
                    OCs = np.load(ocs_output_file,allow_pickle=True)['OCs'][()]
               else:
                    logging.info(f'  >> determining ocs....')
                    obs_images = ih.get_images_paths(obs_paths)
                    OCs = ih.separate_ocs(obs_images, verbose = False)
                    np.savez(ocs_output_file,OCs=OCs)

               # DATA PROCESSING 2) process data. Can be one file or many in parallel

               process_ocs = list(OCs.keys())
               logging.info(f'  >> available ocs {len(process_ocs)}')
               if len(process_ocs) == 0:
                    logging.error(f'  >> Error. No ocs found for obs_ID: {obs_ID}')
                    exit()
               # process_ocs = process_ocs[parse_range(cfg.process_ocs)]
               logging.info(f'  >> process ocs {process_ocs}')

               process_ocs = [process_ocs[i] for i in parse_range(cfg.process_ocs,max_value = len(process_ocs)-1)]

               logging.info(f'  >> process ocs {process_ocs}')

               # def _entry(oc, OCs, cfg, dc_real, ff_data, obs_ID, ff_paths, dc_paths, process_line_index):
               # # delega al original:
               #      return reduce_image_0_5(
               #           oc, OCs=OCs, cfg=cfg, dc_real=dc_real, ff_data=ff_data,
               #           obs_ID=obs_ID, ff_paths=ff_paths, dc_paths=dc_paths,
               #           process_line_index=process_line_index
               #      )

               # Ensure process_ocs is always a list
               if not isinstance(process_ocs, list):
                    process_ocs = [process_ocs]

               reduce_partial = partial(
                    reduce_image_0_5,
                    OCs=OCs,
                    cfg=cfg_dict,
                    dc_real=dc_real,
                    ff_data=ff_data,
                    obs_ID=obs_ID,
                    ff_paths=ff_paths,
                    dc_paths=dc_paths,
                    process_line_index=process_line_index
                    )

               if len(process_ocs) > 1 and cfg.parallel:

                    # with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    #      futures = {
                    #           executor.submit(
                    #                _entry, oc, OCs, cfg_dict, dc_real, ff_data, obs_ID, ff_paths, dc_paths, process_line_index
                    #           ): oc
                    #           for oc in process_ocs
                    #      }
                    #      for fut in as_completed(futures):
                    #           oc = futures[fut]
                    #           try:
                    #                fut.result()
                    #           except Exception:
                    #                logging.exception(f"Falló OC={oc}")

                    with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                         executor.map(reduce_partial, process_ocs)

               elif len(process_ocs) > 1 and not cfg.parallel:
                    for i in process_ocs:
                         reduce_partial(i)
                    # for oc in process_ocs:
                    #      reduce_image_0_5(
                    #           oc, OCs=OCs, cfg=cfg_dict, dc_real=dc_real, ff_data=ff_data,
                    #           obs_ID=obs_ID, ff_paths=ff_paths, dc_paths=dc_paths,
                    #           process_line_index=process_line_index
                    #      )
               else:
                    reduce_partial(process_ocs[0])
                    # reduce_image_0_5(
                    #      process_ocs[0], OCs=OCs, cfg=cfg_dict, dc_real=dc_real, ff_data=ff_data,
                    #      obs_ID=obs_ID, ff_paths=ff_paths, dc_paths=dc_paths,
                    #      process_line_index=process_line_index
                    # )

     logging.info('-----------------------------------')
     try:
          deconvolve_first = cfg.force_redo['level0_6']
          from_label = "LV_0.5"
          to_label = "LV_0.6"
     except:
          deconvolve_first = False

     if deconvolve_first:

          logging.info(f'  >> Procesing level 0.6')
          cfg.force_redo['level1_1'] = False

          dataid = obs_ID+"_TM_"+cf.om_config[cfg.process_line]["name"]+'_'+str(cf.om_config[cfg.process_line]["Nlambda"])+"_"
          ext_ = f"_{from_label}_v{cfg.proc_version}.fits"

          directory = cfg.output_folder+obs_ID+'/'
          files_list = sorted(listdir(directory))
          files_list = [f for f in files_list if f.endswith('.fits') and not f.startswith('._')]
          files_list = [f for f in files_list if dataid in f ]
          files_list = [f for f in files_list if ext_ in f ]
          files_list = [path.join(directory, f) for f in files_list]
          indices = parse_range(cfg.process_files,max_value = len(files_list)-1)
          if len(files_list) >= 1 and len(indices) <= len(files_list):
               files = [files_list[i] for i in indices]
               logging.info(f'  >> available files {len(files)}')
          if len(files_list) == 0:
               logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
               exit()
          if len(indices) > len(files_list):
               logging.error(f'  >> More indices than files {indices} len of file list {len(files_list)}')
               exit()

          # Ensure process_ocs is always a list
          if not isinstance(files, list):
               files = [files]

          print('Importing zernikes: ', cfg.zernike_id)
          BASE_DIR = Path(__file__).resolve().parent
          file_path_csv = BASE_DIR / "TuMag_PD_results_All_filters_clean.csv"          
          zk = phased.import_zernikes(cfg.zernike_id,csv_path=file_path_csv) #06_SPOT_Fe2.02_0,06_SPOT_Mg1_0,06_SPOT_Fe2.02_1,06_SPOT_Mg1_1

          reduce_partial = partial(
               reduce_image_0_6,
               cfg=cfg_dict,
               zk = zk,
               from_label = from_label, to_label = to_label
               )

          if len(files) > 1 and cfg.parallel:
               with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    executor.map(reduce_partial, files)
          elif len(files) > 1 and not cfg.parallel:
               for i in files:
                    reduce_partial(i)
          else:
               reduce_partial(files[0])

     logging.info('-----------------------------------')

     try:
          if cfg.force_redo['use_pd']:
               cfg.force_redo['level1_1'] = False
               from_label = "LV_0.6"
               to_label = "LV_0.8"
          else:
               from_label = "LV_0.5"
               to_label = "LV_0.7"
     except:
          from_label = "LV_0.5"
          to_label = "LV_0.7"
   
     if cfg.force_redo['level0_7']:

          logging.info(f'  >> processing level 0.7 (alignment and demodulation)')

          dataid = obs_ID+"_TM_"+cf.om_config[cfg.process_line]["name"]+'_'+str(cf.om_config[cfg.process_line]["Nlambda"])+"_"
          ext_ = f"_{from_label}_v{cfg.proc_version}.fits"

          directory = cfg.output_folder+obs_ID+'/'
          files_list = sorted(listdir(directory))
          files_list = [f for f in files_list if f.endswith('.fits') and not f.startswith('._')]
          files_list = [f for f in files_list if dataid in f ]
          files_list = [f for f in files_list if ext_ in f ]
          files_list = [path.join(directory, f) for f in files_list]
          indices = parse_range(cfg.process_files,max_value = len(files_list)-1)
          if len(files_list) >= 1 and len(indices) <= len(files_list):
               files = [files_list[i] for i in indices]
               logging.info(f'  >> available files {len(files)}')
          if len(files_list) == 0:
               logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
               exit()
          if len(indices) > len(files_list):
               logging.error(f'  >> More indices than files {indices} len of file list {len(files_list)}')
               exit()

          BASE_DIR = Path(__file__).resolve().parent
          file_path_csv = BASE_DIR / "CALDATA" / f"{cfg.level_07['align_rot_data_filter']}"          

          df = read_csv(file_path_csv)  # <-- make sure the file path is correct
          # Orden temporal
          df = df.sort_values('timestamp').reset_index(drop=True)
          
          # Ensure process_ocs is always a list
          if not isinstance(files, list):
               files = [files]

          reduce_partial = partial(
               reduce_image_0_7,
               cfg=cfg_dict,
               df=df,
               line = cf.om_config[cfg.process_line],
               from_label = from_label, to_label = to_label
               )

          if len(files) > 1 and cfg.parallel:
               with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    executor.map(reduce_partial, files)
          elif len(files) > 1 and not cfg.parallel:
               for i in files:
                    reduce_partial(i)
          else:
               reduce_partial(files[0])

     logging.info('-----------------------------------')
     if cfg.force_redo['level1_0']:

          try:
               if cfg.force_redo['use_pd']:
                    cfg.force_redo['level1_1'] = False
                    from_label = "LV_0.8"
                    to_label = "LV_1.2"
               else:
                    from_label = "LV_0.7"
                    to_label = "LV_1.0"

          except:
               from_label = "LV_0.7"
               to_label = "LV_1.0"
     
          logging.info(f'  >> processing level 1.0')
          # cmatrix = cfg.mmatrix 

          dataid = obs_ID+"_TM_"+cf.om_config[cfg.process_line]["name"]+'_'+str(cf.om_config[cfg.process_line]["Nlambda"])+"_"
          ext_ = f"_{from_label}{cfg.level_07['add_level_07_label']}_v{cfg.proc_version}.fits"

          directory = cfg.output_folder+obs_ID+'/'
          files = sorted(listdir(directory))
          files = [f for f in files if f.endswith('.fits') and not f.startswith('._')]
          files = [f for f in files if dataid in f ]
          files = [f for f in files if ext_ in f ]
          files = [path.join(directory, f) for f in files]
          if len(files) > 1:
               files = [files[i] for i in parse_range(cfg.process_files,max_value = len(files)-1)]
          logging.info(f'  >> available files {len(files)}')
          if len(files) == 0:
               logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
               exit()

          # Ensure process_ocs is always a list
          if not isinstance(files, list):
               files = [files]

          reduce_partial = partial(
               reduce_image_1_0,
               cfg=cfg_dict,
               from_label = from_label,
               to_label = to_label
               )

          if len(files) > 1 and cfg.parallel:
               with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    executor.map(reduce_partial, files)
          elif len(files) > 1 and not cfg.parallel:
               for i in files:
                    reduce_partial(i)
          else:
               reduce_partial(files[0])

     logging.info('-----------------------------------')
     if cfg.force_redo['level1_1']:

          logging.info(f'  >> processing level 1.1')
          # cmatrix = cfg.mmatrix 

          dataid = obs_ID+"_TM_"+cf.om_config[cfg.process_line]["name"]+'_'+str(cf.om_config[cfg.process_line]["Nlambda"])+"_"
          ext_ = f"_LV_1.0{cfg.level_07['add_level_07_label']}_v{cfg.proc_version}.fits"

          directory = cfg.output_folder+obs_ID+'/'
          files = sorted(listdir(directory))
          files = [f for f in files if f.endswith('.fits') and not f.startswith('._')]
          files = [f for f in files if dataid in f ]
          files = [f for f in files if ext_ in f ]
          files = [path.join(directory, f) for f in files]
          print(files)
          if len(files) > 1:
               files = [files[i] for i in parse_range(cfg.process_files,max_value = len(files)-1)]
          logging.info(f'  >> available files {len(files)}')
          if len(files) == 0:
               logging.error(f'  >> Error. No files found for obs_ID: {obs_ID}')
               exit()

          # Ensure process_ocs is always a list
          if not isinstance(files, list):
               files = [files]
          ######## 
          print('Importing zernikes: ', cfg.zernike_id)
          BASE_DIR = Path(__file__).resolve().parent
          file_path_csv = BASE_DIR / "TuMag_PD_results_All_filters_clean.csv"          
          zk = phased.import_zernikes(cfg.zernike_id,csv_path=file_path_csv) #06_SPOT_Fe2.02_0,06_SPOT_Mg1_0,06_SPOT_Fe2.02_1,06_SPOT_Mg1_1

          reduce_partial = partial(
               reduce_image_1_1,
               cfg=cfg_dict,
               zk = zk
               )

          if len(files) > 1 and cfg.parallel:
               with ProcessPoolExecutor(max_workers=cfg.max_workers) as executor:
                    executor.map(reduce_partial, files)
          elif len(files) > 1 and not cfg.parallel:
               for i in files:
                    reduce_partial(i)
          else:
               reduce_partial(files[0])



