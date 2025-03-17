# ---------------------------- DESCRIPTION --------------------------------------- #

"""
This module provides functionality to handle and process .img files. 
It includes classes for reading image headers and metadata and managing observations.

author: Pablo Santamarina Guerrero (pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import time
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime

# Own libs
from utils import read_Tumag
import config as cf

# Config
current_dir = os.path.dirname(os.path.abspath(__file__))
Organization_folder_files = os.path.join(current_dir, "Organized_files")

# ------------------------------  CODE  ------------------------------------------ # 

# Class to create header object for images from the raw values. 
class raw_header:

    # Transfer function for hvps commanded from counts to volts 
    def hvps_commanded_2_volts(self, counts : int, sign : int) -> float:
        return -1 ** (sign + 1) * (counts * 4999 / (2 ** 12 - 1))
    
    # Transfer function for hvps read from counts to volts
    def hvps_read_2_volts(self, counts : int) -> float:
        return 5000 * (((2 * counts) / 4095) - 1)
    
    # Transfer function for lcvr 1 from counts to volts
    def lcvr1_counts_2_volts(self,counts : int) -> float:
        return counts * 0.00179439 
    
    # Transfer function for lcvr 2 from counts to volts
    def lcvr2_counts_2_volts(self, counts : int) -> float:
        return counts * 0.00179507 

    def __init__(self, camID : int, om : int, nAcc : int, Roix : int, Roiy : int, 
                 Roix_offset : int, Roiy_offset : int, oc : int, fw1 : int, fw2 : int,
                 hvps_counts : int, hvps_sign : int, lcvr1_counts : int, lcvr2_counts : int,
                 hvps_read_counts : int, lcvr1_read_counts : int, lcvr2_read_counts : int, image_name : str) -> None:
        # Initialize information dictionary
        self.info = {
            "cam" : camID,
            "ObservationMode_index" : om,
            "ObservationMode" : cf.observation_modes[om],
            "nAcc" : nAcc,
            "Roix" : Roix,
            "Roiy" : Roiy,
            "Roix_offset" : Roix_offset,
            "Roiy_offset" : Roiy_offset,
            "ObservationCounter" : oc,
            "FW1_ind" : fw1,
            "FW2_ind" : fw2,
            "FW1" : cf.filter_wheel["1"][fw1],
            "FW2" : cf.filter_wheel["2"][fw2],
            "hvps_counts" : hvps_counts,
            "hvps_sign" : hvps_sign,
            "lcvr1_counts" : lcvr1_counts,
            "lcvr2_counts" : lcvr2_counts,
            "hvps_read_counts" : hvps_read_counts,
            "lcvr1_read_counts" : lcvr1_read_counts,
            "lcvr2_read_counts" : lcvr2_read_counts,
            "t_exp" : 42,
            "image_name" : image_name,
            "Date" : get_time_from_filename(image_name)

        }
        # Compute hvps commanded counts to volts
        self.info["hvps_comm_volts"] = self.hvps_commanded_2_volts(hvps_counts, hvps_sign)
        # Compute hvps read counts to volts
        self.info["hvps_read_volts"] = self.hvps_read_2_volts(hvps_read_counts)
        # Compute lcvr1  counts to volts
        self.info["lcvr1_volts"] = self.lcvr1_counts_2_volts(lcvr1_counts)
        # Compute lcvr2  counts to volts
        self.info["lcvr2_volts"] = self.lcvr2_counts_2_volts(lcvr2_counts)

    # Return dictionary
    def get_info(self):
        return self.info

# Read single images
def read(image_path : str):
    # Read .img file
    img, h = read_Tumag(image_path) 
    # Create header object
    head = raw_header(int(h["CameraID"]), int(h["ObservationMode"]), int(h["nAcc"]), 
                      int(h["Roi_x_size"]), int(h["Roi_y_size"]), int(h["Roi_x_offset"]),
                      int(h["Roi_y_offset"]), int(h["Observation_Counter"]), int(h["FW1"]),
                      int(h["FW2"]), int(h["EtalonDN"]), int(h["EtalonSign"]), int(h["Rocli1_LCVR"]),
                      int(h["Rocli2_LCVR"]), int(h["EtalonVoltsReading"]), int(h["LCVR1_DN_Real"]),
                      int(h["LCVR2_DN_Real"]), os.path.basename(image_path))
    
    return img, head.get_info()

# Class to process the observation mode -> headers and array 
class nominal_observation:

    def __init__(self, om, images_path, dc):

        self.info = {"ObservationMode" : om,
                     "Images_headers" : {}}

        self.data = np.zeros((2, # Number of cameras 
                              cf.om_config[om]["Nlambda"],  # Number of wavelengths
                              cf.om_config[om]["Nmods"],    # Number of Modulations
                              cf.xsize, cf.ysize))  # Size of image (x, y)

        nmods   = cf.om_config[om]["Nmods"]     # N mods from config file
        nlambda = cf.om_config[om]["Nlambda"]   # N wavelengths from config file
        
        images_path_reshaped = np.array(images_path).reshape(nlambda, nmods, 2)

        _, h1 = read(images_path_reshaped[0, 0, 0]) # read first image to get acc
       
        for lambd in range(nlambda):
            print(f"Processing wavelength : {lambd + 1} / {nlambda}")
            self.info["Images_headers"][f"wv_{lambd}"] = {}
            for mod in range(nmods):
                # Reading each image
                im0, head0 = read(images_path_reshaped[lambd, mod, 0]) # Cam 1
                im1, _ = read(images_path_reshaped[lambd, mod, 1]) # Cam 2
                # Saving images header except for CameraID entry
                self.info["Images_headers"][f"wv_{lambd}"][f"M{mod}"] = {}
                for key in head0:
                    if key == "cam":
                        pass
                    else:
                        self.info["Images_headers"][f"wv_{lambd}"][f"M{mod}"][key] = head0[key]
            
                # Correct dark-current and save image data into main data array
                self.data[0, lambd, mod] = im0 - (dc[0] * head0["nAcc"])
                self.data[1, lambd, mod] = np.flip(im1, axis = -1) - (dc[1] * head0["nAcc"]) # Flip cam 2 image. 
        
        # Completing info of Observation Mode with info from header
        self.info["nAcc"] = head0["nAcc"]
        self.info["Roix"] = head0["Roix"]
        self.info["Roiy"] = head0["Roiy"]
        self.info["Roix_offset"] = head0["Roix_offset"]
        self.info["Roiy_offset"] = head0["Roiy_offset"]
            
        # Saving info from config file into observtaion mode info
        for entry in cf.om_config[head0["ObservationMode"]]:
            self.info[entry] = cf.om_config[head0["ObservationMode"]][entry]

    # Get information of the observation mode and individual headers
    def get_info(self):
        return self.info
    
    # Get the data array
    def get_data(self):
        return self.data
    
# Class to process flat-field modes -> headers and array 
class nominal_flat:

    # Process the observations
    def __init__(self, om, images_path, nreps, dc, lambda_repeat = 4, verbose = False):

        print(f"Processing images...")

        self.info = {"ObservationMode" : om,
                     "Images_headers" : {}}

        self.data = np.zeros((2, # Number of cameras 
                              cf.om_config[om]["Nlambda"],  # Number of wavelengths
                              cf.om_config[om]["Nmods"],    # Number of Modulations
                              cf.xsize, cf.ysize))  # Size of image (x, y)

        nmods   = cf.om_config[om]["Nmods"]     # N mods from config file
        nlambda = cf.om_config[om]["Nlambda"]   # N wavelengths from config file

        images_path_reshaped = np.array(images_path).reshape(nreps, nlambda, lambda_repeat, nmods, 2)

        
        for rep in range(nreps):
            for lambd in range(nlambda):
                if f"wv_{lambd}" not in self.info["Images_headers"]:
                    self.info["Images_headers"][f"wv_{lambd}"] = {}
                for lambd_rep in range(lambda_repeat):

                    for mod in range(nmods):
                        if f"Mod_{mod}" not in self.info["Images_headers"][f"wv_{lambd}"]:
                            self.info["Images_headers"][f"wv_{lambd}"][f"Mod_{mod}"] = {}
                        # Reading each image
                        im0, head0 = read(images_path_reshaped[rep, lambd, lambd_rep, mod, 0]) # Cam 1
                        im1, _ = read(images_path_reshaped[rep, lambd, lambd_rep, mod, 1]) # Cam 2
                        # Saving images header except for CameraID entry
                        for key in head0:
                            if key == "cam":
                                pass
                            else:
                                if key not in self.info["Images_headers"][f"wv_{lambd}"][f"Mod_{mod}"]:
                                    self.info["Images_headers"][f"wv_{lambd}"][f"Mod_{mod}"][key] = []
                                self.info["Images_headers"][f"wv_{lambd}"][f"Mod_{mod}"][key].append(head0[key])
            
                        # Sving image data into main data array
                        self.data[0, lambd, mod] += im0 - (dc[0] * head0["nAcc"])
                        self.data[1, lambd, mod] += np.flip(im1, axis = -1) - (dc[1] * head0["nAcc"]) # Flip cam 2 image. 
        
        self.data /= (nreps * lambda_repeat)
        self.data /= (nreps * lambda_repeat)

        # Completing info of Observation Mode with info from header
        self.info["nAcc"] = head0["nAcc"]
        self.info["Roix"] = head0["Roix"]
        self.info["Roiy"] = head0["Roiy"]
        self.info["Roix_offset"] = head0["Roix_offset"]
        self.info["Roiy_offset"] = head0["Roiy_offset"]
            
        # Saving info from config file into observtaion mode info
        for entry in cf.om_config[head0["ObservationMode"]]:
            self.info[entry] = cf.om_config[head0["ObservationMode"]][entry]

    # Get information of the observation mode and individual headers
    def get_info(self):
        return self.info
    
    # Get the data array
    def get_data(self):
        return self.data

def get_images_paths(queries):

    """
    Queries have to be in format: "DXX-start-end"]
    start and end are integers. 
    "DXX" has to be one of the observation days -> D09 - D16 
    """
    
    "Allowing for various queries in case observation changed day"

    if isinstance(queries, list):
        selection = []
        for qry in queries:
            parsed = qry.split("-")
            day = parsed[0]
            start = int(parsed[1])
            end = int(parsed[2])

            if end < start:
                raise Exception(f"Query : {qry} not valid. Please prove end of quey larger than start.")
            if day not in ["D09", "D10", "D11","D12","D13","D14","D15","D16"]:
                raise Exception(f"Query : {qry} not valid. Please prove a day within the list: D09, D10, D11, D12, D13,D14, D15, D16")
            df = pd.read_csv(f"{Organization_folder_files}/{day}.csv", index_col=False, header=None)
            selection_df = df[(df.iloc[:, 0] >= start) & (df.iloc[:, 0] <= end)]
            selection.append(selection_df.iloc[:, 1].tolist())

        selection = selection[0] + selection[1] # Concatenating both lists.
        
    else:
        parsed = queries.split("-")
        day = parsed[0]
        start = int(parsed[1])
        end = int(parsed[2])

        if end < start:
            raise Exception(f"Query : {queries} not valid. Please prove end of quey larger than start.")
        if day not in ["D09", "D10", "D11","D12","D13","D14","D15","D16"]:
            raise Exception(f"Query : {queries} not valid. Please prove a day within the list: D09, D10, D11, D12, D13,D14, D15, D16")
        df = pd.read_csv(f"{Organization_folder_files}/{day}.csv", index_col=False, header=None)
        selection_df = df[(df.iloc[:, 0] >= start) & (df.iloc[:, 0] <= end)]
        selection = selection_df.iloc[:, 1].tolist()

    return selection
        
def read_ID(image_index, plotflag = False, verbose = False, header = False, binning = False):
    """
    Reads an image with a given ID.
    Inputs: 
        - image_index : format DXX-YYYY
        - plotflag : Boolean to select plotting
        - verbose : Boolean info on terminal
        - header : Boolean to print header. 
        - Binning : Boolean to bin for faster reading
    Outputs:
        - I : Image data
        - H : Header data
    """
    
    day = image_index[:3]
    index = int(image_index[4:])

    def bin_image(image, bin_size = 4):
    
        Nx, Ny = image.shape
        # Reshape and bin the image by averaging over bin_size x bin_size blocks
        binned_image = image.reshape(Nx // bin_size, bin_size, Ny // bin_size, bin_size).mean(axis=(1, 3))
        
        return binned_image

    df = pd.read_csv(f"{Organization_folder_files}/{day}.csv")
    row = df[df.iloc[:, 0] == index]
    print(row.iloc[0, 1])
    I, H = read(row.iloc[0, 1])

    if binning:
        I = bin_image(I)

    if verbose:
        print("OC", H["ObservationCounter"])
        print("OM", H["ObservationMode"])

    if header:
        print(H)

    if plotflag:
        plt.figure(figsize = (10, 10))
        plt.imshow(I, cmap = "gray")
        plt.show()

    return I, H

def separate_ocs(paths, verbose = True, flat_fieldmode = False):

    """
    
    
    """


    print(f"\nSeparating Observation counters...")
    tic = time.time()
    OCs = {}

    if flat_fieldmode:
        mult = 4
    else:
        mult = 1

    completed_ocs = []

    def labelling(oc, completed_ocs):
        annex = 0
        flag = True
        label = oc
        while flag:
            if label not in completed_ocs:
                flag = False
            else:
                label = f"{oc}_{annex}"
                annex += 1
        return label

    for ind, im in enumerate(paths):
        
        print(f"{ind}/{len(paths)} read.")
        _, H = read(im)

        oc = H['ObservationCounter']
        oc_ind = labelling(oc, completed_ocs)

        if oc_ind not in OCs:
            OCs[oc_ind] = {}
            OCs[oc_ind]["OM"] = H["ObservationMode"]
            OCs[oc_ind]["ims"] = []
            OCs[oc_ind]["empty"] = True
            OCs[oc_ind]["ims"].append(im)

            if H["ObservationMode"] in cf.om_config:
                OCs[oc_ind]["Expected Nim"] = cf.om_config[H["ObservationMode"]]["images_per_mode"]  * mult
            else:
                OCs[oc_ind]["Expected Nim"] = 999

        elif OCs[oc_ind]["empty"]:
            OCs[oc_ind]["ims"].append(im)
            if len( OCs[oc_ind]["ims"]) == OCs[oc_ind]["Expected Nim"]:
                OCs[oc_ind]["empty"] = False
                completed_ocs.append(oc_ind)

    if verbose:
        ocs = [x for x in OCs]
        print(f"Images processed in {round(time.time() - tic, 2)}s.")
        print(f"{len(ocs)} found.")
        print(f"Images por mode:")
     
        for OC in OCs:
            if OCs[OC]["empty"]:
                state = "Incomplete"
            else:
                state = "Complete"
                
            print(f"OC : {OC} - Obs Mode : {OCs[OC]['OM']} - Nims : {len(OCs[OC]['ims'])} - {state}")

    return OCs

def get_time_from_filename(filename):
    split = [int(x) for x in filename[:-4].split("_")]
    return datetime(split[0], split[1], split[2], split[3], split[4], split[5])

def obs_mode_separator(paths, verbose = False):

    obs = {}
    for ind, im in enumerate(paths):
        print(f"Ordering ims - {ind}/{len(paths)}")
        _, H = read(im)

        oc = H["ObservationCounter"]
        wave = H["hvps_comm_volts"]
        mod = f"{H['lcvr1_counts']}_{H['lcvr2_counts']}"

        if oc not in obs:
            obs[oc] = {}

        if wave not in obs[oc]:
            obs[oc][wave] = {}

        if mod not in obs[oc][wave]:
            obs[oc][wave][mod] = []
        
        obs[oc][wave][mod].append(im)

    if verbose:
        for oc in obs:
            print(f"\n Obs Count: {oc}\n-------------")
            for wave in obs[oc]:
                print(f"\n Wave: {wave}")
                for mod in obs[oc][wave]:
                    print(f"   - Mod: {mod} -> Nims : {len(obs[oc][wave][mod])}")

    return obs

def check_timestamps(paths, verbose = True):

    intervals = []

    prev = get_time_from_filename(os.path.basename(paths[0]))

    times = []

    for ind, file in enumerate(paths[1:]):

        filename = os.path.basename(file)

        time = get_time_from_filename(filename)

        times.append(time)

        intervals.append(time - prev)

        prev = time

    intervals = [ x.total_seconds() for x in intervals]

    fig, axs = plt.subplots(figsize = (10, 5))

    axs.plot(times, intervals, c = 'indigo', lw = 1)

    axs.grid(True, c = 'k', alpha = 0.3)

    axs.set_ylabel("Interval between consecutive images [s]")
    axs.set_xlabel("Time of image accquisition.")
    plt.tight_layout()
    plt.show()

def snapshot_processing(paths, dc, verbose = True):

    if verbose:
        print(f"Proccessing snapshot mode...")

    c1 = []
    c2 = []
    for ind, im in enumerate(paths):

        I, H = read(im)

        # Correct dark-current and save image data into main data array
        if H["cam"] == 0:
            c1.append(I - dc[0] * H["nAcc"])
        else:
            c2.append(np.flip(I, axis = -1) - dc[1] * H["nAcc"])# Flip cam 2 image. 

    return np.array([c1, c2])

def polarizers_parser(paths, filt, dc):

    Nfilts = 3
    Nmods = 4
    Ncams = 2

    if filt == "525.02":
        nfilt = 0
    elif filt == "517":
        nfilt = 1
    elif filt == "525.06":
        nfilt = 2
    elif filt == "all":
        nfilt = -1
    else:
        raise Exception("Please provide filt string within: 525.02 / 517 / 525.06 for the micropols")

    reshaped = np.array(paths).reshape(Nfilts, Nmods, Ncams)

    if nfilt >= 0:
        micros = np.zeros((Ncams, Nmods, 2016, 2016))    
        for mod in range(Nmods):
            for cam in range(Ncams):
                print(f"Mod: {mod} - Cam: {cam}")
                im0, H = read(reshaped[nfilt, mod, 0]) # Cam 1
                im1, _ = read(reshaped[nfilt, mod, 1]) # Cam 2
                
                micros[0,  mod] = im0 - (dc[0] * H["nAcc"])
                micros[1,  mod] = np.flip(im1, axis = -1) - (dc[1] * H["nAcc"]) # Flip cam 2 image. 
    else:
        micros = np.zeros((Ncams, Nfilts, Nmods, 2016, 2016))   
        for filt in range(Nfilts): 
            for mod in range(Nmods):
                for cam in range(Ncams):
                    print(f"Mod: {mod} - Cam: {cam}")
                    im0, H = read(reshaped[filt, mod, 0]) # Cam 1
                    im1, _ = read(reshaped[filt, mod, 1]) # Cam 2
                    
                    micros[0, filt, mod] = im0 - (dc[0] * H["nAcc"])
                    micros[1, filt, mod] = np.flip(im1, axis = -1) - (dc[1] * H["nAcc"]) # Flip cam 2 image. 
                
    return micros
        
    



