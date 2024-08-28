# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Information of the headers of the .img files and header creation of reduced
observations.

author: Pablo Santamarina Guerrero(pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Own libs
from utils import read_Tumag
import config as cf

# Config
Organization_folder_files = "/home/users/dss/orozco/Tumag/PabloTests/Organized_files"

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
                 hvps_read_counts : int, lcvr1_read_counts : int, lcvr2_read_counts : int) -> None:
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
                      int(h["LCVR2_DN_Real"]) )
    
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

        for lambd in range(nlambda):
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
            
                # Sving image data into main data array
                self.data[0, lambd, mod] = im0 - dc[0]
                self.data[1, lambd, mod] = np.flip(im1, axis = -1) - dc[1] # Flip cam 2 image. 
        
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

class nominal_flat:

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
                for lambd_rep in range(lambda_repeat):
                    for mod in range(nmods):
                        self.info["Images_headers"][f"Mod_{mod}"] = {}
                        # Reading each image
                        im0, head0 = read(images_path_reshaped[rep, lambd, lambd_rep, mod, 0]) # Cam 1
                        im1, _ = read(images_path_reshaped[rep, lambd, lambd_rep, mod, 1]) # Cam 2
                        # Saving images header except for CameraID entry
                        self.info["Images_headers"][f"Mod_{mod}"][f"wave_{lambd}"] = {}
                        for key in head0:
                            if key == "cam":
                                pass
                            else:
                                self.info["Images_headers"][f"Mod_{mod}"][f"wave_{lambd}"][key] = head0[key]
            
                        # Sving image data into main data array
                        self.data[0, lambd, mod] += im0 - dc[0]
                        self.data[1, lambd, mod] += np.flip(im1, axis = -1) - dc[1] # Flip cam 2 image. 
        
        self.data[0, lambd, mod] /= (nreps * lambda_repeat)
        self.data[1, lambd, mod] /= (nreps * lambda_repeat)

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
            print(qry)

            parsed = qry.split("-")
            day = parsed[0]
            start = int(parsed[1])
            end = int(parsed[2])

            if end < start:
                raise Exception(f"Query : {qry} not valid. Please prove end of quey larger than start.")
            if day not in ["D09", "D10", "D11","D12","D13","D14","D15","D16"]:
                raise Exception(f"Query : {qry} not valid. Please prove a day within the list: D09, D10, D11, D12, D13,D14, D15, D16")
            df = pd.read_csv(f"{Organization_folder_files}/{day}.csv")
            selection_df = df[(df.iloc[:, 0] >= start) & (df.iloc[:, 0] <= end)]

            selection.append(selection_df.iloc[:, 1].tolist())
    else:
        parsed = queries.split("-")
        day = parsed[0]
        start = int(parsed[1])
        end = int(parsed[2])

        if end < start:
            raise Exception(f"Query : {queries} not valid. Please prove end of quey larger than start.")
        if day not in ["D09", "D10", "D11","D12","D13","D14","D15","D16"]:
            raise Exception(f"Query : {queries} not valid. Please prove a day within the list: D09, D10, D11, D12, D13,D14, D15, D16")
        df = pd.read_csv(f"{Organization_folder_files}/{day}.csv")
        selection_df = df[(df.iloc[:, 0] >= start) & (df.iloc[:, 0] <= end)]

        selection = selection_df.iloc[:, 1].tolist()

    return selection
        
def read_ID(image_index, plotflag = False, verbose = False, header = False):
    
    day = image_index[:3]
    index = int(image_index[4:])

    df = pd.read_csv(f"{Organization_folder_files}/{day}.csv")
    row = df[df.iloc[:, 0] == index]
    print(row.iloc[0, 1])
    I, H = read(row.iloc[0, 1])

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
