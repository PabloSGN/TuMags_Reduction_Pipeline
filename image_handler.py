# ---------------------------- DESCRIPTION --------------------------------------- #

"""
Information of the headers of the .img files and header creation of reduced
observations.
"""

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in Libs
from astropy.io import fits

# Own libs
from utils import read_Tumag

# ------------------------------ CONFIG ------------------------------------------ #

# RAW Header information 
filter_wheel = {
    "1" : { # First filter wheel
            0 : "Micropolarizers",
            1 : "Linear Polarizer",
            2 : "F4",
            3 : "Phase Diversity",
            4 : "Pinholes"
        },
    "2" : { # Second filter wheel
            0 : "517",
            1 : "525.02",
            2 : "525.06",
            3 : "Empy",
            4 : "Dummy" 
        }
    } 

# Observation modes and inex correspondance
observation_modes = {
    0 : "0s",
    1 : "0p",
    2 : "1",
    3 : "2.02",
    4 : "2.06",
    5 : "3.02",
    6 : "3.06", 
    7 : "4",
    8 : "5.02",
    9 : "5.06",
    64 : "PD_calibration",
    65 : "Spectral_calibration",
    66 : "Polarimetric_calibration"}

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
            "ObservationMode" : observation_modes[om],
            "nAcc" : nAcc,
            "Roix" : Roix,
            "Roiy" : Roiy,
            "Roix_offset" : Roix_offset,
            "Roiy_offset" : Roiy_offset,
            "ObservationCounter" : oc,
            "FW1_ind" : fw1,
            "FW2_ind" : fw2,
            "FW1" : filter_wheel["1"][fw1],
            "FW2" : filter_wheel["2"][fw2],
            "hvps_counts" : hvps_counts,
            "hvps_sign" : hvps_sign,
            "lcvr1_counts" : lcvr1_counts,
            "lcvr2_counts" : lcvr2_counts,
            "hvps_read_counts" : hvps_read_counts,
            "lcvr1_read_counts" : lcvr1_read_counts,
            "lcvr2_read_counts" : lcvr2_read_counts,
        }
        print(self.info["hvps_counts"], self.info["hvps_sign"])
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