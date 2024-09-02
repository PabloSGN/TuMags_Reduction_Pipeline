
import sys
import os
import glob
import matplotlib.pyplot as plt
import numpy as np


# Own Libs
sys.path.append("/home/users/dss/orozco/Tumag/PabloTests")
import config as cf
from utils import read_Tumag
from field_stop_finder import compute_alignment, apply_fieldstop_and_align_array
from master_dark import compute_master_darks
from master_flatfield import compute_master_flat_field
import image_handler as ih
from demodulation import demodulate
from phase_diversity import pd_observation_parser





pd_indexes = "D10-6194-6513"


pd_observation_parser(pd_indexes)