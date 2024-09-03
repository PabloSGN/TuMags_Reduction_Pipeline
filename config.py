# ---------------------------- DESCRIPTION --------------------------------------- #

"""
File with all the pre-established configurations for the observations. 

author: Pablo Santamarina Guerrero (pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ CONFIG ------------------------------------------ #

# Exposure time
tauexp = 42 # ms
xsize = 2016 # Size in pixels X-dimension
ysize = 2016 # Size in pixels Y-dimension

# Etalon tuning constant
tuning_ct = {"517" : 3.1E-4, # A/V
             "525" : 2.9E-4  # A/V
            }

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

# Observation modes and index correspondance
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

# Observation Modes configuration:
om_config  = {
    "0s" : {"Nlambda" : 15,
            "line" : "517",
            "lambda_array" : [-400, -300, -200, -100, 0, 100, 200, 300, 400, 500, 600, 650],
            "V_array" : [-3576, -3243, -2910, -2577, -2244, -1911, -1578, -1245, -912, -579, -246, -80],
            "Nmods" : 1,
            "lvcr_mode" : "vectorial",
            "images_per_mode" : 15 * 1 * 2},

    "0p" : {"Nlambda" : 15,
            "line" : "517",
            "lambda_array" : [-400, -300, -200, -100, 0, 100, 200, 300, 400, 500, 600, 650],
            "V_array" : [-3576, -3243, -2910, -2577, -2244, -1911, -1578, -1245, -912, -579, -246, -80],
            "Nmods" : 4,
            "lvcr_mode" : "vectorial",
            "images_per_mode" : 15 * 4 * 2},
            
    "1" : {"Nlambda" : 10,
           "line" : "517",
            "lambda_array" : [-300, -200, -100, -50, 0, 50, 100, 200, 300, 650],
            "V_array" : [-3243, -2910, -2577, -2411, -2244, -2078, -1911, -1578, -1245, -80],
            "Nmods" : 4,
            "lvcr_mode" : "vectorial",
            "images_per_mode" : 10 * 4 * 2},

    "2.02" : {"Nlambda" : 8,
              "line" : "525.02",
              "lambda_array" : [-120, -80, -40, 0, 40, 80, 120, 227],
              "V_array" : [1729, 1863, 1996, 2129, 2262, 2395, 2529, 2885],
              "Nmods" : 4,
              "lvcr_mode" : "vectorial",
              "images_per_mode" : 8 * 4 * 2},

    "2.06" : {"Nlambda" : 8,
              "line" : "525.06",
              "lambda_array" : [-120, -80, -40, 0, 40, 80, 120, 227],
              "V_array" : [-2907, -2773, -2640, -2507, -2374, -2241, -2107, -1751],
              "Nmods" : 4,
              "lvcr_mode" : "vectorial",
              "images_per_mode" : 8 * 4 * 2},

    "3.02" : {"Nlambda" : 5,
              "line" : "525.02",
              "lambda_array" : [-80, -40, 40, 80, 227],
              "V_array" : [1863, 1996, 2262, 2395, 2885],
              "Nmods" : 2,
              "lvcr_mode" : "longitudinal",
              "images_per_mode" : 5 * 2 * 2},

    "3.06" : {"Nlambda" : 5,
              "line" : "525.06",
              "lambda_array" : [-80, -40, 40, 80, 227],
              "V_array" : [-2773, -2640, -2374, -2241, -1751],
              "Nmods" : 2,
              "lvcr_mode" : "longitudinal"}, 
              "images_per_mode" : 5 * 2 * 2,

    "4" : {"Nlambda" : 3,
           "line" : "517",
           "lambda_array" : [-80, 0, 80],
           "V_array" : [-2510, -2244, -1978],
           "Nmods" : 4,
           "lvcr_mode" : "vectorial"},
           "images_per_mode" : 3 * 4 * 2,

    "5.02" : {"Nlambda" : 3,
              "line" : "525.02",
              "lambda_array" : [-40, 0, 40],
              "V_array" : [1996, 2129, 2262],
              "Nmods" : 4,
              "lvcr_mode" : "vectorial"},
              "images_per_mode" : 3 * 4 * 2,

    "5.06" : {"Nlambda" : 3,
              "line" : "525.06",
              "lambda_array" : [-40, 0, 40],
              "V_array" : [-2640, -2507, -2374],
              "Nmods" : 4,
              "lvcr_mode" : "vectorial",
              "images_per_mode" : 3 * 4 * 2}
}





