# ---------------------------- DESCRIPTION --------------------------------------- #

"""
File with all the pre-established configurations for the observations. 

author: Pablo Santamarina Guerrero (pablosantamarinag@gmail.com) 
Instituto de Astrofísica de Andalucía (IAA-CSIC) 
"""

# ------------------------------ CONFIG ------------------------------------------ #

# Exposition time
tauexp = 42 # ms

# Etalon tuning constant
tuning_ct = {"517" : 3.1E-4, # A/V
             "525" : 2.9E-4  # A/V
            }

# Observation Modes configuration:
om_config  = {
    "0s" : {"Nlambda" : 15,
            "line" : "517",
            "lambda_array" : [-400, -300, -200, -100, 0, 100, 200, 300, 400, 500, 600, 650],
            "V_array" : [-3576, -3243, -2910, -2577, -2244, -1911, -1578, -1245, -912, -579, -246, -80],
            "lvcr_mode" : "vectorial"},
    "1" : {"Nlambda" : 10,
           "line" : "517",
            "lambda_array" : [-300, -200, -100, -50, 0, 50, 100, 200, 300, 650],
            "V_array" : [-3243, -2910, -2577, -2411, -2244, -2078, -1911, -1578, -1245, -80],
            "lvcr_mode" : "vectorial"},
    "2.02" : {"Nlambda" : 8,
              "line" : "525.02",
              "lambda_array" : [-120, -80, -40, 0, 40, 80, 120, 227],
              "V_array" : [1729, 1863, 1996, 2129, 2262, 2395, 2529, 2885],
              "lvcr_mode" : "vectorial"},
    "2.06" : {"Nlambda" : 8,
              "line" : "525.06",
              "lambda_array" : [-120, -80, -40, 0, 40, 80, 120, 227],
              "V_array" : [-2907, -2773, -2640, -2507, -2374, -2241, -2107, -1751],
              "lvcr_mode" : "vectorial"},
    "3.02" : {"Nlambda" : 5,
              "line" : "525.02",
              "lambda_array" : [-80, -40, 40, 80, 227],
              "V_array" : [1863, 1996, 2262, 2395, 2885],
              "lvcr_mode" : "longitudinal"},
    "3.06" : {"Nlambda" : 5,
              "line" : "525.06",
              "lambda_array" : [-80, -40, 40, 80, 227],
              "V_array" : [-2773, -2640, -2374, -2241, -1751],
              "lvcr_mode" : "longitudinal"},
    "4" : {"Nlambda" : 3,
           "line" : "525.06",
           "lambda_array" : [-80, 0, 80],
           "V_array" : [-2510, -2244, -1978],
           "lvcr_mode" : "longitudinal"},
    "5.02" : {"Nlambda" : 3,
              "line" : "525.02",
              "lambda_array" : [-40, 0, 40],
              "V_array" : [1996, 2129, 2262],
              "lvcr_mode" : "vectorial"},
    "5.06" : {"Nlambda" : 3,
              "line" : "525.06",
              "lambda_array" : [-40, 0, 40],
              "V_array" : [-2640, -2507, -2374],
              "lvcr_mode" : "vectorial"},









}



