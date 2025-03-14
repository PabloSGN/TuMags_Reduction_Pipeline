# ---------------------------- DESCRIPTION --------------------------------------- #

# ------------------------------ IMPORTS ----------------------------------------- #

# Built-in 
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle

from functools import partial

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mp
from astropy.io import fits
from scipy.optimize import minimize
from scipy.integrate import simpson
from scipy.interpolate import interp1d


# Own-libs
import image_handler as ih
import config as cf

# ------------------------------ CONFIG ------------------------------------------ #

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

# Etalon params 
"""Et = {
    'R' : 0.75,
    'n' : 2.56,
    'd' : 281e-6, #0.0003159677641232576,
    'theta' : 0}"""


Config = {
'517' : {'Gamma'    : 0.47,
         'wls_norm' : 0,
         'Pend'     : 0.00030907042253499933,
         'Ord'      : 5173.432608450703,
         'pref_b'    : 5172.5,
         'Min_wvl'  : 5170,
         'Max_wvl'  : 5176,
         'R'        : 0.81,
         'n'        : 2.56,
         'd'        : 281e-6,
         'theta'    : 0,
         "pref_c" : 0.5},
'525.02': {
        'Gamma'    : 0.98,
        'wls_norm' : -5,
        'Pend'     : 0.0002957121398329138,
        'Ord'      : 5249.543594995222,
        'pref_b'    : 5250.5,
        'Min_wvl'  : 5246,
        'Max_wvl'  : 5255,
        'R'        : 0.77,
        'n'        : 2.56,
        'd'        : 281e-6,
        'theta'    : 0,
        "pref_c" : 0.65},
'525.06' : {
        'Gamma'    : 1,
        'wls_norm' : 7,
        'Pend'     : 0.000288733333332857,
        'Ord'      : 5251.371833333332,
        'pref_b'    : 5250.5,  
        'Min_wvl'  : 5246,
        'Max_wvl'  : 5255,
        'R'        : 0.75,
        'n'        : 2.56,
        'd'        : 281e-6,
        'theta'    : 0,
        "pref_c" : 0.65 }}

# ------------------------------  AUX FUNCTIONS  --------------------------------- # 


def Finesse(R, dR):
    return 4 * (R * dR) / (1 - R * dR) ** 2

def thickness_tunning(n, d, theta, l0):
    
    m  = round(2 * n * d * np.cos(theta) / l0)   # Order of the resonance peak
    dh = (m * l0 - 2 * n * d * np.cos(theta)) / (2 * n)          # Variation of thickness to tune again to wvl0
    thick = d + dh

    return thick

def Etalon(lambd, l0, R, dR, n, d, theta):

    thickness = thickness_tunning(n, d, theta, l0)
    a = n * thickness * np.cos(theta)

    return 1 / (1 + Finesse(R, dR) * np.sin(2 * np.pi * a / lambd) ** 2)

def fts_spectra(wli, wlf):

    """'Kitt Peak FTS-Spectral-Atlas'"""
    # print('Enter end wavelength (3290 - 12508 A)')
    #path = '/Users/dorozco/Dropbox (IdAdA)/Python/'
    file = os.path.join(__location__, 'fts.npz')
    # np.savez('fts.npz', fts=fts, fts_w=fts_w)
    data = np.load(file)
    fts = data['fts']
    fts_w = data['fts_w']

    indx = np.where((fts_w > int(wli)) & (fts_w < int(wlf)))

    return fts[indx], fts_w[indx]

def prefilter(wavelengths, a, b, c):
    exp = np.exp(-((wavelengths - b) ** 2 / (2 * c ** 2)))    
    return a * exp

def Int(Wavelengths, l0, Spectrum, a, b, c, R, dR, n, d, theta):
  
    I = simpson(Spectrum * prefilter(Wavelengths, a, b, c) * Etalon(Wavelengths * 1E-10, l0 * 1E-10,  R, dR, n, d, theta) ** 2, x = Wavelengths)
    
    return I

def Prof(Wls, Wavelengths, Spectrum, a, b, c, R, dR, n, d, theta, wl_norm, Gamma):
    
    Scan = [Int(Wavelengths, wl, Spectrum, a, b, c, R, dR, n, d, theta) for wl in Wls] 
    
    Scan = np.array(Scan)
    
    Scan = Gamma * (Scan / Scan[wl_norm])
    
    return Scan 

def diff(I, Wvls, wavelengths, Spectrum, a, b, c, gamma, R, n, d, theta, config):

    diff_profs = (I - Prof(Wvls, Wavelengths = wavelengths, Spectrum=Spectrum, 
                           a=a, b=b, c=c, R=R, dR=1, n = n, d=d, 
                           theta= theta, wl_norm=config["wls_norm"], Gamma = gamma))

    return np.sum([x** 2 for x in diff_profs])

def compute_profile(wls, wavelengths, Spectrum, a, b, c, gamma, R, n, d, theta, config):
    
    return Prof(wls, Wavelengths = wavelengths, Spectrum=Spectrum, a=a, b=b, c=c, R=R, dR=1, n=n, d=d, theta = theta, wl_norm=config["wls_norm"], Gamma = gamma)


def parse_prefilter_scan(images_paths, filt, verbose = True, cam = 0):

    if verbose:
        print(f"Parsing prefilter scan images...")

    volts = []
    intensity = []

    found_filters = []
    found_cams = []
    for ind, img_path in enumerate(images_paths):
        I, H = ih.read(img_path)

        if H["FW2"] not in found_filters:
            found_filters.append(H["FW2"])
        if H["cam"] not in found_cams:
            found_cams.append(H["cam"])
        if H["cam"] == cam and H["FW2"] == filt:
            volts.append(H["hvps_read_volts"])
            intensity.append(np.mean(I[500:-500, 500:-500]))
        else:
            pass

    if len(volts) < 1:
        raise Exception(f"No images found for filt {filt} and cam:{cam}.\nPlease select one of:\nFound cameras : {found_cams}\nFound filters : {found_filters}")

    volts = np.array(volts)
    intensity = np.array(intensity)

    sorted_indices = np.argsort(volts)  # Get sorting indices based on arr1
    volts_sorted = volts[sorted_indices]
    intensity_sorted = intensity[sorted_indices]

    return volts_sorted, intensity_sorted

# ------------------------------ MAIN CODE ------------------------------------------ #


def loss_r_n_b_c(params, wvls_measured, int_measured, wavelengths_spectrum, spectrum, 
               a, gamma, d, theta, config):
    R, n, b, c = params
    guess  = compute_profile(wvls_measured, wavelengths_spectrum, spectrum, a, b, c, gamma, R, n, d, theta, config)
    error = np.sum((guess - int_measured) ** 2)  # Sum of squared errors
    return error


def volts_2_lambda(volts, config):
    return config['Pend'] * volts + config['Ord']
    
def lambda_2_volts(lambd, config):
    return (lambd - config["Ord"]) / config["Pend"]

def fit_prefilter(V, I, filt, save_flag = False, bounds = None, plot_flag = False, plot_filename = "Prefilter_fitting.png"):

    a = 0.8

    config = Config[filt]

    # Normalization
    I /= I[config['wls_norm']]
    I *= config['Gamma']

    Spectrum, wavelengths = fts_spectra(config['Min_wvl'], config['Max_wvl'])
    Spectrum = Spectrum / 10000 # Continuum normalization
    Wvls = volts_2_lambda(np.array(V), config)

    G = config["Gamma"]

    initial_guess = [config["R"], config["n"], config["pref_b"], config["pref_c"]]

    if bounds is None:
        bounds = [(0.6, 0.81), 
                  (2.4, 2.7),
                  (config["pref_b"] - 0.1, config["pref_b"] + 0.1),
                  (0.4, 0.6)]

    result = minimize(loss_r_n_b_c, initial_guess,
                      args=(Wvls, I, wavelengths, Spectrum, a, G, config["d"], config["theta"], config), 
                      method='L-BFGS-B', bounds=bounds)

    fitted_r = result.x[0]
    fitted_n = result.x[1]
    fitted_b = result.x[2]
    fitted_c = result.x[3]

    print(f"FITTED MODEL : \nR: {fitted_r}, n: {fitted_n}, b: {fitted_b}, c: {fitted_c}")

    whole_volts_range = np.linspace(-4000, 4000, 1000)

    fitted_prefilter = compute_profile(volts_2_lambda(whole_volts_range, config), wavelengths,
                                       np.ones(len(wavelengths)), a, fitted_b, fitted_c, 
                                       G, fitted_r, fitted_n, config["d"], config["theta"],
                                       config)

    fitted_model = interp1d(whole_volts_range, fitted_prefilter / np.max(fitted_prefilter), kind='cubic')  

    if save_flag:
        # Store necessary data
        data = {'x': whole_volts_range, 'y': fitted_prefilter / np.max(fitted_prefilter)}
        # Save using pickle
        with open('prefilter_model.pkl', 'wb') as file:
            pickle.dump(data, file)

    if plot_flag:
        fig, axs = plt.subplots(figsize = (13, 5))
        axs.plot(wavelengths, Spectrum, c = 'k', lw = 2)
        axs.plot(wavelengths, prefilter(wavelengths, a, fitted_b, fitted_c), c = 'dodgerblue', 
                 lw = 2, label = "Fitted Prefilter")
        axs.plot(Wvls, I, c = "darkorange", ls ='', marker = ".", markersize = 10, label = 'Measure')
        axs.plot(Wvls, compute_profile(Wvls, wavelengths, Spectrum, a, fitted_b,  fitted_c, G, 
                                       fitted_r, fitted_n, config["d"], config["theta"],
                                       config), c = 'indigo', lw = 3, label = "Fitted profile")

        pref_effect = compute_profile(Wvls, wavelengths, np.ones(len(wavelengths)), a, fitted_b,
                                       config["pref_c"], G, fitted_r, fitted_n, config["d"], 
                                       config["theta"], config) 
        pref_effect /= np.max(pref_effect)
        axs.plot(Wvls, pref_effect, c = 'forestgreen', lw = 3, label = "Fitted prefilter effect [Normalized]")
        axs.plot(Wvls, compute_profile(Wvls, wavelengths, Spectrum, a, fitted_b, 
                                       config["pref_c"], G, fitted_r, fitted_n, config["d"], 
                                       config["theta"], config) / pref_effect, 
                                         c = 'deeppink', lw = 3, label = "Corrected profile")
        axs.legend(edgecolor = 'k')
        axs.set_ylim(-0.05, 1.15)
        axs.set_xlim(min(Wvls) - 0.2, max(Wvls) + 0.2)
        axs.grid(True, c = 'k', alpha = 0.2)
        plt.tight_layout()
        fig.savefig(plot_filename, bbox_inches="tight")
        #plt.close(fig)

    return fitted_model

def load_model(filename):
    with open(filename, 'rb') as file:
        data = pickle.load(file)

    model = interp1d(data['x'], data['y'], kind='cubic')

    return model

def remove_line_from_flat_fields(ff, om, pref_model, volts = None, verbose = False):

    if verbose:
        print(f"\nRemoving prefilter effect from flat-fields...")

    if isinstance(pref_model, str):
        pref_model = load_model(pref_model)

    shape = np.shape(ff)
    nlambda = shape[1]

    if volts is None:
        volts = cf.om_config[om]["V_array"]

    new_flat_field = np.zeros(shape)

    for lambd in range(nlambda):
        factor = pref_model(volts[lambd])
        new_flat_field[:, lambd] = ff[:, lambd] * factor

    if verbose:
        print(f"done.\n")

    return new_flat_field

