import numpy as np
import scipy
from pathlib import Path
from scipy.special import erf
from packaging.version import parse
if parse(scipy.__version__) >= parse("1.7"):
    from scipy.integrate import simpson as simps
else:
    from scipy.integrate import simps

MODULE_DIR = Path(__file__).resolve().parent
cal_file = MODULE_DIR / "CALDATA" / f"tumag_calibration_data.npz"
load_ = np.load(cal_file)
blueshift = load_["blueshift"]
reflectividad = load_["reflectividad"]

def read_grid(grid_file):

    with open(grid_file) as f:
        for line in f:
            if ":" not in line:
                continue
            try:
                _, rhs = line.split(":")
                ini, step, end = [
                    float(x.strip().replace(",", ""))
                    for x in rhs.split()
                ]
                n = int(round((end - ini)/step)) + 1
                return np.linspace(ini, end, n)
            except Exception:
                pass
    raise ValueError("No se encontró una definición válida de grid")

grid = read_grid(MODULE_DIR/ "CALDATA" / 'wave.grid')
central_wave = 5172.6843
scan_wls_A=central_wave-np.array([-300, -200, -100, -50, 0, 50, 100, 200, 300, 650])*1e-3,

def single():
    cfg = Config["517"]

    trans = generate_transmission(
        blueshift[100,100],
        reflectividad[100,100],
        scan_wls_A = scan_wls_A,
        wavelengths_A=central_wave-grid*1e-3,
        config=cfg
    )
    return trans

def save_profile(trans,grid, scan_wls_A, filename):
    np.savez(filename, trans=trans, grid=grid, scan_wls_A=scan_wls_A)

Config={
    '517':{
        'wls_norm':0,'Pend':0.0003262,
        'Ord': 5173.4722,
        'pref_b':5172.648,'Min_wvl':5170,
        'Max_wvl':5176,'R':0.898,'n':2.60,'d':281e-6,'theta':0,
        'pref_c':0.5}, #0.5},
    '525.02':{
        'wls_norm':-5,'Pend':0.0002957121398329138,
        'Ord':5249.543594995222,'pref_b':5250.5,'Min_wvl':5246,
        'Max_wvl':5255,'R':0.77,'n':2.56,'d':281e-6,'theta':0,
        'pref_c':0.65},
    '525.06':{
        'wls_norm':7,'Pend':0.000288733333332857,
        'Ord':5251.371833333332,'pref_b':5250.5,'Min_wvl':5246,
        'Max_wvl':5255,'R':0.75,'n':2.56,'d':281e-6,'theta':0,
        'pref_c':0.65}
}


# Índice refracción: 2.50
# Espesor del etalón: 281 micras
# Reflectividad efectiva (del doble paso): 0.90
# Constante sintonización: 0.0003311 Angstrom/V
# Offset del ajuste de la constante: 5173.374645886531 Angstrom
# Longitud central del prefiltro: 517.265 nm
# Sigma del prefiltro (Gaussiano): 0.05 nm


def finesse(R):
    return 4*R/(1-R)**2

def tune_thickness(n,d,theta,l0_m):
    m=round(2*n*d*np.cos(theta)/l0_m)
    return d+(m*l0_m-2*n*d*np.cos(theta))/(2*n)


def etalon(wavelengths_A,l0_A,R,n,d,theta):
    wl=wavelengths_A*1e-10
    l0=l0_A*1e-10
    d_t=tune_thickness(n,d,theta,l0)
    F=finesse(R)
    opd=n*d_t*np.cos(theta)
    return 1/(1+F*np.sin(2*np.pi*opd/wl)**2)

# def prefilter(wavelengths_A,a,b,c):
#     return a*np.exp(-(wavelengths_A-b)**2/(2*c**2))
def prefilter(wavelengths_A, a, b, c, alpha=0.0):
    return prefilter_supergaussian(
        wavelengths_A,
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        p=2.25
    )

def prefilter_supergaussian(wavelengths_A, a, b, c, alpha=0.0, p=2.0):
    """
    Supergaussiana con asimetría (skew).

    wavelengths_A : eje en Angstrom
    a             : amplitud
    b             : centro en Angstrom
    c             : anchura tipo sigma en Angstrom
    alpha         : parámetro de asimetría (0 = simétrica, >0 cola dcha, <0 cola izq)
    p             : exponente. p=2 es gaussiana normal.
                    Para el 517 PC#2 de la Excel, p~2.25 va mejor.
    """
    # 1. Componente supergaussiana base
    base_supergauss = a * np.exp(-0.5 * np.abs((wavelengths_A - b) / c)**p)
    
    # 2. Factor de asimetría usando la función de error
    # Se normaliza con sqrt(2)*c para mantener la consistencia física con la anchura
    skew_factor = 1.0 + erf(alpha * (wavelengths_A - b) / (np.sqrt(2) * c))
    
    return base_supergauss * skew_factor

def transmission(wavelengths_A,l0_A,a,b,c,R,n,d,theta,alpha=0.0):
    return prefilter(wavelengths_A,a,b,c,alpha=alpha) * etalon(wavelengths_A,l0_A,R,n,d,theta)

def integrate_profile(wavelengths_A,spectrum,l0_A,a,b,c,R,n,d,theta,alpha=0.0):

    T = transmission(
        wavelengths_A,
        l0_A,
        a,b,c,
        R,n,d,theta,alpha=alpha
    )

    return simps(
        spectrum*T,
        x=wavelengths_A
    )

def generate_profile(scan_wls_A,wavelengths_A,spectrum,config,pref_amp=1.):

    prof=np.array([
        integrate_profile(
            wavelengths_A,
            spectrum,
            l0,
            pref_amp,
            config["pref_b"],
            config["pref_c"],
            config["R"],
            config["n"],
            config["d"],
            config["theta"]
        )
        for l0 in scan_wls_A
    ])
    # continuum_idx = -1
    # prof /= prof[continuum_idx]
    return prof

def volts_2_lambda(volts, config):
    return config['Pend'] * volts + config['Ord']
    
def generate_transmission(
        blueshift,
        reflectivity,
        scan_wls_A,
        wavelengths_A,
        config):

    t= transmission(
        wavelengths_A,
        scan_wls_A[0],
        1,
        config["pref_b"] + blueshift,
        config["pref_c"],
        reflectivity,
        config["n"],
        config["d"],
        config["theta"]
    )
    trs = np.empty((len(scan_wls_A), len(wavelengths_A)))
    for index,l0 in enumerate(scan_wls_A):
        print(l0,scan_wls_A[0])
        trs[index] = transmission(
            wavelengths_A,
            l0,
            1.0,
            config["pref_b"] + blueshift,
            config["pref_c"],
            reflectivity,
            config["n"],
            config["d"],
            config["theta"]
        )
        

    return trs

def inspect_scan_position(
    wavelengths_A,
    spectrum,
    l0_A,
    config,
    pref_amp=1.,
):

    E = etalon(
        wavelengths_A,
        l0_A,
        config["R"],
        config["n"],
        config["d"],
        config["theta"]
    )

    PF = prefilter(
        wavelengths_A,
        pref_amp,
        config["pref_b"],
        config["pref_c"]
    )

    T = PF*E

    weighted_spectrum = spectrum*T

    num = simps(weighted_spectrum,x=wavelengths_A)
    den = simps(T,x=wavelengths_A)

    integrated_value = num/den

    return {
        "wavelengths":wavelengths_A,
        "input_spectrum":spectrum,
        "etalon":E,
        "etalon2":E,
        "prefilter":PF,
        "transmission":T,
        "weighted_spectrum":weighted_spectrum,
        "integrated_value":integrated_value
    }

def apply_transmission(data_tumag,blueshift,reflectividad,grid,Config):
    central_wave = 5172.6843
    cfg = Config["517"]

    scan_wls_A=central_wave+np.array([-300, -200, -100, -50, 0, 50, 100, 200, 300, 650])*1e-3
    wavelengths_A=central_wave-grid*1e-3

    data_tumag_conv = np.empty((data_tumag.shape[0], len(scan_wls_A), data_tumag.shape[-2], data_tumag.shape[-1]))

    nx = data_tumag.shape[-1]
    ny = data_tumag.shape[-2]
    for i in range(nx):
        print(i)
        for j in range(ny):

            T = generate_transmission(
                blueshift[i,j],
                reflectividad[i,j],
                scan_wls_A,
                wavelengths_A,
                cfg
            )
            for stokes in range(4):
                for l in range(len(scan_wls_A)):
                    data_tumag_conv[stokes,l,i,j] = simps(
                        data_tumag[stokes,:,i,j] * T[l]
                        )

    data_tumag_conv = data_tumag_conv / np.mean(data_tumag_conv[0,-1,0:100,0:100])

def example(data_tumag,blueshift,reflectividad,grid,Config):

    # usando pares unicos

    central_wave = 5172.6843
    cfg = Config["517"]

    scan_wls_A = central_wave + np.array([-300, -200, -100, -50, 0, 50, 100, 200, 300, 650]) * 1e-3
    wavelengths_A = central_wave + grid * 1e-3

    data_tumag_conv = np.empty((data_tumag.shape[0], len(scan_wls_A),data_tumag.shape[-2], data_tumag.shape[-1]))
    data_tumag_conv_clean = np.empty_like(data_tumag_conv)

    nx = data_tumag.shape[-1]
    ny = data_tumag.shape[-2]

    # =====================================================
    # Discretización 50x50
    # =====================================================

    NBS = 50
    NREF = 50

    bs_edges = np.linspace(blueshift.min(), blueshift.max(), NBS + 1)
    ref_edges = np.linspace(reflectividad.min(), reflectividad.max(), NREF + 1)

    bs_idx = np.clip(np.digitize(blueshift, bs_edges) - 1, 0, NBS - 1)
    ref_idx = np.clip(np.digitize(reflectividad, ref_edges) - 1, 0, NREF - 1)

    # pares únicos realmente presentes
    unique_pairs = np.unique(np.column_stack((bs_idx.ravel(), ref_idx.ravel())),axis=0)

    print(f"{len(unique_pairs)} pares únicos de un máximo de {NBS*NREF}")

    # =====================================================
    # Calcula T una sola vez por par único
    # =====================================================

    T_dict = {}

    for ibs, iref in unique_pairs:
        bs_center = 0.5 * (bs_edges[ibs] + bs_edges[ibs + 1])
        ref_center = 0.5 * (ref_edges[iref] + ref_edges[iref + 1])
        T_dict[(ibs, iref)] = generate_transmission( bs_center,ref_center,scan_wls_A,wavelengths_A,cfg)

    # =====================================================
    # Convoluciones
    # =====================================================

    nn = 0
    for i in range(nx):
        print(i)
        for j in range(ny):
            T = T_dict[(bs_idx[i, j],ref_idx[i, j])]
            if nn == 0:
                nn = np.sum(T[l])
            for stokes in range(4):
                profile = data_tumag[stokes, :, i, j]
                for l in range(len(scan_wls_A)):
                    # Tn = T[l] / np.sum(T[l])
                    Tn = T[l] / nn
                    data_tumag_conv[stokes, l, i, j] = simps(profile * Tn)
                    Tnc = np.where(np.abs(wavelengths_A - scan_wls_A[l]) >= 1.0,0,T[l])
                    #Tnc /= np.sum(Tnc)
                    Tnc /= nn
                    data_tumag_conv_clean[stokes, l, i, j] = simps(profile * Tnc)
    