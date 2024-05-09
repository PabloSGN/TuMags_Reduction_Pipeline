# ============================ IMPORTS ====================================== #

# Built-in libs
import os
import time
import glob
import numpy as np
from shutil import copyfile

# Only needed if Fits_Flag is set to true
from astropy.io import fits

# Own libs
import images_reader as imr

# ============================= CONFIG ====================================== #

tic = time.time()





# ============================ LOG READING ================================== #

def LogFinder(Folder : str) -> list:
    """
    This finction identifies the Log files within a folder.

    Parameters
    ----------
    Folder : str
        Path to the data folder containing images and logs.

    Returns
    -------
    valid_log_files : list
        List containing all log files paths.
    """

    # Serach for all files containing "log" in the Data folder
    all_log_files = glob.glob(os.path.join(Folder, '*log*'))
    valid_log_files = [] # variable to store valid logfiles

    # Looping over all files with "log" in the name
    for file in all_log_files:

        # Checking if 'image_name_str' is on the first line of the file
        open_file = open(file, 'r')
        if 'image_name_str' in open_file.readlines()[0]:
            valid_log_files.append(file)

    return valid_log_files

def LogReader(DATA : dict, LogPath : str, data_folder : str):
    """
    Function that given a log File, reads it and stores image names and properties
    in a dictionary (DATA)

    Parameters
    ----------
    DATA : dict
        Dictionary used to store images and info.
    LogPath : str
        Path to the logfile to be read.
    data_folder : str
        Path to the folder that contains the log file.

    Returns
    -------
    None.
    """

    file = open(LogPath, 'r') # Opening Log File
    logName = os.path.basename(LogPath)

    # Looping over all lines
    for line_idx, line in enumerate(file.readlines()):

        line = line.rstrip() # Removing end-line breaks

        # Header identification
        if line_idx == 0:
            header = line.split(',') # Separating columns
            header = [hd.replace(" ", "") for hd in header] # Removing spaces

        # Images Data
        else:
            cols = line.split(',')  # Separating columns
            cols = [cl.replace(" ", "") for cl in cols] # Values without spaces
            name = os.path.basename(cols[0]) # Image name
            DATA[name] = {} # Creating entry in dictionary
            DATA[name]['Log'] = logName # Saving the log name
            DATA[name]['Folder'] = data_folder
            # Looping and storing over all info
            for col, val in enumerate(cols):
                DATA[name][header[col]] = val

    file.close()

# ============================= LOG SORTER ================================== #

def Organizer(Data : dict, Obs_mode : str, Wheel_Info : dict):

    """
    Function that given the dictionary containing all images and their properties,
    sorts them following a pre-determined criteria for each calibration mode.

    Parameters
    ----------
    Data : dict
        Dictionary containing images names and properties found in the logs.
    Obs_mode : str
        "Spectral", "Polarimetric" or "PD" for each of the calibrations.
    Wheel_Info : dict,
        Info of the wheel positions.

    Returns
    -------
    Spec_data / Pol_data / PD_data : dict
        Dictionary containing the sorted images as a function of their properties.

        Format:
            Spectral : DATA -> Filter position -> LCVR Modul -> Etalon Volt -> Cam
            Polarimetric : DATA -> Log -> Filter position -> LCVR Modul -> Etalon Volt -> Cam
            PD : DATA -> Mode(Focused/Defocused) ->  Filter Position -> LCVR Modul -> Cam
    """

    if Obs_mode == 'Spectral':

        Spec_data = {}

        # Looping over all images in the Logs
        for img in Data:

            F1 = Wheel_Info['FW1'][Data[img]['fw1_pos']]
            F2 = Wheel_Info['FW2'][Data[img]['fw2_pos']]
            wheel_pos = F1 + '_' + F2

            # Create entry to store new Mod data
            if not wheel_pos in Spec_data:
                Spec_data[wheel_pos] = {}

            # Mod Voltage LCVR1 + _ + Voltage LCVR2
            Mod = str(int(Data[img]['rocli_lcvr1_dn']) * 0.00179439) + '_' \
                + str(int(Data[img]['rocli_lcvr2_dn']) * 0.00179507)

            # Create entry to store new Mod data
            if not Mod in Spec_data[wheel_pos]:
                Spec_data[wheel_pos][Mod] = {}

            # Etalon's voltage
            Etalon_volts = (-1) ** (int(Data[img]['etalon_volt_sign']) + 1) \
                * (int(Data[img]['etalon_volt_dn']) * 4999 / (2 ** 12 - 1))

            # Crerating entries to store data
            if not Etalon_volts in Spec_data[wheel_pos][Mod]:
                Spec_data[wheel_pos][Mod][Etalon_volts] = {}
                Spec_data[wheel_pos][Mod][Etalon_volts]['0'] = [] # Camera 1
                Spec_data[wheel_pos][Mod][Etalon_volts]['1'] = [] # Camera 2

            Ncamera = Data[img]['camID']
            Spec_data[wheel_pos][Mod][Etalon_volts][Ncamera].append(img) # Storing image name

        return Spec_data

    elif Obs_mode == 'PD':

        PD_data = {}

        # Separating between focused and defocused
        for mode in Data:
            PD_data[mode] = {}

            # Looping over all images in each mode
            for img in Data[mode]:

                F1 = Wheel_Info['FW1'][Data[mode][img]['fw1_pos']]
                F2 = Wheel_Info['FW2'][Data[mode][img]['fw2_pos']]
                wheel_pos = F1 + '_' + F2

                # Create entry to store new Mod data
                if not wheel_pos in PD_data[mode]:
                    PD_data[mode][wheel_pos] = {}

                # Separating images between modulation
                Mod = str(int(Data[mode][img]['rocli_lcvr1_dn']) * 0.00179439) + '_' \
                    + str(int(Data[mode][img]['rocli_lcvr2_dn']) * 0.00179507)

                if not Mod in PD_data[mode][wheel_pos]:
                    PD_data[mode][wheel_pos][Mod] = {}
                    PD_data[mode][wheel_pos][Mod]['0'] = []
                    PD_data[mode][wheel_pos][Mod]['1'] = []

                Ncamera = Data[mode][img]['camID']
                PD_data[mode][wheel_pos][Mod][Ncamera].append(img)

        return PD_data

    elif Obs_mode == 'Polarimetric':

        Pol_data = {}

        # Looping over all images in the Logs
        for img in Data:

            log = Data[img]['Log']

            # Creating log entry in variable
            if not log in Pol_data:
                Pol_data[log] = {}

            # Creating entry for filter wheels position
            F1 = Wheel_Info['FW1'][Data[img]['fw1_pos']]
            F2 = Wheel_Info['FW2'][Data[img]['fw2_pos']]
            wheel_pos = F1 + '_' + F2

            if not wheel_pos in Pol_data[log]:
                Pol_data[log][wheel_pos] = {}

            # Creating entry for Modulation Voltage LCVR1 + _ + Voltage LCVR2
            Mod = str(int(Data[img]['rocli_lcvr1_dn']) * 0.00179439) + '_' \
                + str(int(Data[img]['rocli_lcvr2_dn']) * 0.00179507)

            if not Mod in Pol_data[log][wheel_pos]:
                Pol_data[log][wheel_pos][Mod] = {}

            # Etalon's voltage
            Etalon_volts = (-1) ** (int(Data[img]['etalon_volt_sign']) + 1) \
                * (int(Data[img]['etalon_volt_dn']) * 4999 / (2 ** 12 - 1))

            # Crerating entries to store data
            if not Etalon_volts in Pol_data[log][wheel_pos][Mod]:
                Pol_data[log][wheel_pos][Mod][Etalon_volts] = {}
                Pol_data[log][wheel_pos][Mod][Etalon_volts]['0'] = [] # Camera 1
                Pol_data[log][wheel_pos][Mod][Etalon_volts]['1'] = [] # Camera 2

            # Storing image name
            Ncamera = Data[img]['camID']
            Pol_data[log][wheel_pos][Mod][Etalon_volts][Ncamera].append(img)

        return Pol_data

    else:

        raise Exception('Please enter a valid calibration mode -> ' + \
                        '(Spectral / PD / Polarimetric) (Case sensitive)')

# ========================== IMAGE READING ================================== #

def rebin(arr, size):

    """
    This function rebins an array of a given shape and returns
    the new array with shape given by 'new_shape' (containing the new dimensions).
    The shape of the new array must be a divisor of the old one.
    """

    shape = (size, arr.shape[0] // size,
             size, arr.shape[1] // size)

    return arr.reshape(shape).mean(-1).mean(1)

def ImageReading(paths : list, folder : str, image_size : int, dc : list,
                 flat_field : list, Rebin : dict, FieldStop : dict, Acc_flag = True):

    """
    Function that given a list of paths to images, reads them and computes the
    pertinent calculations.

    Parameters
    ----------
    paths : list
        List containing paths to the images that want to be read.
    image_size : int
        Size of image in pixels.
    dc : np.matrix
        Matrix containg the corresponding dark current.
    flat_field : np.matrix
        Matrix containg the corresponding flat field.
    Rebin : dict
        Dictionary containing boolean variable for rebining and size of the rebinning.
    FieldStop : dict
        Dictionary containing boolean variable and size for FOV selection.
    Acc_flag : bool
        Flag to disable accumulation (defaut is True)
    Returns
    -------
    images : np.array
        Matrix containing data of the images read.
    """

    if Acc_flag:

        if Rebin['Flag']:
            images = np.zeros((Rebin['Size'], Rebin['Size']))

        elif FieldStop['Flag']:
            new_size = FieldStop['Limits'][1] - FieldStop['Limits'][0]
            images = np.zeros((new_size, new_size))

        else:
            images = np.zeros((image_size, image_size))

        for im in paths:

            if Rebin['Flag']:

                if FieldStop['Flag']:

                    # Flat and Dark correction
                    image = imr.image_reader(os.path.join(folder, im), False) - dc
                    image = image / flat_field

                    crop = image[FieldStop['Limits'][0] : FieldStop['Limits'][1],
                                 FieldStop['Limits'][0] : FieldStop['Limits'][1]]

                    binned = rebin(crop, Rebin['Size'])
                    images += binned
                else:
                    image = imr.image_reader(os.path.join(folder, im), False) - dc
                    image / flat_field

                    binned = rebin(image, Rebin['Size'])
                    images += binned
            else:

                if FieldStop['Flag']:
                    image = imr.image_reader(os.path.join(folder, im), False) - dc
                    image = image / flat_field

                    crop = image[FieldStop['Limits'][0] : FieldStop['Limits'][1],
                                 FieldStop['Limits'][0] : FieldStop['Limits'][1]]
                    images += crop
                else:
                    image = imr.image_reader(os.path.join(folder, im), False) - dc
                    image = image / flat_field
                    images += image

    # Accumulation Disabled
    else:

        if FieldStop['Flag']:
            new_size = FieldStop['Limits'][1] - FieldStop['Limits'][0]
            images = np.zeros((len(paths), new_size, new_size))

        else:
            images = np.zeros((len(paths), image_size, image_size))

        for idx, im in enumerate(paths):

            if FieldStop['Flag']:
                image = imr.image_reader(os.path.join(folder, im), False) - dc
                image = image / flat_field

                crop = image[FieldStop['Limits'][0] : FieldStop['Limits'][1],
                             FieldStop['Limits'][0] : FieldStop['Limits'][1]]

                images[idx] = crop
            else:
                image = imr.image_reader(os.path.join(folder, im), False) - dc
                image = image / flat_field
                images[idx] = image

    return images

def Darks_n_flats(Folder : str, image_size : int) -> list:
    """
    Function that given the dark or flat folder reads all images and computes
    the mean for each camera.

    Parameters
    ----------
    Folder : str
        Folder containing all dark images.
    image_size : int
        Size of the images in pixels.

    Returns
    -------
    matrix : list
        Matrix containg the mean for each camera. Size is equal to image.

    """
    # Identifying camera of each image
    LogFiles = LogFinder(Folder)
    C1 = []
    C2 = []
    # Loop over all logs in dark folder
    for log in LogFiles:
        file = open(log, "r")
        cam_idx = 0
        for line_idx, line in enumerate(file.readlines()):
            line = line.rstrip() # Removing end-line breaks
            # Cam ID column identification
            if line_idx == 0:
                header = line.split(',') # Separating columns
                header = [hd.replace(" ", "") for hd in header] # Values without spaces
                for hd_indx, hd in enumerate(header):
                    if hd == 'camID':
                        cam_idx = hd_indx # Storing Number of column
                        break
                    else:
                        pass
            # Images Data
            else:
                cols = line.split(',')  # Separating columns
                cols = [cl.replace(" ", "") for cl in cols] # Values without spaces
                name = os.path.basename(cols[0]) # Image name
                # Separating images depending on the camera
                if cols[cam_idx] == '0':
                    C1.append(name)
                else:
                    C2.append(name)
        file.close()

    # Reading All images of Dark folder as ark frames

    print('Folder: ' + Folder + '\nNº images found for Camera 1: ' \
           + str(len(C1)) + '\nNº images found for Camera 2: ' + str(len(C2)))

    matrix = np.zeros((2, image_size, image_size)) # Dark current storing variable

    # Looping over all dark images
    for img in C1:
        matrix[0] += imr.image_reader(os.path.join(Folder, img), False)
    for img in C2:
        matrix[1] += imr.image_reader(os.path.join(Folder, img), False)

    # Computing the mean
    matrix[0] /= len(C1)
    matrix[1] /= len(C2)

    return matrix


def DataReader(OrganizedData : dict, Obs_mode : str, image_size : int,
               dc : list, flat_field : list, Rebin : dict, FieldStop  : dict,
               Acc_flag : bool, DataFolder : str, PD_directories : dict):
    """
    Function that given an organized dictionary with the data properties, reads
    and stores the images of each of the calibration modes.

    Parameters
    ----------
    OrganizedData : dict
        Dictionary containing the useful properties of the images and already sorted.
    Obs_mode : str
        "Spectral", "Polarimetric" or "PD" for each of the calibrations.
    image_size : int
        Size of the image in pixels.
    dc : list
        Matrix containg the dark current for each camera
    flat_field : list
        Matrix containg the flat field for each camera
    Rebin : dict
        Dictionary containing boolean variable for rebining and size of the rebinning.
    FieldStop : dict
        Dictionary containing boolean variable and size for FOV selection.
    Acc_flag : bool
        Flag used to disable accumulation in PD mode
    DataFolder : str
        Path of the folder containg the images and Logs (only in Spectral and
        Polarimetric).
    PD_directories : dict
        Dictionary containing paths to focused and defocused images. (only in
        PD)
    Returns
    -------
    DATA : np.array
        Matrix containg the data of the images already sorted.

        Format:
            Spectral : (Nfilters x NModulations x NVolts x Ncamera)
            Polarimetric : (Nlogs x Nfilters x NModultaions x NVolts x Ncameras x Image_size x Image_size )
            PD : (Nmodes x NModulations x Ncameras x Nimages x Image_size x Image_size)

    additional : tuple
        Tuple containg additional info, (Voltages, Modulations, etc).

    """

    if Obs_mode == 'Spectral':

        # All filter wheels found in logs
        Wheels = [x for x in OrganizedData]
        # All modulations found in logs
        Modulations = [x for x in OrganizedData[Wheels[0]]]
        # All measured Volts
        Volts = [x for x in OrganizedData[Wheels[0]][Modulations[0]]]

        # Creating data cube to store images
        DATA = np.zeros((len(Wheels), len(Modulations), len(Volts), 2))

        # Printing data summary
        print('\nDATA PROPERTIES\nFilter Positions: ' + str(len(Wheels)) + \
              '\nModulations: ' + str(len(Modulations)) +
              '\nEtalon Volts: ' + str(len(Volts)) +
              '\nImages per Camera: ' + \
                  str(len(OrganizedData[Wheels[0]][Modulations[0]][Volts[0]]['0'])))

        print('\nOutput -> DATA( Nfilters x NModulations x NVolts x Ncamera )')

        print('\n ... Reading images ...')
        # Looping over all filters, modulations and voltages
        for wheel_idx, wheel_pos in enumerate(Wheels):
            for mod_idx, Modulation in enumerate(Modulations):
                for volts_idx, volt in enumerate(Volts):
                    for ncam in range(2):
                        Image = ImageReading(OrganizedData[wheel_pos][Modulation][volt][str(ncam)],
                                             DataFolder, image_size, dc[ncam],
                                             flat_field[ncam], Rebin, FieldStop)

                        DATA[wheel_idx, mod_idx, volts_idx, ncam] = np.mean(Image)

        additional = (Wheels, Modulations, Volts)

        return DATA, additional

    elif Obs_mode == 'PD':

        Modes = [x for x in OrganizedData]
        Wheels = [x for x in OrganizedData[Modes[0]]]
        Modulations = [x for x in OrganizedData[Modes[0]][Wheels[0]]]
        Nimages = len(OrganizedData[Modes[0]][Wheels[0]][Modulations[0]]['0'])

        if FieldStop['Flag']:
            new_size = FieldStop['Limits'][1] - FieldStop['Limits'][0]
            if Acc_flag:
                Data = np.zeros((len(Modes), len(Wheels), len(Modulations), 2, new_size,
                                 new_size))
            else:
                Data = np.zeros((len(Modes), len(Wheels), len(Modulations), 2, Nimages,
                                 new_size, new_size))
        else:
            if Acc_flag:
                Data = np.zeros((len(Modes), len(Wheels), len(Modulations), 2, image_size,
                                 image_size))
            else:
                Data = np.zeros((len(Modes), len(Wheels), len(Modulations), 2, Nimages,
                                 image_size, image_size))

        # Printing data summary
        print('\nDATA PROPERTIES\nModes: ' + str(len(Modes)) + \
              '\nFilter Positions: ' + str(len(Wheels)) + \
              '\nModulations: ' + str(len(Modulations)) +
              '\nImages per Camera: ' + str(Nimages))

        if Acc_flag:
            print('\nAccumulation ENABLED.')
            print('\nOutput -> DATA( Nmodes x N Filter x NModulations x Ncameras x Image_size x Image_size )')
        else:
            print('\nOutput -> DATA( Nmodes x N Filter x NModulations x Ncameras x Nimages x Image_size x Image_size )')

        print('\n ... Reading images ...')
        # Looping over all modulations and voltages
        for mode_idx, mode in enumerate(Modes):
            Wheels = [x for x in OrganizedData[mode]]
            for wheel_idx, wheel_pos in enumerate(Wheels):
                for mod_idx, Modulation in enumerate(Modulations):
                    for ncam in range(2):
                        Image = ImageReading(OrganizedData[mode][wheel_pos][Modulation][str(ncam)],
                                             PD_directories[mode], image_size, dc[ncam],
                                             flat_field[ncam], Rebin, FieldStop, Acc_flag)

                        Data[mode_idx, wheel_idx, mod_idx, ncam] = Image

        additional = (Modes, Wheels, Modulations)
        return Data, additional

    elif Obs_mode == 'Polarimetric':

        # All Logs wheels found in logs
        logs = [x for x in OrganizedData]

        # All filter wheels found in logs
        Wheels = [x for x in OrganizedData[logs[0]]]

        # All modulations found in logs
        Modulations = [x for x in OrganizedData[logs[0]][Wheels[0]]]

        # All measured Volts
        Volts = [x for x in OrganizedData[logs[0]][Wheels[0]][Modulations[0]]]

        # Creating data cube to store images
        if Rebin['Flag']:
            DATA = np.zeros((len(logs), len(Wheels), len(Modulations), len(Volts),
                             2, Rebin['Size'], Rebin['Size']))
        elif FieldStop['Flag']:
            new_size = FieldStop['Limits'][1] - FieldStop['Limits'][0]
            DATA = np.zeros((len(logs), len(Wheels), len(Modulations), len(Volts),
                             2, new_size, new_size))
        else:
            DATA = np.zeros((len(logs), len(Wheels), len(Modulations), len(Volts),
                             2, image_size, image_size))

        # Printing data summary
        print('\nDATA PROPERTIES\nLogs : ' + str(len(logs)) + \
              '\nFilter Positions: ' + str(len(Wheels)) + \
              '\nModulations: ' + str(len(Modulations)) +
              '\nEtalon Volts: ' + str(len(Volts)) +
              '\nImages per Camera: ' + \
                  str(len(OrganizedData[logs[0]][Wheels[0]][Modulations[0]][Volts[0]]['0'])))

        print('\nOutput -> DATA( Nlogs x Nfilters x NModultaions x NVolts x Ncameras x Image_size x Image_size )')

        print('\n ... Reading images ...')
        # Looping over all modulations and voltages
        for log_idx, log in enumerate(logs):
            for wheel_idx, wheel_pos in enumerate(Wheels):
                for mod_idx, Modulation in enumerate(Modulations):
                    for volts_idx, volt in enumerate(Volts):
                        for ncam in range(2):

                            Image = ImageReading(OrganizedData[log][wheel_pos][Modulation][volt][str(ncam)],
                                                 DataFolder, image_size, dc[ncam],
                                                 flat_field[ncam], Rebin, FieldStop)
                            DATA[log_idx, wheel_idx, mod_idx, volts_idx, ncam] = Image

        additional = (logs, Wheels, Modulations, Volts)

        return DATA, additional

    else:

        raise Exception('Please enter a valid calibration mode -> ' + \
                        '(Spectral / PD / Polarimetric) (Case sensitive)')

# =============================== MAIN ====================================== #

def CalibrationReader(Obs_mode : str, DataFolder : str, PD_directories : dict,
                      Darks : dict, Flats : dict, Rebin : dict, FieldStop : dict,
                      Acc_flag : bool, Wheel_Info : dict, Fits : dict):
    """
    Function that given mode and directories, reads the images and sorts them
    following a pre-designed criteria.

    Parameters
    ----------
    Obs_mode : str
        "Spectral", "Polarimetric" or "PD" for each of the calibrations.
    DataFolder : str
        Path of the folder containg the images and Logs (only in Spectral and
        Polarimetric).
    PD_directories : dict
        Dictionary containing paths to focused and defocused images. (only in
        PD)
    Darks : dict
        Path to the folder containing dark images observations.
    Flats : dict
        Path to the folder containing flat images observations.
    Rebin : dict
        Dictionary containing boolean variable for rebining and size of the rebinning.
    FieldStop : dict
        Dictionary containing boolean variable and size for FOV selection.
    Acc_flag : bool
        Flag used in PD mode to disable accumulation
    Wheel_Info : dict
        Info of wheel positions.
    Fits : dict
        Dictionary with fits info (Flag and Overwrite) in case it needs
        to be created'
    Returns
    -------
    Reading : np.array
        Matrix containing the data of the images .
    add : tuple
        Tuple containg additional info, (Voltages, Modulations, etc).
    organized : dict
        Dictionary used for the sorting of the data.

    """

    print('\nSelected mode: ' + Obs_mode + '\n---------------------------------')

    # Spearating observation modes
    if Obs_mode == 'Spectral' or Obs_mode == 'Polarimetric':

        print('\nDATA\nFolder : ' + DataFolder)

        # Finding all logs in folder
        LogFiles = LogFinder(DataFolder)
        print('Found ' + str(len(LogFiles)) + ' log file(s).')
        DATA = {} # Dictionary to store all properties of the images
        # Sorting Logs
        LogFiles = sorted(LogFiles)
        for logFile in LogFiles:
            LogReader(DATA, logFile, DataFolder) # Storing images names and properties

        keys = [x for x in DATA]
        image_size = int(DATA[keys[0]]['roi_x_size'])

    elif Obs_mode == 'PD':

        print('\nDATA\nFolder:\nFocused ---> '  + PD_directories['Focused'] +  \
                             '\nDefocused -> '+ PD_directories['Defocused'])

        DATA = {}
        # Finding all logs in both folders
        FocusedLogs    = LogFinder(PD_directories['Focused'])
        DefocusedLogs  = LogFinder(PD_directories['Defocused'])

        print('\nData Logs Found:\nFocused: '   + str(len(FocusedLogs)) + \
                  '\nDefocused: ' + str(len(DefocusedLogs)) )

        # Storing images names and properties
        focused_data = {}
        for logFile in FocusedLogs:
            LogReader(focused_data, logFile, PD_directories['Focused'])

        # Storing images names and properties
        defocused_data = {}
        for logFile in DefocusedLogs:
            LogReader(defocused_data, logFile, PD_directories['Defocused'])

        DATA['Focused']   = focused_data
        DATA['Defocused'] = defocused_data

        keys = [x for x in focused_data]
        image_size = int(focused_data[keys[0]]['roi_x_size'])

    else:
        raise Exception('Please enter a valid calibration mode -> ' + \
                        '(Spectral / PD / Polarimetric) (Case sensitive)')

    # Commputing Dark current and Flat fields

    if Darks['Flag']:
        print('\nDARKS')
        dc = Darks_n_flats(Darks['Folder'], image_size)
    else:
        dc = np.zeros((2, image_size, image_size))

    if Flats['Flag']:
        print('\nFLATS')
        flat_field = Darks_n_flats(Flats['Folder'], image_size) - dc
    else:
        flat_field = np.ones((2, image_size, image_size))

    # Organize the data following desired criterium for each observation mode
    organized = Organizer(DATA, Obs_mode, Wheel_Info)
    # Reading images and storing the data following the organized criteria
    Reading, add = DataReader(organized, Obs_mode, image_size, dc, flat_field, Rebin,
                              FieldStop, Acc_flag, DataFolder, PD_directories)

    if Fits['Flag']:
        hdu = fits.PrimaryHDU(Reading)


        if Obs_mode == 'Spectral':

            hdu.header['Fil_axis'] = 0
            hdu.header['Mod_axis'] = 1
            hdu.header['Vol_axis'] = 2
            hdu.header['N_filts'] = len(add[0])
            hdu.header['N_mods'] = len(add[1])
            hdu.header['N_volts'] = len(add[2])

            for filt_idx, filt in enumerate(add[0]):
                hdu.header['F_Pos_' + str(filt_idx)] = filt

            for mod_idx, mod in enumerate(add[1]):
                hdu.header['Mod_' + str(mod_idx)] = mod

            hdu.header['Min_Et_V'] = add[2][0]
            hdu.header['Max_Et_V'] = add[2][-1]
            hdu.header['V_step'] = add[2][1] - add[2][0]


        if Obs_mode == 'Polarimetric':

            hdu.header['Log_axis'] = 0
            hdu.header['Fil_axis'] = 1
            hdu.header['Mod_axis'] = 2
            hdu.header['Vol_axis'] = 3
            hdu.header['N_Logs'] = len(add[0])
            hdu.header['N_filts'] = len(add[1])
            hdu.header['N_mods'] = len(add[2])
            hdu.header['N_volts'] = len(add[3])

            hdu.header['log_0'] = add[0][0]
            hdu.header['log_last'] = add[0][-1]


            for filt_idx, filt in enumerate(add[1]):
                hdu.header['F_Pos_' + str(filt_idx)] = filt

            for mod_idx, mod in enumerate(add[2]):
                hdu.header['Mod_' + str(mod_idx)] = mod

            hdu.header['Etal_V'] = add[3][0]


        if Obs_mode == 'PD':

            hdu.header['Mode_ax'] = 0
            hdu.header['Fil_axis'] = 1
            hdu.header['Mod_axis'] = 2
            hdu.header['Cam_axis'] = 3
            hdu.header['Nimag_ax'] = 4
            hdu.header['N_modes'] = len(add[0])
            hdu.header['N_filts'] = len(add[1])
            hdu.header['N_modul'] = len(add[2])
            hdu.header['N_images'] = np.shape(Reading)[4]

            for filt_idx, filt in enumerate(add[1]):
                hdu.header['F_Pos_' + str(filt_idx)] = filt

            for mod_idx, mod in enumerate(add[2]):
                hdu.header['Mod_' + str(mod_idx)] = mod

        hdu.writeto(Fits['FileName'], overwrite = Fits['Overwrite'])

    return Reading, add, organized


def move_folders(Root,folders):
    for fold in folders:
        DataFolder = os.path.join(Root, fold)
        LogFiles = LogFinder(DataFolder)
        DATA = {}

        for log in LogFiles:
            LogReader(DATA, log, DataFolder)

        for img in DATA:
            if DATA[img]['fw2_pos'] == '0':
                Filt = '517'
            elif DATA[img]['fw2_pos'] == '1':
                Filt = '525_02'
            elif DATA[img]['fw2_pos'] == '2':
                Filt = '525_06'
            else:
                raise Exception('Problema con los filtros')

            if not os.path.isdir(os.path.join(DataFolder, Filt)):
                os.mkdir(os.path.join(DataFolder, Filt))
            os.rename(os.path.join(DataFolder, img),
                     os.path.join(DataFolder, Filt, img))
            try:

                os.rename(os.path.join(DataFolder, DATA[img]['Log']),
                         os.path.join(DataFolder, Filt, DATA[img]['Log']))
            except:
                pass
    return
"""
Current_folder = 'Lo_que_Sea'
New_folder = '/media/pablo/TOSHIBA EXT/01122021/Image E2E/Sorted_USAF_throughfocus'
Logs = LogFinder(DataFolder)
DATA = {}
for log in Logs:
    LogReader(DATA, log, DataFolder)

for img in DATA:

    if DATA[img]['fw2_pos'] == '0':
        Filt = '517'

        if DATA[img]['fw1_pos'] == '2':
            Fold = 'F4'

        elif DATA[img]['fw1_pos'] == '3':
            Fold = 'PD'
        else:
            raise Exception('Aqui hay gusano')

    elif DATA[img]['fw2_pos'] == '1':

         Filt = '525_02'
         if DATA[img]['fw1_pos'] == '2':
             Fold = 'F4'
         elif DATA[img]['fw1_pos'] == '3':
             Fold = 'PD'
         else:
             raise Exception('Aqui hay gusano')

    elif DATA[img]['fw2_pos'] == '2':

        Filt = '525_06'
        if DATA[img]['fw1_pos'] == '2':
            Fold = 'F4'
        elif DATA[img]['fw1_pos'] == '3':
            Fold = 'PD'
        else:
            raise Exception('Aqui hay gusano')

    else:
        raise Exception('Aqui hay gusano')

    copyfile(os.path.join(DataFolder, img), os.path.join(New_folder, Filt + '_' + Fold, img))
    copyfile(os.path.join(DataFolder, DATA[img]['Log']), os.path.join(New_folder,  Filt + '_' + Fold, DATA[img]['Log']))
"""
