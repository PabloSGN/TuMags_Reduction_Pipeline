"""
.. py:module:: TuMAGtools.utils
.. module:: utils
        :platform: Unix
        :synopsis: function for reading TuMAG images and headers
.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""

# ============================ IMPORTS ====================================== #

import numpy as np
import struct
import matplotlib.pyplot as plt
from astropy.io import fits

# ============================= CONFIG ====================================== #

# ------------------------------ SIZE IN BYTES ------------------------------ #

# Size of bytes of the different fields in the header of the image
nBytesMpsHeader = 20
nBytesLengthImageStruct = 2
nBytesLengthSSt = 2
nBytesLengthTail = 8  # SSC + CRC

# COMMON sizes
# cameraId                uint8_t     --> struct.unpack('B')
nBytesCameraId = 1
# timeStamp_start         Uint64_t    --> struct.unpack('Q')
nBytesTimeStamp_start = 8
# timeStamp_end           Uint64_t    --> struct.unpack('Q')
nBytesTimeStamp_end = 8
# imageSize               uint32_t    --> struct.unpack('L')
nBytesImageSize = 4
# observationMode         uint8_t     --> struct.unpack('B')
nBytesObservationMode = 1
# componentID             uint8_t     --> struct.unpack('B')
nBytesComponentID = 1
# pipelineConfig          uint8_t     --> struct.unpack('B')
nBytesPipelineConfig = 1
# NAcc                    uint16_t    --> struct.unpack('H')
nBytesnAcc = 2
# imageIndex              uint32_t    --> struct.unpack('L')
nBytesImageIndex = 4
# imageIndex_end          uint32_t    --> struct.unpack('L')
nBytesImageIndex_end = 4
# roi_x_offset            uint16_t    --> struct.unpack('H')
nBytesRoi_x_offset = 2
# roi_x_size              uint16_t    --> struct.unpack('H')
nBytesRoi_x_size = 2
# roi_y_offset            uint16_t    --> struct.unpack('H')
nBytesRoi_y_offset = 2
# roi_y_size              uint16_t    --> struct.unpack('H')
nBytesRoi_y_size = 2
# Image type              uint8_t     --> struct.unpack('B')
nBytesImageType = 1
# OM_counter              uint16_t    --> struct.unpack('H')
nBytesOMCounter = 2

# TuMAG sizes
# position FW1            uint8_t     --> struct.unpack('B')
nBytesFW1 = 1
# position FW2            uint8_t     --> struct.unpack('B')
nBytesFW2 = 1
# voltage etalon DN       uint16_t    --> struct.unpack('H')
nBytesEtalonDN = 2
# sign etalon             uint8_t     --> struct.unpack('B')
nBytesEtalonSign = 1
# Rocli1_LCVR             uint16_t    --> struct.unpack('H')
nBytesRocli1_LCVR = 2
# Rocli2_LCVR             uint16_t    --> struct.unpack('H')
nBytesRocli2_LCVR = 2
# Etalon Volts reading    uint16_t    --> struct.unpack('H')
nBytesEtalonVoltLecture = 2
# Real FW1 Pos            uint8_t     --> struct.unpack('B')
nBytesfw1PosReal = 1
# Real FW2 Pos            uint8_t     --> struct.unpack('B')
nBytesfw2PosReal = 1
# Counts lcvr1            uint16_t    --> struct.unpack('H')
nByteslcvr1DNReal = 2
# Counts lcvr2            uint16_t    --> struct.unpack('H')
nByteslcvr2DNReal = 2

# ------------------------------- THUMBNAILS -------------------------------- #

# Image types in HEADER VERSION 5
Thumb_type_5 = {
    0: {'type': 'Image',                            'sf': 1},
    10: {'type': 'No Data',                         'sf': 1},
    11: {'type': 'Th. Binning Los',                 'sf': 4},
    12: {'type': 'Th. Binning No Los',              'sf': 8},
    13: {'type': 'Th. Full Test',                   'sf': 1},
    14: {'type': 'Th. Cropped Centered Los',        'sf': 4},
    15: {'type': 'Th. Cropped Centered No Los',     'sf': 8},
    16: {'type': 'Th. Cropped Up left Los',         'sf': 4},
    17: {'type': 'Th. Cropped Up left No Los',      'sf': 8},
    18: {'type': 'Th. Cropped Up Center Los',       'sf': 4},
    19: {'type': 'Th. Cropped Up Center No Los',    'sf': 8},
    20: {'type': 'Th. Cropped Up Right Los',        'sf': 4},
    21: {'type': 'Th. Cropped Up Right No Los',     'sf': 8},
    22: {'type': 'Th. Cropped Center left Los',     'sf': 4},
    23: {'type': 'Th. Cropped Center left no Los',  'sf': 8},
    24: {'type': 'Th. Cropped Center Right Los',    'sf': 4},
    25: {'type': 'Th. Cropped Center Right no Los', 'sf': 8},
    26: {'type': 'Th. Cropped Down left Los',       'sf': 4},
    27: {'type': 'Th. Cropped Down left no Los',    'sf': 8},
    28: {'type': 'Th. Cropped Down center Los',     'sf': 4},
    29: {'type': 'Th. Cropped Down center no Los',  'sf': 8},
    30: {'type': 'Th. Cropped Down right Los',      'sf': 4},
    31: {'type': 'Th. Cropped Down right no Los',   'sf': 8}
}


# ------------------------------- THUMBNAILS -------------------------------- #

# Image types in HEADER VERSION 6
Thumb_type_6 = {
    0:  {'type': 'Image',          'bin': 1},
    10: {'type': 'No Data',        'bin': 1},
    13: {'type': 'Full Image',     'bin': 1},
    32: {'type': 'Th. Binning 2',  'bin': 2},
    33: {'type': 'Th. Binning 4',  'bin': 4},
    34: {'type': 'Th. Binning 8',  'bin': 8},
    35: {'type': 'Th. Binning 16', 'bin': 16},
    36: {'type': 'Th. Binning 32', 'bin': 32},
    37: {'type': 'Th. Cropped 2',  'bin': 2},
    38: {'type': 'Th. Cropped 4',  'bin': 4},
    39: {'type': 'Th. Cropped 8',  'bin': 8},
    40: {'type': 'Th. Cropped 16', 'bin': 16},
    41: {'type': 'Th. Cropped 32', 'bin': 32},
}

# =========================================================================== #

def GetDatafromHeader(receivedHeader):
    """
    Function that reads the bytes in the header and extracts the differetn fields, 
    storing them in a dictionary.
    Parameters
    ----------
    receivedHeader : bytes
        Bytes object containing the header of the image.

    Returns
    -------
    Header : Dict
        Dictionary containing the different fields included in the header.
    """

    dataLine = []

    # Keys of the variables that get stored in the header dictionary
    Keys = ['CameraID', 'TimeStamp_start', 'TimeStamp_end', 'ImageSize',
            'ObservationMode', 'PipelineConfig', 'nAcc', 'ImageIndex', 'ImageIndex_end',
            'Roi_x_offset', 'Roi_x_size', 'Roi_y_offset', 'Roi_y_size', 'ImageType', 'Observation_Counter',
            'FW1', 'FW2', 'EtalonDN', 'EtalonSign', 'Rocli1_LCVR', 'Rocli2_LCVR', 'EtalonVoltsReading',
            'FW1_Real', 'FW2_Real', 'LCVR1_DN_Real', 'LCVR2_DN_Real']

    # Reading size of Headers from file
    positionStartSstLength = nBytesMpsHeader
    nBytesSSt = struct.unpack('H', receivedHeader[positionStartSstLength: positionStartSstLength + nBytesLengthSSt])[0]
    positionStartImageStructLength = positionStartSstLength + nBytesLengthSSt + nBytesSSt

    # COMMON positions in bytes
    positionStartCameraId        = positionStartImageStructLength + nBytesLengthImageStruct
    positionStartTimeStamp_start = positionStartCameraId + nBytesCameraId
    positionStartTimeStamp_end   = positionStartTimeStamp_start + nBytesTimeStamp_start
    positionStartImageSize       = positionStartTimeStamp_end + nBytesTimeStamp_end
    positionStartObservationMode = positionStartImageSize + nBytesImageSize
    positionStartPipelineConfig  = positionStartObservationMode + nBytesObservationMode
    positionStartnAcc            = positionStartPipelineConfig + nBytesPipelineConfig
    positionStartImageIndex      = positionStartnAcc + nBytesnAcc
    positionStartImageIndex_end  = positionStartImageIndex + nBytesImageIndex
    positionStartRoi_x_offset    = positionStartImageIndex_end + nBytesImageIndex_end
    positionStartRoi_x_size      = positionStartRoi_x_offset + nBytesRoi_x_offset
    positionStartRoi_y_offset    = positionStartRoi_x_size + nBytesRoi_x_size
    positionStartRoi_y_size      = positionStartRoi_y_offset + nBytesRoi_y_offset
    positionStartImageType       = positionStartRoi_y_size + nBytesRoi_y_size
    positionStartOmCounter       = positionStartImageType + nBytesImageType

    # TuMAG positions in bytes
    positionFW1               = positionStartOmCounter + nBytesOMCounter
    positionFW2               = positionFW1 + nBytesFW1
    positionEtalonDN          = positionFW2 + nBytesFW2
    positionEtalonSign        = positionEtalonDN + nBytesEtalonDN
    positionRocli1_LCVR       = positionEtalonSign + nBytesEtalonSign
    positionRocli2_LCVR       = positionRocli1_LCVR + nBytesRocli1_LCVR
    positionEtalonVoltReading = positionRocli2_LCVR + nBytesRocli2_LCVR
    positionFW1Real           = positionEtalonVoltReading + nBytesEtalonVoltLecture
    positionFW2Real           = positionFW1Real + nBytesfw1PosReal
    positionLCVR1_DN_real     = positionFW2Real + nBytesfw2PosReal
    positionLCVR2_DN_real     = positionLCVR1_DN_real + nByteslcvr1DNReal

    # Extract COMMON Data
    dataLine.append(struct.unpack('B', receivedHeader[positionStartCameraId: positionStartCameraId + nBytesCameraId])[0])
    dataLine.append(struct.unpack('Q', receivedHeader[positionStartTimeStamp_start: positionStartTimeStamp_start + nBytesTimeStamp_start])[0])
    dataLine.append(struct.unpack('Q', receivedHeader[positionStartTimeStamp_end: positionStartTimeStamp_end + nBytesTimeStamp_end])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageSize: positionStartImageSize + nBytesImageSize])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartObservationMode: positionStartObservationMode + nBytesObservationMode])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartPipelineConfig: positionStartPipelineConfig + nBytesPipelineConfig])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartnAcc: positionStartnAcc + nBytesnAcc])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageIndex: positionStartImageIndex + nBytesImageIndex])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageIndex_end: positionStartImageIndex_end + nBytesImageIndex_end])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_x_offset: positionStartRoi_x_offset + nBytesRoi_x_offset])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_x_size: positionStartRoi_x_size + nBytesRoi_x_size])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_y_offset: positionStartRoi_y_offset + nBytesRoi_y_offset])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_y_size: positionStartRoi_y_size + nBytesRoi_y_size])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartImageType: positionStartImageType + nBytesImageType])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartOmCounter: positionStartOmCounter + nBytesOMCounter])[0])

    # Extract TuMAG Data
    dataLine.append(struct.unpack('B', receivedHeader[positionFW1: positionFW1 + nBytesFW1])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionFW2: positionFW2 + nBytesFW2])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionEtalonDN: positionEtalonDN + nBytesEtalonDN])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionEtalonSign: positionEtalonSign + nBytesEtalonSign])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionRocli1_LCVR: positionRocli1_LCVR + nBytesRocli1_LCVR])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionRocli2_LCVR: positionRocli2_LCVR + nBytesRocli2_LCVR])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionEtalonVoltReading: positionEtalonVoltReading + nBytesEtalonVoltLecture])[0])

    # Temporal Parameters
    dataLine.append(struct.unpack('B', receivedHeader[positionFW1Real: positionFW1Real + nBytesfw1PosReal])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionFW2Real: positionFW2Real + nBytesfw2PosReal])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionLCVR1_DN_real: positionLCVR1_DN_real + nByteslcvr1DNReal])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionLCVR2_DN_real: positionLCVR2_DN_real + nByteslcvr2DNReal])[0])

    # Saving data in dictionary
    Header = {}
    for ind, key in enumerate(Keys):
        Header[key] = dataLine[ind]

    return Header

def read_Tumag(file, write_fits = False, fits_file_name = 'Image.fits', 
               plot_flag = False, vmin = 0, vmax = 4096, onlyheader = False):
   
    """
    Function that reads TuMag images and thumbnails and return both the data 
    and the header. 

    Parameters
    ----------
    file : str
        String containing the path to an image.
    write_fits : Boolean, optional
       Boolean variable that selects the opction of writing a fits. 
    fits_file_name : str, optional
       String containing the name of the fits file if write_fits = True.
       The default is 'Image.fits'.
    plot_flag : Boolean, optional
        Boolean variable that selects the option of plotting the image. 
        The default is False.
    vmin : int, optional
        Minimum value of the map if plt_flag = True. The default is 0.
    vmax : int, optional
        Minimum value of the map if plt_flag = True. The default is 4096.
    onlyheader : Boolean, optional
        Boolean variable that selects the option of returning only heade.
        Default is False.
        
    Returns
    -------
    H : dict
        Dictionary containing the info read in the header of the image.
    Image : np.array
        Array containing the data of the image. Onlye returned if image contains
        data.

    """

    # Process the MPS Science Header (nBytesMpsHeader bytes --> 20 bytes)
    nBytesSync = 2              # Identify start of record --> 0x44 0x54
    nBytesSystemId = 1          # System ID of generated data
    nBytesDataId = 1            # Description of data type
    # Total data packet length (including header)
    nBytesTotalLength = 4
    nBytesTime = 8              # Data acquisition start time
    nBytesSensorId = 1          #
    nBytesHeaderVersion = 1     # Header Version
    # Total length of header (in bytes). Take it as an offset to Data
    nBytesHeaderLength = 2

    # Compute start position in bytes of total length, Header version and header length
    positionStartTotalLength = nBytesSync + nBytesSystemId + nBytesDataId
    positionStartHeaderVersion = nBytesSync + nBytesSystemId + \
        nBytesDataId + nBytesTotalLength + nBytesTime + nBytesSensorId
    positionStartHeaderLength = nBytesSync + nBytesSystemId + nBytesDataId + \
        nBytesTotalLength + nBytesTime + nBytesSensorId + nBytesHeaderVersion

    # Start reading the file
    file = open(file, "rb")
    fullReceivedImage = file.read()

    # Extract mps header and needed fields
    mpsHeader = fullReceivedImage[0: nBytesMpsHeader]
    totalLength   = struct.unpack('i', mpsHeader[positionStartTotalLength: positionStartTotalLength + nBytesTotalLength])[0]
    headerVersion = struct.unpack('B', mpsHeader[positionStartHeaderVersion: positionStartHeaderVersion + nBytesHeaderVersion])[0]
    headerLength  = struct.unpack('H', mpsHeader[positionStartHeaderLength: positionStartHeaderLength + nBytesHeaderLength])[0]

    # Separate Header from image
    receivedHeader = fullReceivedImage[0:headerLength]
    receivedImage  = fullReceivedImage[headerLength: totalLength - nBytesLengthTail]

    # Extract TuMag Header info
    H = GetDatafromHeader(receivedHeader)

    if onlyheader:
        return H 
    
    else:

        # Default values
        bytesppx = 2
        dtypef = '<i2'
        # Bytes per pixel (2 or 4)
        if H['PipelineConfig'] == 0:
            bytesppx = 2
            dtypef = '<i2'
        if H['PipelineConfig'] == 1:
            bytesppx = 4
            dtypef = '<i4'
        if H['PipelineConfig'] == 2:
            bytesppx = 4
            dtypef = '<i4'
    
        ImageFlag = False
    
        # Cheking the size of the image
    
        # PREVIOUS VERSION
        if headerVersion == 5:
    
            # Compute Size of image
            # Normal Image
            if H['ImageType'] == 0:
                ImageFlag = True
                width, height = H['Roi_x_size'], H['Roi_y_size']
            elif H['ImageType'] == 10:
                print('Thumbnail containing only header')
                return H 
            # Thumbnails
            else:
                ImageFlag = True
                factor_bytesppx = bytesppx // 2
                Size_factor = Thumb_type_5[H['ImageType']]['sf'] # Get the binning
    
                # Real width and height
                width = H['Roi_x_size'] // (Size_factor * factor_bytesppx)
                height = H['Roi_y_size'] // (Size_factor * factor_bytesppx)
            
        # CURRENT VERSION -- Header Version 6
        else:
            # Compute Size of image
    
            # Normal Image
            if H['ImageType'] == 0:
                ImageFlag = True
                width, height = H['Roi_x_size'], H['Roi_y_size']
                H["Thumbnail"] = False
            
            elif H['ImageType'] == 10:
                
                return H 
    
            # Thumbnails
            else:
                
                ImageFlag = True
                H["Thumbnail"] = True

                Binning = Thumb_type_6[H['ImageType']]['bin'] # Get the binning
                H["Thumbnail_binning"] = Binning

                # Real width and height
                width = H['Roi_x_size'] // Binning
                height = H['Roi_y_size'] // Binning
    
        # If the file contained an image (not header only thumbnail type)
        if ImageFlag:
            
            # Read the image from bytes format
            Image = np.frombuffer(receivedImage, dtype=dtypef).reshape([height, width]).astype(np.uint16)
            
            # If a fits file is wanted to be written
            if write_fits:
                
                FITS = fits.PrimaryHDU(Image)
                head = FITS.header
                
                # Cattegories compatible with fits file formatting
                ReducedKeys = ['CameraID', 'T_start', 'T_end', 'Img_size', 'OM', 
                               'Pipeconf', 'nacc', 'Img_idx', 'Img_idxe', 'Roix_off',
                               'Roi_X', 'Roiy_off', 'Roi_y', 'ImgType', 'OM_Count', 
                               'FW1', 'FW2',  'EtalonDN', 'Etal_sig','LCVR1',
                               'LCVR2', 'Et_real', 'FW1_real', 'FW2_real', 
                               'LCVR1_re', 'LCVR2_re']
                
                # Description of the fields
                Comments = ['Camera number', 'Timestamp start', 'Timestamp end', 'Real size (bytes)',
                            'Observation Mode', 'Pipeline Config', 'Accumulations number', 
                            'Index of image', 'End of image index', 'ROI X Offset', 
                            'ROI X', 'ROI Y OFFSET', 'ROY', 'Image type', 'Observation Mode Counter',
                            'Filter wheel 1 pos', 'Filter wheel 2 pos', 'Etalon Volts (Counts)', 
                            'Sign of etalon volts', 'LCVR1 volts (Counts)', 'LCVR2 volts (Counts)', 
                            'Measured value for etalon volts (counts)', 'Measured pos of Filter wheel 1', 
                            'Measured pos of Filter Wheel 2', 'Measured volts for LCVR1 (Counts)',
                            'Measured volts for LCVR2 (Counts)']
                
                # Saving the header
                for ind, key in enumerate(H):
                    head.append((ReducedKeys[ind], H[key], Comments[ind]))
                    
                # Fits writing
                FITS.writeto(fits_file_name, overwrite=False)
                    
            # If a plot is wanted to be shown
            if plot_flag:
                plt.imshow(Image, cmap = 'inferno', vmin = vmin, vmax = vmax)
                plt.colorbar()
                plt.show()
        
            return Image, H

