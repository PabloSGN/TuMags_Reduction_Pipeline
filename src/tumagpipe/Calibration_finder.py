# IMPORTS
import pandas as pd
from astropy.io import fits
import image_handler as ih

# CONFIG 
dark_csv = "Darks_fits_finder.csv"  # CSV file containing dark current data
flat_csv = "Flats_fits_finder.csv"  # CSV file containing flat field data
pd.options.mode.chained_assignment = None  # Disable chained assignment warning for cleaner terminal output

# MAIN CODE
def find_closest_date(df, target_date):
    """
    Finds the row with the closest date to the target date.

    Inputs: 
    df : pandas DataFrame
        dark or flat info DataFrames
    target_date : datetime
        The target date to find the closest match.

    Outputs:
    pandas Series
        The row in df with the closest date to target_date.
    """

    df["Time Difference"] = abs(pd.to_datetime(df["Date"]) - target_date)
    closest_row = df.loc[df["Time Difference"].idxmin()] # Find closest date
    df.drop(columns=["Time Difference"], inplace=True) # Add time difference for more complete info
    return closest_row

def display_row(row):
    """
    Displays the contents in terminal of a row in a readable format.
    Inputs: 
    row : pandas Series
        A row from a DataFrame.
    """
    print("Selected file:\n")
    for col, value in row.items():
        print(f"    - {col}: {value}")

def get_closest_dark_current(target_date, verbose=True):
    """
    Gets the closest dark current to the target date.

    Inputs:
    target_date : datetime ("%m/%d/%Y, %H:%M:%S" format)
        The target date to search for.
    verbose : bool, optional (default=True)
        If True, displays details of the selected file.

    Outputs:
    numpy.ndarray
        FITS data from the closest matching file.
    """
    df = pd.read_csv(dark_csv)
    closest_file = find_closest_date(df, target_date)

    if verbose:
        print(f"\nClosest dark to {target_date}:")
        display_row(closest_file)

    return fits.getdata(closest_file["Filepath"])

def get_dark_index(target_indexes, verbose=True):
    """
    Gets the dark current by its indexes.

    Inputs:
    target_indexes : str  (DXX-nnnn-nnnn)
        The target indexes.
    verbose : bool, optional (default=True)
        If True, displays details of the selected file.
    Outputs:
    numpy.ndarray
        FITS data from the matching file.
    """
    dummy_filename = f"dc_{target_indexes}.fits" # Creating a dummy strg to search for
    df = pd.read_csv(dark_csv)
    matching_row = df.loc[df["Filename"] == dummy_filename]

    if matching_row.empty:
        raise Exception(f"Dark with indexes {target_indexes} not found.")
    else:
        if verbose:
            matching_row = matching_row.iloc[0]
            print(f"\nFound dark with indexes: {target_indexes}:")
            display_row(matching_row)

        return fits.getdata(matching_row["Filepath"])

def get_flat_index(target_indexes, verbose=True):
    """
    Gets a flat field by its indexes.

    Inputs:
    target_indexes : str (DXX-nnnn-nnnn) / (DYY-nnnn-nnnn)
        The target indexes used to identify the file.
    verbose : bool, optional (default=True)
        If True, displays details of the selected file.

    Outputs:
    tuple (numpy.ndarray, astropy.io.fits.Header)
        FITS data and header from the matching file.
    """
    df = pd.read_csv(flat_csv)
    matching_row = df[df["Indexes"] == target_indexes]

    if matching_row.empty:
        raise Exception(f"\nFlat with indexes {target_indexes} not found.")
    else:
        if verbose:
            matching_row = matching_row.iloc[0]
            print(f"\nFound flat with indexes: {target_indexes}:")
            display_row(matching_row)

        return fits.getdata(matching_row["Filepath"]), fits.getheader(matching_row["Filepath"])

def get_flat_date(target_date, obs_mode, verbose=True):
    """
    Gets the closest flat field to a target date and observation mode.

    Inputs:
    target_date : datetime ("%m/%d/%Y, %H:%M:%S" format)
        The target date to search for.
    obs_mode : str
        The observation mode to filter by.
    verbose : bool, optional (default=True)
        If True, displays details of the selected file.

    Outputs:
    tuple (numpy.ndarray, astropy.io.fits.Header)
        FITS data and header from the closest matching file.
    """
    df = pd.read_csv(flat_csv)
    filtered_df = df[df["Obs_modes"] == obs_mode]

    matching_row = find_closest_date(filtered_df, target_date)

    if verbose:
        print(f"\nFound closest flat to: {target_date} for obs mode:{obs_mode}:")
        display_row(matching_row)         

    return fits.getdata(matching_row["Filepath"]), fits.getheader(matching_row["Filepath"])


def find_calibration(image_index, dark_index = None, flat_index = None, verbose = True):
    """
    Finds the calibration files for a given image index.

    Inputs:
    image_index : str
        The index of the image to find calibration for.
    dark_index : str, optional (default=None)
        The index of the dark current to use. If None, the closest dark current is used.
    flat_index : str, optional (default=None)
        The index of the flat field to use. If None, the closest flat field is used.

    Outputs:
    tuple (numpy.ndarray, numpy.ndarray, astropy.io.fits.Header)
        Dark current data, flat field data, flat header.
    """
    # Get image data
    _, H = ih.read_ID(image_index)

    date = H["Date"]
    obs_mode = H["ObservationMode"]

    if verbose:
        print(f"Finding calibration for Obsmode: {obs_mode} on {date}.")
       

    # Get dark current
    if dark_index is not None:
        print(f"Searching darks with index {dark_index}")
        dark_data = get_dark_index(dark_index, verbose=verbose)
    else:
        print(f"Automatic dark search")
        dark_data = get_closest_dark_current(date, verbose=verbose)

    # Get flat field
    if flat_index is not None:
        print("Searching flats with index {flat_index}")
        flat_data, flat_header = get_flat_index(flat_index, verbose=verbose)
    else:
        print("Automatic flat search")
        flat_data, flat_header = get_flat_date(date, obs_mode, verbose=verbose)

    return dark_data, flat_data, flat_header