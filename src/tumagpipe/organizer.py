import csv
import numpy as np

import glob
import os
from datetime import datetime
import bisect

data_main_folder = "/work/obs/TuMAG_data"

day_folders = sorted(os.listdir(data_main_folder))
day_folders = [x for x in day_folders if os.path.isdir(f"{data_main_folder}/{x}")]

D09_paths = []
D10_paths = []
D11_paths = []
D12_paths = []
D13_paths = []
D14_paths = []
D15_paths = []
D16_paths = []

D09_imagenames = []
D10_imagenames = []
D11_imagenames = []
D12_imagenames = []
D13_imagenames = []
D14_imagenames = []
D15_imagenames = []
D16_imagenames = []

for primary_folder in day_folders:
    print(f"Processing folder: {primary_folder}")
    hour_folders = os.listdir(f"{data_main_folder}/{primary_folder}")

    for secondary_folder in hour_folders:
        all_files = sorted(glob.glob(f"{data_main_folder}/{primary_folder}/{secondary_folder}/*img"))

        for image in all_files:

            image_name = os.path.basename(image)

            day = image_name[8:10]

            if day == "09":
                D09_paths.append(image)
                D09_imagenames.append(image_name)
            elif day == "10":
                D10_paths.append(image)
                D10_imagenames.append(image_name)

            elif day == "11":
                D11_paths.append(image)
                D11_imagenames.append(image_name)

            elif day == "12":
                D12_paths.append(image)
                D12_imagenames.append(image_name)

            elif day == "13":
                D13_paths.append(image)
                D13_imagenames.append(image_name)

            elif day == "14":
                D14_paths.append(image)
                D14_imagenames.append(image_name)

            elif day == "15":
                D15_paths.append(image)
                D15_imagenames.append(image_name)

            elif day == "16":
                D16_paths.append(image)
                D16_imagenames.append(image_name)

            else:
                print("eeeeeeeh raro")
                raise Exception("quieto parao")


sort_indices = np.argsort(D09_imagenames)
D09_paths = np.array(D09_paths)
sorted_paths = D09_paths[sort_indices]
csvfile = open("Organized_files_v2/D09.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D09_paths)}")
csvfile.close()


sort_indices = np.argsort(D10_imagenames)
D10_paths = np.array(D10_paths)
sorted_paths = D10_paths[sort_indices]
csvfile = open("Organized_files_v2/D10.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D10_paths)}")
csvfile.close()

sort_indices = np.argsort(D11_imagenames)
D11_paths = np.array(D11_paths)
sorted_paths = D11_paths[sort_indices]
csvfile = open("Organized_files_v2/D11.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D11_paths)}")
csvfile.close()

sort_indices = np.argsort(D12_imagenames)
D12_paths = np.array(D12_paths)
sorted_paths = D12_paths[sort_indices]
csvfile = open("Organized_files_v2/D12.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D12_paths)}")
csvfile.close()

sort_indices = np.argsort(D13_imagenames)
D13_paths = np.array(D13_paths)
sorted_paths = D13_paths[sort_indices]
csvfile = open("Organized_files_v2/D13.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D13_paths)}")
csvfile.close()

sort_indices = np.argsort(D14_imagenames)
D14_paths = np.array(D14_paths)
sorted_paths = D14_paths[sort_indices]
csvfile = open("Organized_files_v2/D14.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D14_paths)}")
csvfile.close()

sort_indices = np.argsort(D15_imagenames)
D15_paths = np.array(D15_paths)
sorted_paths = D15_paths[sort_indices]
csvfile = open("Organized_files_v2/D15.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D15_paths)}")
csvfile.close()

sort_indices = np.argsort(D16_imagenames)
D16_paths = np.array(D16_paths)
sorted_paths = D16_paths[sort_indices]
csvfile = open("Organized_files_v2/D16.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)
for ind, file in enumerate(sorted_paths):
    csv_writer.writerow([ind, file])
    print(f"{ind}/{len(D16_paths)}")
csvfile.close()








"""



csvfile = open("Organized_files/D09.csv", 'w', newline='')
csv_writer = csv.writer(csvfile)

day = "09"
counter = 0
for folder in day_folders:
    print(f"Organizing folder: {folder}")
    
    timestamp = get_time(folder[:19])
    folder_day = folder[8:10]
    print(f"Detected day: {folder_day}")

    if folder_day != day:
        csvfile.close()
        day = folder_day
        print(f"Creating file: 'Organized_files/D{day}.csv'") 
        csvfile = open(f"Organized_files/D{day}.csv", 'w', newline='')
        csv_writer = csv.writer(csvfile)
        counter = 0
    
    hour_folders = sorted(os.listdir(f"{data_main_folder}/{folder}"))
    hour_folders = [x for x in hour_folders if os.path.isdir(f"{data_main_folder}/{folder}/{x}")]

    for hour_fold in hour_folders:

        all_files = sorted(glob.glob(f"{data_main_folder}/{folder}/{hour_fold}/*img"))
        for file in all_files:
            csv_writer.writerow([counter, file])
            counter += 1
        print(counter)

"""