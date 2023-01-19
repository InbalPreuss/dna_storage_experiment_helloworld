import csv
from typing import List, Dict
import re
import os

import pandas as pd


def sorted_human(iterable: List[str]) -> List[str]:
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(iterable, key=alphanum_key)


def is_dir_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def write_dict_to_csv(dict, csv_path):
    with open(csv_path, 'w') as f:  # You will need 'wb' mode in Python 2.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)


def open_csv_file_to_df(file_name):
    df = pd.read_csv(file_name)
    return df


def hamming_distance(string1, string2):
    # Start with a distance of zero, and count up
    distance = 0
    # Loop over the indices of the string
    L = len(string1)
    for i in range(L):
        # Add 1 to the distance if these two characters are not equal
        if string1[i] != string2[i]:
            distance += 1
    # Return the final count of differences
    return distance


def dict_to_csv(dict_dist: Dict, file_name: str):
    df = pd.DataFrame(dict_dist)
    df.to_csv(file_name)
