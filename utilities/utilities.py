import csv
from pathlib import Path
from typing import List, Dict, Union
import re
import os
import subprocess
import pandas as pd
import time
from functools import wraps
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.monotonic()
        result = func(*args, **kwargs)
        end_time = time.monotonic()
        print(f"{func.__name__},{(end_time - start_time):.2f},seconds to run.")

        return result

    return wrapper


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


def write_dict_to_csv_as_dict(dict: Dict, csv_path: Path):
    with open(csv_path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in dict.items():
            writer.writerow([key, value])


def dict_to_csv(dict: Dict, file_name: str):
    df = pd.DataFrame(dict)
    df.to_csv(file_name)


def write_list_to_csv(data: List, file_name: Union[Path, str]) -> None:
    with open(file_name, 'a', encoding='UTF8') as f:
        writer = csv.writer(f)
        writer.writerow(data)


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


def open_fasta(input_file: Path) -> List[str]:
    with open(input_file, 'r') as inf:
        reads = list(SeqIO.parse(inf, 'fasta'))

    return reads


def open_fastq(input_file: Path) -> List[str]:
    with open(input_file, 'r') as inf:
        reads = list(SeqIO.parse(inf, 'fastq'))

    return reads


def open_fastq_yield(input_file: Union[str, Path]):
    with open(str(input_file), 'r') as inf:
        for record in SeqIO.parse(inf, 'fastq'):
            yield record


def get_amount_of_reads_from_file(file_path: Union[Path, str]):
    # Convert the string to a Path object
    file_path = Path(file_path)
    # Determine the file format based on the file extension
    if file_path.suffix == '.fastq':
        file_format = 'fastq'
    elif file_path.suffix == '.fasta':
        file_format = 'fasta'
    else:
        raise ValueError("Unsupported file format. The file must be either .fastq or .fasta")

    # Open the file and create an index
    index = SeqIO.index(str(file_path), file_format)

    # Get the number of sequences in the index
    num_sequences = len(index)
    return num_sequences


def run_command_in_git_bash(command):
    git_bash_path = "C:\\Program Files\\Git\\bin\\bash.exe"
    return subprocess.run([git_bash_path, "-c", command], capture_output=True)


def fasta_to_fastq(input_file_or_folder, default_quality=40):
    # Ensure the input is a directory
    if not Path(input_file_or_folder).is_dir():
        raise ValueError("The input must be a directory.")

    # Ensure the input is a directory
    if not Path(input_file_or_folder).is_dir():
        raise ValueError("The input must be a directory.")

    # Iterate over files in the directory
    for file_path in Path(input_file_or_folder).iterdir():
        if file_path.suffix == '.fasta':
            # Define the output FASTQ file name
            fastq_file = file_path.with_suffix('.fastq')

            # Read FASTA file
            fasta_sequences = SeqIO.parse(str(file_path), "fasta")

            # Create FASTQ records with default quality scores
            fastq_records = []
            for fasta_record in fasta_sequences:
                # Create a string of default quality scores
                quality_scores = [default_quality] * len(fasta_record.seq)

                # Create a SeqRecord object for FASTQ
                fastq_record = SeqRecord(Seq(fasta_record.seq),
                                         id=fasta_record.id,
                                         description=fasta_record.description,
                                         letter_annotations={"phred_quality": quality_scores})

                fastq_records.append(fastq_record)

            # Write FASTQ file
            with open(fastq_file, "w") as output_handle:
                SeqIO.write(fastq_records, output_handle, "fastq")