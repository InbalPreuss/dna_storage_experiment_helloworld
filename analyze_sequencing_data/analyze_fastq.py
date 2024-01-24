import ast
import collections
import heapq
import os
import shutil

from pathlib import Path
from typing import Union, Dict, List, Tuple, Any

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from Levenshtein import distance as lev
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt, patches
from pandas import DataFrame
from pandas.io.parsers import TextFileReader
from tqdm import tqdm
import csv
import subprocess
import itertools
import re

import utilities.utilities as uts


def compare_with_errors(read, seq, max_dist=6):
    return sum(r != s for r, s in zip(read, seq)) <= max_dist


class AnalyzeFastqData:
    def __init__(self, input_file_or_folder: Union[Path, str],
                 const_design_file: Union[Path, str],
                 payload_design_file: Union[Path, str],
                 barcodes_design_file: Union[Path, str],
                 len_reads_hist_output_file: Union[Path, str],
                 results_good_reads_file: Union[Path, str],
                 results_good_reads_with_len_bigger_then_y: Union[Path, str],
                 results_most_common_file: Union[Path, str],
                 design_simulation_file: Union[Path, str],
                 design_results_only_z_file: Union[Path, str],
                 compare_design_to_experiment_results_output_file: Union[Path, str],
                 foreach_bc_payload_count_file: Union[Path, str],
                 heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str],
                 count_reads_for_each_bc_file: Union[Path, str],
                 missing_bcs_file: Union[Path, str],
                 output_hist_folder: Union[Path, str],
                 output_folder: Union[Path, str],
                 design_folder: Union[Path, str],
                 output_graphs_folder: Union[Path, str],
                 output_csv_folder: Union[Path, str],
                 output_heatmap_folder: Union[Path, str],
                 hist_per_bc_file: Union[Path, str],
                 hist_foreach_bc_read_count_file: Union[Path, str],
                 hist_foreach_error_count_of_bc_file: Union[Path, str],
                 hist_foreach_read_count_count_bc_file: Union[Path, str],
                 len_reads_to_retrieve: int,
                 amount_of_bc: int,
                 barcode_len: int,
                 subset_size: int,
                 amount_of_payloads: int,
                 z_to_k_mer_representative: Dict,
                 k_mer_representative_to_z: Dict,
                 payload_pos: List,
                 uni_start_pos: List,
                 sampling_rate_from_good_reads_graph: Union[Path, str],
                 output_line_graphs_folder: Union[Path, str],
                 cycles_array: List,
                 bc_cycles_array: List,
                 universal_len: int,
                 payload_len: int,
                 five_prime_len: int,
                 three_prime_len: int,
                 th_minimum_len_reads_to_analyse: int,
                 general_information_file: Union[Path, str],
                 count_reads_len_file: Union[Path, str],
                 algo_approach: str,
                 output_fasta_folder: Union[Path, str],
                 universals_fasta_format_file: Union[Path, str],
                 reads_chunk_to_fasta_format_file: Union[Path, str],
                 amount_of_universals: int,
                 blast_database_folder: Union[Path, str],
                 blast_db: Union[Path, str],
                 is_sampling_rate: bool,
                 sampling_rate_array: List,
                 output_all_sampling_rate_folder: Union[Path, str],
                 bc_list: List,
                 amount_of_cycles: int
                 ):
        self.input_file_or_folder = input_file_or_folder
        self.const_design_file = const_design_file
        self.payload_design_file = payload_design_file
        self.barcodes_design_file = barcodes_design_file
        self.results_good_reads_file = results_good_reads_file
        self.results_good_reads_with_len_bigger_then_y = results_good_reads_with_len_bigger_then_y
        self.len_reads_hist_output_file = len_reads_hist_output_file
        self.results_most_common_file = results_most_common_file
        self.design_simulation_file = design_simulation_file
        self.compare_design_to_experiment_results_output_file = compare_design_to_experiment_results_output_file
        self.foreach_bc_payload_count_file = foreach_bc_payload_count_file
        self.heatmap_foreach_bc_and_x_count_with_most_common_file = heatmap_foreach_bc_and_x_count_with_most_common_file
        self.count_reads_for_each_bc_file = count_reads_for_each_bc_file
        self.missing_bcs_file = missing_bcs_file
        self.output_hist_folder = output_hist_folder
        self.output_folder = output_folder
        self.design_folder = design_folder
        self.output_csv_folder = output_csv_folder
        self.output_graphs_folder = output_graphs_folder
        self.output_heatmap_folder = output_heatmap_folder
        self.hist_per_bc_file = hist_per_bc_file
        self.hist_foreach_bc_read_count_file = hist_foreach_bc_read_count_file
        self.hist_foreach_error_count_of_bc_file = hist_foreach_error_count_of_bc_file
        self.hist_foreach_read_count_count_bc_file = hist_foreach_read_count_count_bc_file
        self.len_reads_to_retrieve = len_reads_to_retrieve
        self.amount_of_bc = amount_of_bc
        self.barcode_len = barcode_len
        self.subset_size = subset_size
        self.amount_of_payloads = amount_of_payloads
        self.z_to_k_mer_representative = z_to_k_mer_representative
        self.k_mer_representative_to_z = k_mer_representative_to_z
        self.payload_pos = payload_pos
        self.uni_start_pos = uni_start_pos
        self.sampling_rate_from_good_reads_graph = sampling_rate_from_good_reads_graph
        self.output_line_graphs_folder = output_line_graphs_folder
        self.cycles_array = cycles_array
        self.bc_cycles_array = bc_cycles_array
        self.universal_len = universal_len
        self.payload_len = payload_len
        self.three_prime_len = three_prime_len
        self.five_prime_len = five_prime_len
        self.th_minimum_len_reads_to_analyse = th_minimum_len_reads_to_analyse
        self.design_results_only_z_file = design_results_only_z_file
        self.general_information_file = general_information_file
        self.count_reads_len_file = count_reads_len_file
        self.algo_approach = algo_approach
        self.output_fasta_folder = output_fasta_folder
        self.universals_fasta_format_file = universals_fasta_format_file
        self.reads_chunk_to_fasta_format_file = reads_chunk_to_fasta_format_file
        self.amount_of_universals = amount_of_universals
        self.blast_database_folder = blast_database_folder
        self.blast_db = blast_db
        self.sampling_rate_array = sampling_rate_array
        self.is_sampling_rate = is_sampling_rate
        self.output_all_sampling_rate_folder = output_all_sampling_rate_folder
        self.bc_list = bc_list
        self.amount_of_cycles = amount_of_cycles

    @uts.timer
    def run(self):
        self.create_folders()

        # upload design
        const_design_pd, payload_design_pd, barcodes_design_pd = self.upload_design()

        self.analyze_data_with_the_chosen_approach(algo_approach=self.algo_approach,
                                                   const_design_pd=const_design_pd,
                                                   payload_design_pd=payload_design_pd,
                                                   barcodes_design_pd=barcodes_design_pd)

        input_csv_path = self.results_good_reads_with_len_bigger_then_y
        # input_csv_path_sampling_percentage_100 = 'output_all_sampling_rate\\output_100\\csv\\results_good_reads_with_len_bigger_then_y.csv'
        for sampling_percentage in self.sampling_rate_array:
            for iteration in range(10):
                # if sampling_percentage == 100:
                #     input_csv_path_sampling_percentage_100 = 'output_all_sampling_rate/output_100/iter_0/csv/results_good_reads_with_len_bigger_then_y.csv'
                #     continue
                if sampling_percentage != 100:
                    input_csv_path = self.create_new_input_csv_according_to_sampling_rate(sampling_rate=sampling_percentage,
                                                                                          input_csv_path=input_csv_path_sampling_percentage_100,
                                                                                          iteration=iteration)
                # Find most common for each bc and for every cycle in that bc in results of good reads
                self.find_most_common(input_csv_path=input_csv_path,
                                      foreach_bc_payload_count_file_dict_to_csv=self.foreach_bc_payload_count_file,
                                      most_common_dict_to_csv_path=self.results_most_common_file,
                                      compare_design_to_experiment_results_output_file=self.compare_design_to_experiment_results_output_file,
                                      design_simulation_file=self.design_simulation_file,
                                      heatmap_foreach_bc_and_x_count_with_most_common_file=self.heatmap_foreach_bc_and_x_count_with_most_common_file)

                # For each bc count amount of reads sequenced
                self.hist_reads(csv_output_file=self.count_reads_for_each_bc_file,
                                input_csv_path=input_csv_path)

                # Create graph with sampling rate
                self.create_sampling_rate_from_good_reads_graph(input_csv_path=input_csv_path)

                # Move the data into the designated folder for each sampling rate
                output_folder = self.save_output_into_separate_folder(sampling_percentage=sampling_percentage, iteration=iteration)
                if sampling_percentage == 100:
                    input_csv_path_sampling_percentage_100 = output_folder + '/csv/results_good_reads_with_len_bigger_then_y.csv'

            # Delete all csvs from the output folder
            self.delete_csv_files_from_folder()

        self.graph_of_all_sampling_most_common_x_to_each_bc()

    ''' 
    identify_payload_by_pos_universals
    To identify the payloads, we only look at uni3-uni10
    '''

    def identify_payload_by_pos_universals(self, read: SeqIO, payload_design: pd.DataFrame, uni_pos_list: list,
                                           dist_option: str) -> \
            List[int]:
        res = [0, ] * len(self.payload_pos)
        # To identify the payloads, we only look at uni1
        for p_idx, pos in enumerate(uni_pos_list):
            if len(read) < (pos + self.payload_len):
                return res
            pos = pos + self.universal_len
            for s_idx, s in payload_design.iterrows():
                if dist_option == 'hamming':
                    if compare_with_errors(read[pos:pos + self.payload_len], s['Seq']):
                        res[p_idx] = s_idx
                        break
                else:  # Levenshtein
                    if lev(read[pos:pos + self.payload_len], s['Seq']) <= 6:
                        res[p_idx] = s_idx
                        break
        return res

    def identify_oligo_by_pos_universals(self, read: SeqIO, payload_design: pd.DataFrame, barcodes_design: pd.DataFrame,
                                         uni_pos_list: list, dist_option='hamming') -> List[int]:
        # identify_barcode
        # u2 location - bc len
        oligo_start_pos = uni_pos_list[0] - self.barcode_len
        bc = self.identify_barcode(read=read, barcodes_design=barcodes_design, oligo_start_pos=oligo_start_pos,
                                   dist_option=dist_option)
        # identify_payload
        payload = self.identify_payload_by_pos_universals(read=read, payload_design=payload_design,
                                                          uni_pos_list=uni_pos_list,
                                                          dist_option=dist_option)
        payload.insert(0, bc)
        return payload

    def identify_oligo(self, read: SeqIO, payload_design: pd.DataFrame, barcodes_design: pd.DataFrame,
                       oligo_start_pos=0, dist_option='hamming') -> List[int]:
        bc = self.identify_barcode(read=read, barcodes_design=barcodes_design, oligo_start_pos=oligo_start_pos,
                                   dist_option=dist_option)
        payload = self.identify_payload(read=read, payload_design=payload_design, oligo_start_pos=oligo_start_pos,
                                        dist_option=dist_option)
        payload.insert(0, bc)
        return payload

    # Identify which payload in each position
    def identify_payload(self, read: SeqIO, payload_design: pd.DataFrame, oligo_start_pos: int, dist_option: str) -> \
            List[int]:
        res = [0, ] * len(self.payload_pos)
        for p_idx, pos in enumerate(self.payload_pos):
            for s_idx, s in payload_design.iterrows():
                pos = pos + oligo_start_pos
                if dist_option == 'hamming':
                    if compare_with_errors(read[pos:pos + self.payload_len], s['Seq']):
                        res[p_idx] = s_idx
                        break
                else:  # Levenshtein
                    if lev(read[pos:pos + self.payload_len], s['Seq']) <= 6:
                        res[p_idx] = s_idx
                        break
        return res

    # Identify which barcode in each position
    def identify_barcode(self, read: SeqIO, barcodes_design: pd.DataFrame, oligo_start_pos: int,
                         dist_option: str) -> int:
        # pos = self.universal_len + self.five_prime_len + oligo_start_pos
        pos = oligo_start_pos
        for s_idx, s in barcodes_design.iterrows():
            if dist_option == 'hamming':
                if compare_with_errors(read[pos:pos + self.barcode_len], s['Seq']):
                    return s_idx
            else:  # Levenshtein
                if lev(read[pos:pos + self.barcode_len], s['Seq']) <= 6:
                    return s_idx

        return 0

    def verify_const_universal(self, read: SeqIO, const_design: pd.DataFrame, oligo_start_pos: int,
                               dist_option: str) -> bool:
        for i, (s_idx, s) in enumerate(const_design.iterrows()):
            pos = s['Pos'] + oligo_start_pos
            if dist_option == 'hamming':
                if not compare_with_errors(read[pos:pos + self.universal_len], s['Seq']):
                    return False
            else:  # Levenshtein
                if not lev(read[pos:pos + self.universal_len], s['Seq']) <= 6:
                    return False
        return True

    def verify_const_universal_and_reverse_complement(self, read: SeqIO, const_design: pd.DataFrame, oligo_start_pos=0,
                                                      dist_option='hamming') -> SeqIO:
        if self.verify_const_universal(read=read, const_design=const_design, oligo_start_pos=oligo_start_pos,
                                       dist_option=dist_option):
            return read
        return SeqIO.SeqRecord("None")

    # def reads_chunk_to_fasta_format(self, reads, idx) -> str:
    #     # Write the sequences to the output file in FASTA format
    #     output_fasta_file = self.reads_chunk_to_fasta_format_file + str(idx) + '.fasta'
    #     with open(output_fasta_file, "w") as output_handle:
    #         SeqIO.write(reads.values(), output_handle, "fasta")
    #
    #     return output_fasta_file

    def reads_chunk_to_fasta_format(self, reads, idx) -> str:
        # Write the sequences and their reverse complements to the output file in FASTA format
        output_fasta_file = self.reads_chunk_to_fasta_format_file + str(idx) + '.fasta'
        with open(output_fasta_file, "w") as output_handle:
            for read_id, read_seq in reads.items():
                # Write the original sequence
                original_seq_record = SeqRecord(read_seq.seq,
                                                id=read_id,
                                                description="original")
                SeqIO.write(original_seq_record, output_handle, "fasta")

                # Generate and write the reverse complement
                reverse_complement_seq = str(read_seq.seq.reverse_complement())
                reverse_complement_record = SeqRecord(Seq(reverse_complement_seq),
                                                      id=read_id,
                                                      description="reverse_complement")
                SeqIO.write(reverse_complement_record, output_handle, "fasta")

        return output_fasta_file

    def create_blast_database(self):
        makeblastdb_cmd = f"makeblastdb -in {self.universals_fasta_format_file} -dbtype nucl -out {self.blast_db}"
        try:
            subprocess.run(makeblastdb_cmd, shell=True, check=True)
            print(f"BLAST database '{self.blast_db}' created successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error creating BLAST database: {e}")

    def convert_design_file_to_fasta(self):

        # Open the CSV file and read the data
        with open(self.const_design_file, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  # skip the header row
            sequences = []
            for row in csvreader:
                # Extract the sequence ID, sequence, and position from the CSV row
                seq_id = row[0]
                seq = row[1]
                pos = int(row[2])

                # Add the sequence to the list of sequences
                sequences.append((seq_id, seq, pos))

        # Write the sequences to a multi-FASTA file
        with open(self.universals_fasta_format_file, 'w') as outfile:
            for seq_id, seq, pos in sequences:
                # Write the sequence header in FASTA format
                header = f'>{seq_id} position={pos}\n'
                outfile.write(header)

                # Write the sequence in FASTA format, breaking it into 80-character lines
                seq_lines = [seq[i:i + 80] for i in range(0, len(seq), 80)]
                outfile.write('\n'.join(seq_lines) + '\n')

    @uts.timer
    def run_blastn_to_find_location_of_universals_in_reads_chunks(self, reads_chunk_i: int, blast_database_name: str,
                                                                  output_fasta_file: str) -> Dict:
        """qseqid: Query Seq-id
            sseqid: Subject Seq-id
            pident: Percentage of identical matches
            length: Alignment length
            mismatch: Number of mismatches
            gapopen: Number of gap openings
            qstart: Start of alignment in query
            qend: End of alignment in query
            sstart: Start of alignment in subject
            send: End of alignment in subject
            evalue: Expect value
            bitscore: Bit score
            btop: Blast traceback operations"""

        blastn_cmd = 'blastn -query ' + output_fasta_file + \
                     ' -db ' + self.blast_db + \
                     ' -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore btop qseq sseq" -word_size 7 -gapopen 5 -gapextend 2 -reward 1 -penalty -3 -min_raw_gapped_score 18 -num_threads 6'

        # Run the command and capture the output
        try:
            query_results = subprocess.check_output(blastn_cmd, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e.output.decode())

        blastn_output_idx_file = self.blast_database_folder + 'blastn_output_' + str(reads_chunk_i) + '.txt'
        with open(blastn_output_idx_file, "wb") as binary_file:

            # Write bytes to file
            binary_file.write(query_results)
        query_results = {}
        with open(blastn_output_idx_file) as blastn_output_file:
            for line in blastn_output_file:
                fields = line.strip().split("\t")
                read_id = fields[0]
                uni_id = fields[1]
                # percentage_of_identical_matches = fields[2]
                alignment_length = fields[3]
                number_of_mismatches = fields[4]
                number_of_gap_openings = fields[5]
                start_of_alignment_in_read = fields[6]
                end_of_alignment_in_read = fields[7]
                # start_of_alignment_in_uni = fields[8]
                # end_of_alignment_in_uni = fields[9]
                # expect_value = fields[10]
                bit_score = fields[11]
                # blast_traceback_operations = fields[12]
                # aligned_part_of_query_sequence = fields[13]
                # aligned_part_of_subject_sequence = fields[14]
                # alignment_location = f"{uni_id}:{start_of_alignment_in_read}-{end_of_alignment_in_read}"
                if read_id not in query_results:
                    query_results[read_id] = {}
                    query_results[read_id][uni_id] = []
                if uni_id not in query_results[read_id]:
                    query_results[read_id][uni_id] = []
                query_results[read_id][uni_id].append((start_of_alignment_in_read, end_of_alignment_in_read, bit_score,
                                                       alignment_length, number_of_mismatches, number_of_gap_openings))

        return query_results

    def create_new_input_csv_according_to_sampling_rate(self, sampling_rate: float,
                                                        input_csv_path: str,
                                                        iteration: int) -> str:
        # Set the chunk size and sampling fraction
        chunk_size = 100000  # 100,000 rows per chunk
        sampling_fraction = float(sampling_rate / 100)

        # Create a directory to store the sampled CSV files
        output_sampling_rate_dir = self.output_all_sampling_rate_folder + f'output_{sampling_rate}/iter_{iteration}/'
        os.makedirs(output_sampling_rate_dir, exist_ok=True)

        # Iterate over the CSV file in chunks
        sampled_filepaths = []
        for i, chunk in enumerate(pd.read_csv(input_csv_path, chunksize=chunk_size)):
            # Randomly sample rows from the chunk with the specified fraction
            sampled_chunk = chunk.sample(frac=sampling_fraction)
            # Write the sampled chunk to a new CSV file
            chunk_filename = f'sampled_chunk_{i}_iter_{iteration}.csv'
            chunk_filepath = os.path.join(output_sampling_rate_dir, chunk_filename)
            sampled_chunk.to_csv(chunk_filepath, index=False)
            sampled_filepaths.append(chunk_filepath)

        # Concatenate the sampled files into a single output file
        output_filepath = os.path.join(output_sampling_rate_dir,
                                       f'sampled_{sampling_rate}_iter_{iteration}.csv')
        write_header = True
        with open(output_filepath, 'w+') as output_file:
            for filepath in sampled_filepaths:
                with open(filepath, 'r') as input_file:
                    if write_header:
                        output_file.write(input_file.readline())  # Write the header
                        write_header = False
                    else:
                        input_file.readline()  # Skip the header

                    # Copy the rest of the file content
                    shutil.copyfileobj(input_file, output_file)
                os.remove(filepath)

        # Change the permissions of the output file to 0o777
        os.chmod(output_filepath, 0o777)

        return output_filepath

    def find_u2_location_in_read(self, read: SeqIO, u2: str) -> (int, bool):
        # TODO: what should be the maximum dist to find the universal?
        max_distance = 6
        best_location = 0
        best_distance = float('inf')
        is_found_u2 = False

        for i in range(len(read) - len(u2) + 1):
            substring = read[i:i + len(u2)]
            distance = lev(substring, u2)
            if distance <= max_distance:
                if distance < best_distance:
                    best_distance = distance
                    best_location = i
        if best_location != 0:
            is_found_u2 = True
            return (best_location - self.barcode_len), is_found_u2
        return best_location, is_found_u2

    def custom_sort_key(self, key):
        # Get the numeric suffix of the key
        num = int(key.replace('Universal', ''))
        # Return the numeric suffix as the sorting key
        return num

    def convert_uni_pos_in_read_to_list(self, unis_pos_in_read: Dict) -> List[int]:
        unis_pos_in_read_sorted = collections.OrderedDict(
            sorted(unis_pos_in_read.items(), key=lambda t: self.custom_sort_key(t[0]), reverse=False))
        experiment_uni_pos = [0] * self.amount_of_universals
        best_uni_info_prev = -1
        for uni_name, uni_info in unis_pos_in_read_sorted.items():
            # experiment_uni_pos[]
            uni_idx = int(uni_name.replace('Universal', ''))
            best_bit_score = -1
            best_uni_info = -1
            for info in uni_info:
                if float(info[2]) > best_bit_score:
                    best_bit_score = float(info[2])
                    best_uni_info = info

            if best_uni_info_prev < int(best_uni_info[0]):
                experiment_uni_pos[uni_idx - 1] = int(best_uni_info[0])
                best_uni_info_prev = int(best_uni_info[0])

        # blast give the positions with +1, so we fix it by subtracting it by 1
        experiment_uni_pos_blast_fix_minus_1 = [x - 1 if x != 0 else x for x in experiment_uni_pos]
        return experiment_uni_pos_blast_fix_minus_1

    def merge_design_pos_and_pos_in_read_list(self, experiment_uni_pos: List[int]) -> List[int]:
        merged_design_uni_pos_and_pos_in_read_list = [-1] * self.amount_of_universals
        if experiment_uni_pos[0] == 0:
            merged_design_uni_pos_and_pos_in_read_list[0] = self.uni_start_pos[0]
        else:
            merged_design_uni_pos_and_pos_in_read_list[0] = experiment_uni_pos[0]

        i = 1
        for experiment_uni_pos_i, design_uni_start_pos in zip(experiment_uni_pos[1:], self.uni_start_pos[1:]):
            if experiment_uni_pos_i == 0:
                merged_design_uni_pos_and_pos_in_read_list[i] = merged_design_uni_pos_and_pos_in_read_list[i - 1] + (
                        self.uni_start_pos[i] - self.uni_start_pos[i - 1])
            else:
                merged_design_uni_pos_and_pos_in_read_list[i] = experiment_uni_pos_i
            i += 1

        return merged_design_uni_pos_and_pos_in_read_list

    def find_pos_in_read_list(self, unis_pos_in_read: Dict) -> List[int]:
        experiment_uni_pos_list = self.convert_uni_pos_in_read_to_list(unis_pos_in_read=unis_pos_in_read)
        merged_design_uni_pos_and_pos_in_read_list = self.merge_design_pos_and_pos_in_read_list(
            experiment_uni_pos=experiment_uni_pos_list)
        return merged_design_uni_pos_and_pos_in_read_list

    def find_unis_locations_in_read(self, read: SeqIO, const_design: pd.DataFrame) -> Dict:
        # TODO: what should be the maximum dist to find the universal?
        max_distance = 6
        best_distance = float('inf')
        prev_best_location = 0
        curr_best_location = prev_best_location
        uni_pos_dict = {}

        for uni_name, uni_seq in const_design['Seq'].items():
            for i in range(curr_best_location, len(read) - len(uni_seq) + 1):
                substring = read[i:i + len(uni_seq)]
                distance = lev(substring, uni_seq)
                if distance <= max_distance and distance < best_distance:
                    best_distance = distance
                    curr_best_location = i
            if curr_best_location != prev_best_location:
                uni_pos_dict[uni_name] = (curr_best_location, True)
                prev_best_location = curr_best_location

        return uni_pos_dict

    def extract_start_position_and_reads_results_to_csv(self, reads: List[str],
                                                        const_design: pd.DataFrame,
                                                        payload_design: pd.DataFrame,
                                                        barcodes_design: pd.DataFrame,
                                                        dist_option: str,
                                                        output_csv_path: Path) -> None:
        reads = self.retrieve_reads_in_specific_len_at_least_list(reads=reads,
                                                                  length=self.th_minimum_len_reads_to_analyse)
        res = list()
        failed = 0
        read_idx = 0
        none = SeqIO.SeqRecord("None")
        u2 = const_design.loc[const_design.index == "Universal1", "Seq"]['Universal1']
        uni_pos_list = [0] * self.amount_of_universals
        count_found_u2 = 0
        with open(output_csv_path, "ab") as f:
            cols_names = [self.bc_cycles_array]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        for read_idx, read in enumerate(reads):
            if read_idx % 1000 == 999:
                uts.write_list_to_csv(
                    ['processed:', read_idx + 1, 'reads:', failed, str(100 * failed / (read_idx + 1)) + '%', 'failed'],
                    self.general_information_file)
                with open(output_csv_path, "ab") as f:
                    np.savetxt(f, res, fmt='%i', delimiter=",")
                res = list()
            u2_location, is_found_u2 = self.find_u2_location_in_read(read=read, u2=u2)
            read = self.verify_const_universal_and_reverse_complement(read=read,
                                                                      const_design=const_design,
                                                                      oligo_start_pos=u2_location,
                                                                      dist_option=dist_option)
            if is_found_u2:
                count_found_u2 += 1
            if read.seq == none.seq:
                failed += 1
                continue

            res.append(self.identify_oligo(read=read,
                                           payload_design=payload_design,
                                           barcodes_design=barcodes_design,
                                           oligo_start_pos=u2_location,
                                           dist_option=dist_option))
            if res[-1].__contains__(0):
                failed += 1

        uts.write_list_to_csv(
            ['processed:', read_idx + 1, 'reads:', failed, str(100 * failed / (read_idx + 1)) + '%', 'failed'],
            self.general_information_file)
        uts.write_list_to_csv(
            ['count_found_u2:', count_found_u2],
            self.general_information_file)
        with open(output_csv_path, "ab") as f:
            np.savetxt(f, res, fmt='%i', delimiter=",")

    def process_fastq_files(self):
        if os.path.isfile(self.input_file_or_folder) and self.input_file_or_folder.endswith('.fastq'):
            fastq_files = [Path(self.input_file_or_folder)]
        elif os.path.isdir(self.input_file_or_folder):
            fastq_files = [os.path.join(Path(self.input_file_or_folder), Path(file_path)) for file_path in
                           os.listdir(self.input_file_or_folder) if
                           file_path.endswith('.fastq')]
        else:
            raise ValueError("The input must be a .fastq file or a directory containing .fastq files.")

        return fastq_files

    def prepare_run_extract_all_pos_and_reads_results_to_csv(self, output_csv_path):
        self.convert_design_file_to_fasta()
        blast_database_name = self.blast_database_folder + "blast_database"
        self.create_blast_database()

        is_run_loop = True
        reads_length_counts = {'length': 'count'}

        threads = 1

        with open(output_csv_path, "ab") as f:
            cols_names = [self.bc_cycles_array]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        number_of_reads = 0

        return number_of_reads, is_run_loop, reads_length_counts, blast_database_name

    @uts.timer
    def extract_all_pos_and_reads_results_to_csv(self, payload_design, barcodes_design,
                                                 dist_option, output_csv_path, number_of_reads, is_run_loop,
                                                 reads_length_counts, blast_database_name):

        chunk_size = 1000
        chunk_end = chunk_size
        reads_chunk_i = 1
        res = list()
        failed = 0

        number_of_reads += uts.get_amount_of_reads_from_file(file_path=self.input_file_or_folder)
        uts.write_list_to_csv(['# reads', number_of_reads], self.general_information_file)

        reads_iter = uts.open_fastq_yield(self.input_file_or_folder)

        ''' Usning blastn to analyze the reads '''
        while is_run_loop:
            if number_of_reads < chunk_end:
                chunk_size = number_of_reads % chunk_size
                is_run_loop = False
            print(chunk_size)
            reads_chunk = list(itertools.islice(reads_iter, chunk_size))
            reads_specific_len, reads_length_counts = self.retrieve_reads_in_specific_len_at_least_dict(
                reads=reads_chunk,
                length=self.th_minimum_len_reads_to_analyse, reads_length_counts=reads_length_counts)
            if len(reads_specific_len) == 0:
                if not is_run_loop:
                    break
                chunk_end += chunk_size
                reads_chunk_i += 1
                continue
            fasta_chunk_idx = chunk_end - (chunk_end - (reads_chunk_i * chunk_size)) + chunk_size
            output_fasta_file = self.reads_chunk_to_fasta_format(reads=reads_specific_len, idx=fasta_chunk_idx)

            query_results = self.run_blastn_to_find_location_of_universals_in_reads_chunks(reads_chunk_i=reads_chunk_i,
                                                                                           blast_database_name=blast_database_name,
                                                                                           output_fasta_file=output_fasta_file)

            for read_idx, (read_id, unis_pos_in_read) in enumerate(query_results.items()):
                read = reads_specific_len[read_id]
                uni_pos_list = self.find_pos_in_read_list(unis_pos_in_read=unis_pos_in_read)

                res.append(self.identify_oligo_by_pos_universals(read=read,
                                                                 payload_design=payload_design,
                                                                 barcodes_design=barcodes_design,
                                                                 uni_pos_list=uni_pos_list,
                                                                 dist_option=dist_option))
                if res[-1].__contains__(0):
                    failed += 1

            uts.write_list_to_csv(
                ['processed reads:', chunk_end, 'failed reads:', failed, str(100 * failed / (chunk_end + 1)) + '%',
                 'failed.', 'query_results: ', len(query_results)],
                self.general_information_file)
            with open(output_csv_path, "ab") as f:
                np.savetxt(f, res, fmt='%i', delimiter=",")
            res = list()

            if not is_run_loop:
                break
            chunk_end += chunk_size
            reads_chunk_i += 1

        uts.write_dict_to_csv_as_dict(reads_length_counts, self.count_reads_len_file)
        self.hist_length_counts_reads()

    def most_frequency_value_in_reads(self, lens: List[int]) -> None:
        length_counts = {}
        for len in lens:
            if len in length_counts:
                length_counts[len] += 1
            else:
                length_counts[len] = 1

        print(length_counts)
        most_common_length = max(length_counts, key=length_counts.get)
        print(
            f'most_frequency_value_in_reads: {most_common_length} with {length_counts[most_common_length]} appearances')

        # Write general information to csv
        uts.write_list_to_csv(
            ['most_frequency_value_in_reads:', most_common_length, 'appearances:', length_counts[most_common_length]],
            self.general_information_file)

    def reads_len_hist(self, reads: List[str]) -> None:
        len_reads = len(reads)
        uts.write_list_to_csv(['# reads', len_reads], self.general_information_file)
        print(f'# reads: {len_reads}')

        lens = [len(r) for r in reads]
        plt.hist(lens, bins=2000)
        plt.xlabel('length')
        plt.ylabel('# reads')
        plt.savefig(self.len_reads_hist_output_file)
        plt.close()

        # most frequency value in reads
        self.most_frequency_value_in_reads(lens=lens)

    def upload_design(self) -> Tuple[Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame]]:
        const_design_pd = pd.read_csv(self.const_design_file, index_col=0, header=0)
        payload_design_pd = pd.read_csv(self.payload_design_file, index_col=0, header=0)
        barcodes_design_pd = pd.read_csv(self.barcodes_design_file, index_col=0, header=0)

        return const_design_pd, payload_design_pd, barcodes_design_pd

    def retrieve_reads_in_specific_len(self, reads: List[str], length: int) -> List[str]:
        reads = [r for r in reads if len(r) == length]
        print(f'{len(reads)} reads in len {length}')
        return reads

    def retrieve_reads_in_specific_len_at_least_list(self, reads: List[str], length: int) -> List[str]:
        reads = [r for r in reads if len(r) >= length]
        print(f'{len(reads)} reads in len at least {length}')
        return reads

    def retrieve_reads_in_specific_len_at_least_dict(self, reads: List[SeqIO.SeqRecord], length: int,
                                                     reads_length_counts: Dict) -> Tuple[Dict, Dict]:
        reads_specific_len = {}
        for r in reads:
            read_length = len(r)
            if read_length >= length:
                reads_specific_len[r.id] = r

                if read_length in reads_length_counts:
                    reads_length_counts[read_length] += 1
                else:
                    reads_length_counts[read_length] = 1
        print(f'{len(reads_specific_len)} reads in len at least {length}')

        return reads_specific_len, reads_length_counts

    def reads_results_to_csv(self, reads: List[str],
                             const_design: pd.DataFrame,
                             payload_design: pd.DataFrame,
                             barcodes_design: pd.DataFrame,
                             output_csv_path: Union[Path, str]) -> None:
        res = list()
        failed = 0
        read_idx = 0
        none = SeqIO.SeqRecord("None")
        u2 = const_design.loc[const_design.index == "Universal2", "Seq"]['Universal2']

        with open(self.output_csv_path, "ab") as f:
            cols_names = [self.bc_cycles_array]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        for read_idx, read in enumerate(reads):
            if read_idx % 1000 == 999:
                print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
                with open(self.output_csv_path, "ab") as f:
                    np.savetxt(f, res, fmt='%i', delimiter=",")
                res = list()

            read = self.verify_const_universal_and_reverse_complement(read=read, const_design=const_design)
            if read.seq == none.seq:
                failed += 1
                continue

            res.append(self.identify_oligo(read=read, payload_design=payload_design, barcodes_design=barcodes_design))
            if res[-1].__contains__(0):
                failed += 1

        print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
        with open(output_csv_path, "ab") as f:
            np.savetxt(f, res, fmt='%i', delimiter=",")

    def hist_per_bc(self, dict_bc_payload, bc, payload):
        plt.bar(dict_bc_payload.keys(), dict_bc_payload.values())
        amount_of_payloads = self.amount_of_payloads + 1
        plt.xticks(range(amount_of_payloads))
        plt.xlim(0.5, amount_of_payloads)
        plt.title('bc=' + str(bc) + 'payload=' + payload)
        uts.is_dir_exists(self.hist_per_bc_file)
        plt.savefig(self.hist_per_bc_file + '/bc=' + str(bc) + '_payload=' + payload + '_hist.png')
        plt.close()

    def analyze_results_good_reads(self, input_csv_path: Union[Path, str],
                                   dict_to_csv: Union[Path, str]) -> Dict:
        df = pd.read_csv(input_csv_path)
        dict_bc = {}

        amount_of_bc = self.amount_of_bc + 1
        for i in range(amount_of_bc):
            dict_bc_i = {}
            for payload in self.cycles_array:
                payload_dict = {}
                for p in range(self.amount_of_payloads + 1):
                    payload_dict[p] = 0
                dict_p = {payload: payload_dict}
                dict_bc_i.update(dict_p)
            dict_bc[i] = dict_bc_i
        df_bc = df.sort_values(self.bc_cycles_array, ascending=True, key=np.sin)
        df_bc.fillna(0)
        for row_idx, row in tqdm(df_bc.iterrows()):
            for location_payload, payload in row[1:].items():
                dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

        uts.write_dict_to_csv(dict_bc, dict_to_csv)

        return dict_bc

    @uts.timer
    def most_common_for_each_bc(self, dict_bc: Dict, dict_to_csv_path: Union[Path, str]) -> None:
        dict_bc_most_common = {}
        amount_of_bc = self.amount_of_bc + 1

        # init dict_bc_most_common_i
        for i in range(amount_of_bc):
            dict_bc_most_common_i = {}
            for payload in self.cycles_array:
                dict_p = {payload: []}
                dict_bc_most_common_i.update(dict_p)
            dict_bc_most_common[i] = dict_bc_most_common_i

        # find the config['subset_size'] most common
        for bc_i in range(1, amount_of_bc):
            for payload in self.cycles_array:
                del dict_bc[bc_i][payload][0]

                most_common = heapq.nlargest(self.subset_size, dict_bc[bc_i][payload],
                                             key=dict_bc[bc_i][payload].get)
                dict_bc_most_common[bc_i][payload] = most_common

                print(f'bc = {bc_i}, cycle = {payload}, {self.subset_size} most common = {most_common}')
                self.hist_per_bc(dict_bc_payload=dict_bc[bc_i][payload], bc=bc_i, payload=payload)

        uts.write_dict_to_csv(dict=dict_bc_most_common, csv_path=dict_to_csv_path)

    def convert_most_common_to_letters_in_new_alphabet(self, results_most_common_file: Union[Path, str]) -> tuple[
        dict[int, list[Any]], dict[int, list[Any]]]:
        df = pd.read_csv(results_most_common_file, index_col=0)
        dict_most_common = df.to_dict("list")
        result_payload = {}
        result_only_z = {}
        print(dict_most_common)

        for bc_i in range(1, (self.amount_of_bc + 1)):
            result_payload[bc_i] = []
            result_only_z[bc_i] = []
            for payload in self.cycles_array:
                bc_i_str = str(bc_i)
                d = ast.literal_eval(dict_most_common[bc_i_str][0])
                d = d[payload]
                dict_convert_to_x = []
                for i in d:
                    dict_convert_to_x.append('X' + str(i))
                dict_convert_to_x = tuple(uts.sorted_human(dict_convert_to_x))

                try:
                    z = self.k_mer_representative_to_z[dict_convert_to_x]
                except KeyError:
                    z = 'Z1'  # The X tuple is out of range
                result_payload[bc_i].append({z: dict_convert_to_x})
                result_only_z[bc_i].append(z)

        return result_payload, result_only_z

    @uts.timer
    def compare_most_common_to_design(self, result_payload: Dict,
                                      compare_design_to_experiment_results_output_file: Union[Path, str],
                                      design_simulation_file: Union[Path, str]) -> None:
        with open(compare_design_to_experiment_results_output_file, "ab") as f:
            cols_names = [
                ['bc', 'design_z', 'design_x', 'experiment_results_z', 'experiment_results_x', 'is_retrieved_all_z',
                 'count_all_z_mm']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        df_design_simulation_file = pd.read_csv(design_simulation_file)
        df_design_simulation_T_file = df_design_simulation_file.T
        # with open(design_simulation_file, 'r') as inf:
        for bc_i_idx, bc_i_design in enumerate(df_design_simulation_T_file.iterrows()):
            if bc_i_idx == 0:
                continue
            count_all_z_mm = 0
            bc_i_design = bc_i_design[1]
            for z_i_design, z_i_results in zip(bc_i_design, result_payload[bc_i_idx]):
                if z_i_design != list(z_i_results)[0]:
                    count_all_z_mm += 1
            with open(compare_design_to_experiment_results_output_file, "ab") as f:
                design_z = " ".join(bc_i_design)
                design_x = str(
                    [{z: self.z_to_k_mer_representative[z]} for z in bc_i_design]).replace(
                    ",",
                    " ")
                experiment_results_z = " ".join([list(z.keys())[0] for z in result_payload[bc_i_idx]])
                experiment_results_x = str(result_payload[bc_i_idx]).replace(",", " ")
                is_retrieved_all_z = (count_all_z_mm == 0)
                count_all_z_mm = count_all_z_mm
                cols = [[bc_i_idx, design_z, design_x, experiment_results_z, experiment_results_x, is_retrieved_all_z,
                         count_all_z_mm]]
                np.savetxt(f, cols, fmt='%s', delimiter=",")

    def create_heatmap_with_rectangles_on_most_common(self, dict_foreach_bc_and_x_count_all_cycles_matrix: Dict,
                                                      heatmap_foreach_bc_and_x_count_with_most_common_file: Union[
                                                          Path, str],
                                                      ax) -> None:

        # remove col 0 and row 0
        del dict_foreach_bc_and_x_count_all_cycles_matrix[0]
        dict_foreach_bc_and_x_count_all_cycles_matrix.drop(index=['0'])

        # for each bc, normalize the data to 1
        dict_foreach_bc_and_x_count_all_cycles_matrix = dict_foreach_bc_and_x_count_all_cycles_matrix.div(
            dict_foreach_bc_and_x_count_all_cycles_matrix.sum(axis=1), axis=0)

        # fill nan cells with 0
        dict_foreach_bc_and_x_count_all_cycles_matrix = dict_foreach_bc_and_x_count_all_cycles_matrix.fillna(0)

        # create heatmap
        sns.heatmap(dict_foreach_bc_and_x_count_all_cycles_matrix.T)
        plt.xlabel("payloads")
        plt.ylabel("bc")
        # Add lines to distinguish payloads at each cycle
        ax.hlines(
            [self.amount_of_payloads, (self.amount_of_payloads * 2), (self.amount_of_payloads * 3)],
            *ax.get_xlim())

        df_results_most_common = pd.read_csv(self.results_most_common_file)
        del df_results_most_common['0']
        for bc_idx, cycles in df_results_most_common.items():
            if bc_idx == '0':
                continue
            dict_cycles = ast.literal_eval(cycles[0])
            cycle_idx = 1
            for cycle_name, cycle_data in dict_cycles.items():
                for one_of_most_common in cycle_data:
                    ax.add_patch(
                        patches.Rectangle(
                            (int(bc_idx),
                             ((self.amount_of_payloads * (cycle_idx - 1)) + int(one_of_most_common)) - 1),
                            1.0,
                            1.0,
                            edgecolor='red',
                            fill=False,
                            lw=0.5
                        ))
                cycle_idx += 1

        plt.savefig(heatmap_foreach_bc_and_x_count_with_most_common_file, dpi=400)
        plt.close()

    def heatmap_foreach_bc_and_x_count_with_most_common(self,
                                                        heatmap_foreach_bc_and_x_count_with_most_common_file: Union[
                                                            Path, str]) -> None:
        df = pd.read_csv(self.foreach_bc_payload_count_file)
        dict_foreach_bc_and_x_count_str = df.to_dict("list")

        dict_foreach_bc_and_x_count_split_to_cycles = {}
        dict_foreach_bc_and_x_count_all_cycles = {}
        for bc_i_str, payloads in dict_foreach_bc_and_x_count_str.items():
            dict_foreach_bc_and_x_count_all_cycles[bc_i_str] = {}
            dict_foreach_bc_and_x_count_split_to_cycles[bc_i_str] = ast.literal_eval(
                dict_foreach_bc_and_x_count_str[bc_i_str][0])

            cycle_idx = 1
            for cycle_name, cycle_data in dict_foreach_bc_and_x_count_split_to_cycles[bc_i_str].items():
                cycle_data_new = cycle_data
                if cycle_name != 'c1':
                    cycle_idx += 1
                    cycle_data_new = {}
                    for payload_idx, payload_count in cycle_data.items():
                        if payload_idx == 0:
                            continue
                        payload_new_x = (self.amount_of_payloads * (cycle_idx - 1)) + payload_idx
                        cycle_data_new[payload_new_x] = payload_count

                dict_foreach_bc_and_x_count_all_cycles[bc_i_str].update(cycle_data_new)

        fig, ax = plt.subplots(figsize=(10, 7))

        # make the matrix bcs x payloads
        dict_foreach_bc_and_x_count_all_cycles_matrix = pd.DataFrame(dict_foreach_bc_and_x_count_all_cycles).T.fillna(0)

        self.create_heatmap_with_rectangles_on_most_common(dict_foreach_bc_and_x_count_all_cycles_matrix=
                                                           dict_foreach_bc_and_x_count_all_cycles_matrix,
                                                           heatmap_foreach_bc_and_x_count_with_most_common_file=heatmap_foreach_bc_and_x_count_with_most_common_file,
                                                           ax=ax)

    @uts.timer
    def find_most_common(self, input_csv_path: Union[Path, str],
                         foreach_bc_payload_count_file_dict_to_csv: Union[Path, str],
                         most_common_dict_to_csv_path: Union[Path, str],
                         compare_design_to_experiment_results_output_file: Union[Path, str],
                         design_simulation_file: Union[Path, str],
                         heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str]) -> None:
        dict_bc = self.analyze_results_good_reads(input_csv_path=input_csv_path,
                                                  dict_to_csv=foreach_bc_payload_count_file_dict_to_csv)
        self.most_common_for_each_bc(dict_bc=dict_bc, dict_to_csv_path=most_common_dict_to_csv_path)
        result_payload, result_only_z = self.convert_most_common_to_letters_in_new_alphabet(
            results_most_common_file=most_common_dict_to_csv_path)
        self.compare_most_common_to_design(result_payload=result_payload,
                                           compare_design_to_experiment_results_output_file=compare_design_to_experiment_results_output_file,
                                           design_simulation_file=design_simulation_file)

        self.heatmap_foreach_bc_and_x_count_with_most_common(
            heatmap_foreach_bc_and_x_count_with_most_common_file=heatmap_foreach_bc_and_x_count_with_most_common_file)

    def missing_bc_to_csv(self, dict_append_missing_bc):
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        ser_append_missing_bc.to_csv(self.missing_bcs_file, mode='a', header=False)

    def create_sampling_rate_from_good_reads_graph(self, input_csv_path: Path) -> None:
        """
        This function takes in a list of integers and creates a graph with the x-axis being the sampling rate
        and the y-axis being the count of different values in arr[0].
        The sampling rate will be in increments of 10% starting from 0% until 100%.
        """

        df_good_reads = pd.read_csv(input_csv_path)

        sampling_rates = [i / 10 for i in range(11)]  # create a list of sampling rates from 0% to 100%
        counts = []  # list to store counts of different values in arr[0]
        # Iterate over the sampling rates
        for sampling_rate in sampling_rates:
            # Sample the dataframe with the current sampling rate
            sampled_df = df_good_reads.sample(frac=sampling_rate)

            # Count the number of different values in the "bc" column
            count = sampled_df["bc"].nunique()

            # Add the count to the list
            counts.append(count)

        # Plot the graph
        plt.plot(sampling_rates, counts)
        plt.title("sampling_rate_graph")
        plt.xlabel('Sampling Rate %')
        plt.ylabel('Count of Different Values')
        plt.savefig(self.sampling_rate_from_good_reads_graph)
        plt.close()

    def for_each_bc_count_reads_read_to_csv(self, output_file: Union[Path, str], input_csv_path: Path) -> pd.DataFrame:
        with open(output_file, "ab") as f:
            cols_names = [['bc', 'count']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        df = pd.read_csv(input_csv_path)
        forth_column = df.iloc[:, 0]
        counts = forth_column.value_counts()
        count_sorted = counts.sort_index()

        dict_append_missing_bc = {}
        for bc_idx in range(1, self.amount_of_bc + 1):
            if bc_idx not in count_sorted or count_sorted[bc_idx] == 0:
                dict_append_missing_bc[bc_idx] = 0

        self.missing_bc_to_csv(dict_append_missing_bc)
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        count_sorted = count_sorted.append(ser_append_missing_bc)
        count_sorted_with_missing_bc = count_sorted.sort_index()
        count_sorted_with_missing_bc.to_csv(output_file, mode='a', header=False)

        return count_sorted_with_missing_bc

    def calculate_average_coverage(self, csv_file: Union[Path, str]):
            bc_count_df = pd.read_csv(csv_file)

            # filter the DataFrame
            filtered_df = bc_count_df[bc_count_df['bc'].isin(self.bc_list)]

            # calculate the average of the 'count' column
            avg_count = np.mean(filtered_df['count'])

            return avg_count

    def process_csv(self, file_path: Union[Path, str]):
        counts = []
        with open(file_path, "r") as csvfile:
            csvreader = csv.DictReader(csvfile)
            for row_idx, row in enumerate(csvreader):
                if (row_idx + 1) in self.bc_list:
                    design_x = self.parse_custom_format(row['design_x'])
                    experiment_results_x = self.parse_custom_format(row['experiment_results_x'])
                    count_correct_xi = self.compare_rows(design_x, experiment_results_x)
                    counts.append(count_correct_xi)

        return counts

    def compare_rows(self, row1, row2):
        count_correct_xi = []
        for (key1, xi_1), (key2, xi_2) in zip(row1, row2):
            common_xi = set(xi_1).intersection(set(xi_2))
            count_correct_xi.append(len(common_xi))
        return count_correct_xi

    def parse_custom_format(self, input_str: str) -> Dict:
        result = []
        for item_str in re.split("}|{", input_str)[1:-1]:
            item_str = item_str.strip()
            if not item_str:
                continue

            key, values_str = item_str.split(":")
            key = key.strip('\'')
            values = tuple(values_str.strip(" ()").replace("'", "").split())
            result.append((key, values))

        return result

    def hist_foreach_bc_read_count(self, csv_output_file: Union[Path, str]) -> None:
        count_sorted_with_missing_bc = pd.read_csv(csv_output_file)
        x = count_sorted_with_missing_bc["count"].index
        y = count_sorted_with_missing_bc["count"]
        plt.bar(x, y, align='center')
        plt.xlabel('bc')
        plt.ylabel('count')
        plt.savefig(self.hist_foreach_bc_read_count_file)
        plt.close()

    def hist_foreach_error_count_of_bc(self) -> None:
        count_sorted_with_missing_bc = pd.read_csv(self.compare_design_to_experiment_results_output_file)
        x = count_sorted_with_missing_bc["count_all_z_mm"]
        plt.hist(x, bins=len(self.cycles_array))
        plt.xlabel('errors in bc')
        plt.ylabel('count')
        plt.savefig(self.hist_foreach_error_count_of_bc_file)
        plt.close()

    def hist_length_counts_reads(self):

        df_count_reads_len = pd.read_csv(self.count_reads_len_file)
        plt.bar(df_count_reads_len['length'], df_count_reads_len['count'])
        plt.xlabel('length')
        plt.ylabel('# reads')
        plt.savefig(self.len_reads_hist_output_file)
        plt.close()

        # most frequency value in reads
        most_common_length_count = max(df_count_reads_len['count'])
        most_common_length = \
            df_count_reads_len.loc[df_count_reads_len['count'] == most_common_length_count, 'length'].values[0]
        # Write general information to csv
        uts.write_list_to_csv(
            ['most_frequency_value_in_reads:', most_common_length, 'appearances:', most_common_length_count],
            self.general_information_file)

    def hist_foreach_read_count_count_bc(self, csv_output_file: Union[Path, str]) -> None:
        count_sorted_with_missing_bc = pd.read_csv(csv_output_file)
        x = count_sorted_with_missing_bc["count"]
        plt.hist(x, bins=20)
        plt.xlabel('# reads')
        plt.ylabel('bc')
        plt.savefig(self.hist_foreach_read_count_count_bc_file)
        plt.close()

    @uts.timer
    def hist_reads(self, csv_output_file: Union[Path, str], input_csv_path: Path) -> None:
        self.for_each_bc_count_reads_read_to_csv(output_file=csv_output_file, input_csv_path=input_csv_path)
        self.hist_foreach_bc_read_count(csv_output_file)
        self.hist_foreach_read_count_count_bc(csv_output_file=csv_output_file)
        self.hist_foreach_error_count_of_bc()
        # self.hist_length_counts_reads()

    @uts.timer
    def analyze_data_with_the_chosen_approach(self, algo_approach: str,
                                              const_design_pd: Union[TextFileReader, DataFrame],
                                              payload_design_pd: Union[TextFileReader, DataFrame],
                                              barcodes_design_pd: Union[TextFileReader, DataFrame], ) -> None:

        uts.fasta_to_fastq(self.input_file_or_folder)
        fastq_files = self.process_fastq_files()

        if algo_approach == 'use_reads_with_len_bigger_then_y_and_use_all_u_for_pos':
            number_of_reads, is_run_loop, reads_length_counts, blast_database_name = self.prepare_run_extract_all_pos_and_reads_results_to_csv(
                output_csv_path=self.results_good_reads_with_len_bigger_then_y)

        for fastq_file in fastq_files:
            print(f'file = {fastq_file}')
            self.input_file_or_folder = fastq_file

            # Only analyze reads with the good length, therefore all other reads are being ignored
            if algo_approach == 'use_only_good_reads_with_correct_len':
                # reads
                reads = uts.open_fastq(input_file=self.input_file_or_folder)
                # good reads with len 575
                good_reads = self.retrieve_reads_in_specific_len(reads=reads, length=self.len_reads_to_retrieve)

                # Write the good reads with len 575 to results_good_reads.csv
                self.reads_results_to_csv(reads=good_reads,
                                          const_design=const_design_pd,
                                          payload_design=payload_design_pd,
                                          barcodes_design=barcodes_design_pd,
                                          output_csv_path=self.results_good_reads_file)
            # 1. Use all reads with len bigger then some y
            # 2. Extract the universal2 -> extract the rest of the seqs
            elif algo_approach == 'use_reads_with_len_bigger_then_y_and_u2_as_start_pos':
                # reads
                reads = uts.open_fastq(input_file=self.input_file_or_folder)
                self.extract_start_position_and_reads_results_to_csv(reads=reads,
                                                                     const_design=const_design_pd,
                                                                     payload_design=payload_design_pd,
                                                                     barcodes_design=barcodes_design_pd,
                                                                     dist_option='levenshtein',
                                                                     output_csv_path=self.results_good_reads_with_len_bigger_then_y)

            elif algo_approach == 'use_reads_with_len_bigger_then_y_and_use_all_u_for_pos':
                self.extract_all_pos_and_reads_results_to_csv(payload_design=payload_design_pd,
                                                              barcodes_design=barcodes_design_pd,
                                                              dist_option='levenshtein',
                                                              output_csv_path=self.results_good_reads_with_len_bigger_then_y,
                                                              number_of_reads=number_of_reads,
                                                              is_run_loop=is_run_loop,
                                                              reads_length_counts=reads_length_counts,
                                                              blast_database_name=blast_database_name)

    def create_folders(self) -> None:
        uts.is_dir_exists(self.output_hist_folder)
        uts.is_dir_exists(self.output_folder)
        uts.is_dir_exists(self.design_folder)
        uts.is_dir_exists(self.output_heatmap_folder)
        uts.is_dir_exists(self.output_graphs_folder)
        uts.is_dir_exists(self.output_csv_folder)
        uts.is_dir_exists(self.output_line_graphs_folder)
        uts.is_dir_exists(self.output_fasta_folder)
        uts.is_dir_exists(self.blast_database_folder)

    def delete_csv_files_from_folder(self):
        for item in os.listdir(self.output_csv_folder):
            if item.endswith(".csv"):
                file_path = os.path.join(self.output_csv_folder, item)
                try:
                    os.remove(file_path)
                    print(f"Deleted {file_path}")
                except OSError as e:
                    print(f"Error deleting {file_path}: {e}")

    def save_output_into_separate_folder(self, sampling_percentage: float, iteration: int):
        output_folder = f"{self.output_all_sampling_rate_folder}output_{sampling_percentage}/iter_{iteration}"
        os.makedirs(output_folder, exist_ok=True)
        self.move_output_to_sampling_percentage_directory(new_output_dir=output_folder)
        return output_folder

    def move_output_to_sampling_percentage_directory(self, new_output_dir: str):
        src_dir = self.output_folder
        dst_dir = new_output_dir
        # Create the destination directory if it does not exist
        if not os.path.exists(dst_dir):
            try:
                os.makedirs(dst_dir)
            except PermissionError:
                print(f"Unable to create directory {dst_dir}: Permission denied")
                return

        # Change the permissions of the destination directory
        try:
            os.chmod(dst_dir, 0o777)
        except PermissionError:
            print(f"Unable to change permissions for directory {dst_dir}: Permission denied")
            return

        # Move the output files to the destination directory
        for item in os.listdir(src_dir):
            src_path = os.path.join(src_dir, item)
            dst_path = os.path.join(dst_dir, item)
            if os.path.isfile(src_path):
                try:
                    # Check if a file with the same name already exists in the destination folder
                    if os.path.exists(dst_path):
                        # Rename the source file
                        dst_path = os.path.join(dst_dir,
                                                f"{os.path.splitext(item)[0]}_{len(os.listdir(dst_dir))}{os.path.splitext(item)[1]}")
                    shutil.copy2(src_path, dst_path)
                except PermissionError:
                    print(f"Unable to copy file {src_path} to {dst_path}: Permission denied")
            elif os.path.isdir(src_path):
                try:
                    # Check if the destination path already exists and is a directory
                    if os.path.exists(dst_path) and os.path.isdir(dst_path):
                        dst_path = os.path.join(dst_path, os.path.basename(src_path))
                    shutil.copytree(src_path, dst_path)
                except PermissionError:
                    print(f"Unable to copy directory {src_path} to {dst_path}: Permission denied")
            else:
                print(f"Skipping {src_path} as it is not a file or directory")

        # Change the permissions of the destination directory back to its original value
        try:
            os.chmod(dst_dir, 0o755)
        except PermissionError:
            print(f"Unable to change permissions for directory {dst_dir}: Permission denied")

    def graph_of_all_sampling_most_common_x_to_each_bc(self):
        all_results_count = {}
        x_axis = []
        df_design_results_only_z_file = pd.read_csv(self.design_results_only_z_file, index_col=0)

        for folder_name in os.listdir(self.output_all_sampling_rate_folder):
            # if float(folder_name.split('_')[1]) not in [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]:
            #     continue
            # folder_path = os.path.join(self.output_all_sampling_rate_folder, folder_name, "csv")
            sampling_percentage_folder = os.path.join(self.output_all_sampling_rate_folder, folder_name)
            if os.path.isdir(sampling_percentage_folder):
                sampling_percentage = float(folder_name.split("_")[-1])
                all_results_count[sampling_percentage] = {}
                if sampling_percentage not in self.sampling_rate_array:
                    continue

                for iteration_folder in os.listdir(sampling_percentage_folder):
                    # sampling_percentage_folder = "output_all_sampling_rate/output_100"
                    if iteration_folder.startswith("iter_"):
                        all_results_count[sampling_percentage][iteration_folder] = []
                        iter_path = os.path.join(sampling_percentage_folder, iteration_folder)
                        csv_file = os.path.join(iter_path, "csv/foreach_bc_payload_count.csv")

                        if os.path.exists(csv_file):
                            data = pd.read_csv(csv_file, header=None)

                # csv_file = os.path.join(folder_path, "foreach_bc_payload_count.csv")

                            for bc_i, item in data.items():
                                x_count_to_each_bc = []
                                if bc_i == 0 or bc_i not in self.bc_list:
                                    continue
                                # Parse the dictionary from the cell
                                dict_data = ast.literal_eval(item[1])

                                for cycle_i, (cycle_name, payload_count_dict) in enumerate(dict_data.items()):
                                    # Create a Counter from the dictionary
                                    counter = collections.Counter(payload_count_dict)
                                    # Remove if the payload is 0
                                    del counter[0]
                                    # Get the 4 most common items
                                    most_common_payloads = counter.most_common(self.subset_size)
                                    # Filter the payloads that have 0 appearances.
                                    filtered_common = [item for item, count in most_common_payloads if count > 0]
                                    x_count_to_each_bc.append(self.compare_amount_of_x_to_design(set_most_common_x_to_each_bc=set(filtered_common),
                                                                       z=df_design_results_only_z_file[str(bc_i)][cycle_i]))

                                all_results_count[sampling_percentage][iteration_folder].append(x_count_to_each_bc)

        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='errorbar', x_axis_type='sampling_rate')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='boxplot', x_axis_type='sampling_rate')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='errorbar', x_axis_type='bc_avg_coverage')
        self.create_graph_avg_x_to_each_sampling(all_counts=all_results_count, graph_type='boxplot', x_axis_type='bc_avg_coverage')

    def compare_amount_of_x_to_design(self, set_most_common_x_to_each_bc, z):
        z_to_x_tuple = self.z_to_k_mer_representative[z]
        set_of_x_vals = set([int(x[1:]) for x in z_to_x_tuple])

        # Get common elements
        common_elements = set_most_common_x_to_each_bc.intersection(set_of_x_vals)

        # Get count of common elements
        count = len(common_elements)

        return count

    def create_graph_avg_x_to_each_sampling(self, all_counts, graph_type, x_axis_type):
        # x_axis = []
        data = {}
        data_avg = {}

        for sampling_percentage, iterations_data in all_counts.items():
            data[sampling_percentage] = []

            # Determine the number of elements in each iteration (assuming equal for all iterations)
            num_elements = len(
                iterations_data[next(iter(iterations_data))])  # Get the length of one of the iteration lists

            # Iterate over the number of elements
            for i in range(num_elements):
                temp_data = []

                # Iterate over each iteration
                for iteration in iterations_data:
                    # Append the i-th element of each iteration
                    # Check if the i-th element exists in this iteration
                    if i < len(iterations_data[iteration]):
                        temp_data.extend(iterations_data[iteration][i])
                data[sampling_percentage].append(temp_data)  # Extract the value from the list

        for folder_name in os.listdir(self.output_all_sampling_rate_folder):
            # if float(folder_name.split('_')[1]) not in [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]:
            #     continue
            # folder_name = "output_100"
            x_axis = []
            avg_count_per_iter = []
            folder_path = os.path.join(self.output_all_sampling_rate_folder, folder_name)
            sampling_percentage = float(folder_name.split("_")[-1])
            # Add the sampling percentage and counts to the lists
            if x_axis_type == 'sampling_rate':
                x_axis.extend([sampling_percentage] * len(all_counts[sampling_percentage]))
            elif x_axis_type == 'bc_avg_coverage':
                    for iter_i, iterations in all_counts[sampling_percentage].items():
                        csv_file_name = self.count_reads_for_each_bc_file.replace("output/", "")
                        csv_file = os.path.join(folder_path,iter_i ,csv_file_name)
                        try:
                            avg_count = round(self.calculate_average_coverage(csv_file=csv_file))
                        except:
                            print(csv_file)
                        avg_count_per_iter.append([avg_count])
                    x_axis_avg = round(sum([item[0] for item in avg_count_per_iter]) / len(avg_count_per_iter))
                    x_axis.extend([x_axis_avg] * self.amount_of_bc)

                    # Group data points by sampling percentage
                    for bc in range(self.amount_of_bc):
                        data_avg_temp = []
                        for iter_i, iterations in all_counts[sampling_percentage].items():

                            if x_axis_avg not in data_avg:
                                data_avg[x_axis_avg] = []
                            data_avg_temp.extend(all_counts[sampling_percentage][iter_i][bc])
                        data_avg[x_axis_avg].append(data_avg_temp)

        if x_axis_type == 'bc_avg_coverage':
            data = data_avg
        # Calculate mean and standard deviation for each group
        means = []
        std_devs = []
        sorted_keys = sorted(data.keys())
        for key in sorted_keys:
            mean = np.mean(data[key])
            std_dev = np.std(data[key])
            means.append(mean)
            std_devs.append(std_dev)

        # Create a saturation analysis plot
        if graph_type == 'errorbar':
            plt.grid(True)
            plt.errorbar(sorted_keys, means, yerr=std_devs, capsize=5, fmt='o-', markersize=8)
        elif graph_type == 'boxplot':
            # plt.boxplot(sorted_keys, labels=sorted_keys)
            flattened_data = [[item for sublist in lst for item in sublist] for lst in
                              [data[key] for key in sorted_keys]]
            plt.boxplot(flattened_data)
            sorted_keys_xticks_str_arr = [str(round(element,2)) for element in sorted_keys]
            positions = range(1, len(sorted_keys) + 1)
            plt.xticks(positions, sorted_keys_xticks_str_arr, rotation=45)


        # plt.ylim(-0.5, 4.5)
        plt.ylabel('Number Of Observed $s_i$')
        # plt.title('Saturation Analysis of Design X and Experiment Results X')
        plt.subplots_adjust(bottom=0.2)
        if x_axis_type == 'sampling_rate':
            plt.xlabel('Sampling Percentage')
            plt.savefig(self.output_line_graphs_folder + 'sampling_rate_' + graph_type + '.svg')
        elif x_axis_type == 'bc_avg_coverage':
            plt.xlabel('Average Reads Per Strand')
            plt.savefig(self.output_line_graphs_folder + 'bc_avg_coverage_' + graph_type + '.svg')
        plt.show()
        plt.close()
        x=4


