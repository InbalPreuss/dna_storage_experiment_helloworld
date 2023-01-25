import ast
import heapq
from pathlib import Path
from typing import Union, Dict, List, Tuple
from operator import itemgetter
from Levenshtein import distance as lev

import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from matplotlib import pyplot as plt, patches
from pandas import DataFrame
from pandas.io.parsers import TextFileReader
from tqdm import tqdm

import utilities.utilities as utilities


def compare_with_errors(read, seq, max_dist=6):
    return sum(r != s for r, s in zip(read, seq)) <= max_dist


def most_frequency_value_in_reads(lens: Dict) -> None:
    length_counts = {}
    for len in lens:
        if len in length_counts:
            length_counts[len] += 1
        else:
            length_counts[len] = 1

    sorted_lengths = sorted(length_counts.items(), key=itemgetter(1), reverse=True)
    print(length_counts)
    most_common_length = max(length_counts, key=length_counts.get)
    print(
        f'most_frequency_value_in_reads: {most_common_length} with {length_counts[most_common_length]} apperances')


class AnalyzeFastqData:
    def __init__(self, input_file: Union[Path, str],
                 const_design_file: Union[Path, str],
                 payload_design_file: Union[Path, str],
                 barcodes_design_file: Union[Path, str],
                 len_reads_hist_output_file: Union[Path, str],
                 results_good_reads_file: Union[Path, str],
                 results_good_reads_with_len_bigger_then_y: Union[Path, str],
                 results_most_common_file: Union[Path, str],
                 design_simulation_file: Union[Path, str],
                 compare_design_to_experiment_results_output_file: Union[Path, str],
                 foreach_bc_payload_count_file: Union[Path, str],
                 heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str],
                 count_reads_for_each_bc_file: Union[Path, str],
                 missing_bcs_file: Union[Path, str],
                 output_hist_folder: Union[Path, str],
                 output_folder: Union[Path, str],
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
                 sampling_rate_from_good_reads_graph: Union[Path, str],
                 output_line_graphs_folder: Union[Path, str],
                 cycles_array: List,
                 bc_cycles_array: List,
                 universal_len: int,
                 payload_len: int,
                 five_prime_len: int,
                 three_prime_len: int,
                 th_minimum_len_reads_to_analyse: int
                 ):
        self.input_file = input_file
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
        self.sampling_rate_from_good_reads_graph = sampling_rate_from_good_reads_graph
        self.output_line_graphs_folder = output_line_graphs_folder
        self.cycles_array = cycles_array
        self.bc_cycles_array = bc_cycles_array
        self.universal_len = universal_len
        self.payload_len = payload_len
        self.three_prime_len = three_prime_len
        self.five_prime_len = five_prime_len
        self.th_minimum_len_reads_to_analyse = th_minimum_len_reads_to_analyse

    # Verify universal

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
        pos = self.universal_len + self.five_prime_len + oligo_start_pos
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
        for s_idx, s in const_design.iterrows():
            pos = int(s['Pos']) + oligo_start_pos
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

    def open_fastq(self) -> List[str]:
        with open(self.input_file, 'r') as inf:
            reads = list(SeqIO.parse(inf, 'fastq'))

        return reads

    def find_u2_location_in_read(self, read: SeqIO, u2: str) -> int:
        # TODO: what should be the maximum dist to find the universal?
        max_distance = 6
        best_location = 0
        best_distance = float('inf')
        from Levenshtein import distance
        for i in range(len(read) - len(u2) + 1):
            substring = read[i:i + len(u2)]
            distance = lev(substring, u2)
            if distance <= max_distance:
                if distance < best_distance:
                    best_distance = distance
                    best_location = i
        if best_location != 0:
            return best_location - self.barcode_len
        return best_location

    def extract_start_position_and_reads_results_to_csv(self, reads: List[str],
                                                        const_design: pd.DataFrame,
                                                        payload_design: pd.DataFrame,
                                                        barcodes_design: pd.DataFrame,
                                                        dist_option: str,
                                                        output_csv_path: Path) -> None:
        reads = self.retrieve_reads_in_specific_len_at_least(reads=reads, length=self.th_minimum_len_reads_to_analyse)
        res = list()
        failed = 0
        read_idx = 0
        none = SeqIO.SeqRecord("None")
        u2 = const_design.loc[const_design.index == "Universal2", "Seq"]['Universal2']
        with open(output_csv_path, "ab") as f:
            cols_names = [self.bc_cycles_array]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        for read_idx, read in enumerate(reads):
            if read_idx % 1000 == 999:
                print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
                with open(output_csv_path, "ab") as f:
                    np.savetxt(f, res, fmt='%i', delimiter=",")
                res = list()
            u2_location = self.find_u2_location_in_read(read=read, u2=u2)
            read = self.verify_const_universal_and_reverse_complement(read=read,
                                                                      const_design=const_design,
                                                                      oligo_start_pos=u2_location,
                                                                      dist_option=dist_option)
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

        print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
        with open(output_csv_path, "ab") as f:
            np.savetxt(f, res, fmt='%i', delimiter=",")

    def reads_len_hist(self, reads: List[str]) -> None:
        len_reads = len(reads)
        print(f'len_reads: {len_reads}')
        lens = [len(r) for r in reads]
        plt.hist(lens, bins=2000)
        plt.xlabel('length')
        plt.ylabel('# reads')
        plt.savefig(self.len_reads_hist_output_file)
        plt.close()

        # most frequency value in reads
        most_frequency_value_in_reads(lens=lens)

    def upload_design(self) -> Tuple[Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame],
                                     Union[TextFileReader, DataFrame]]:
        const_design_pd = pd.read_csv(self.const_design_file, index_col=0, header=0)
        payload_design_pd = pd.read_csv(self.payload_design_file, index_col=0, header=0)
        barcodes_design_pd = pd.read_csv(self.barcodes_design_file, index_col=0, header=0)

        return const_design_pd, payload_design_pd, barcodes_design_pd

    def retrieve_reads_in_specific_len(self, reads: List[str], length: int) -> List[str]:
        reads = [r for r in reads if len(r) == length]
        print(f'{len(reads)} reads in len {len}')
        return reads

    def retrieve_reads_in_specific_len_at_least(self, reads: List[str], length: int) -> List[str]:
        reads = [r for r in reads if len(r) >= length]
        print(f'{len(reads)} reads in len {len}')
        return reads

    def reads_results_to_csv(self, reads: List[str],
                             const_design: pd.DataFrame,
                             payload_design: pd.DataFrame,
                             barcodes_design: pd.DataFrame) -> None:
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
            u2_location = self.find_u2_location_in_read(read=read, u2=u2)
            read = self.verify_const_universal_and_reverse_complement(read=read, const_design=const_design)
            if read.seq == none.seq:
                failed += 1
                continue

            res.append(self.identify_oligo(read=read, payload_design=payload_design, barcodes_design=barcodes_design))
            if res[-1].__contains__(0):
                failed += 1

        print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
        with open(self.output_csv_path, "ab") as f:
            np.savetxt(f, res, fmt='%i', delimiter=",")

    def hist_per_bc(self, dict_bc_payload, bc, payload):
        plt.bar(dict_bc_payload.keys(), dict_bc_payload.values())
        amount_of_payloads = self.amount_of_payloads + 1
        plt.xticks(range(amount_of_payloads))
        plt.xlim(0.5, amount_of_payloads)
        plt.title('bc=' + str(bc) + 'payload=' + payload)
        utilities.is_dir_exists(self.hist_per_bc_file)
        plt.savefig(self.hist_per_bc_file + '/bc=' + str(bc) + '_payload=' + payload + '_hist.png')
        plt.close()

    def analyze_results_good_reads(self, input_csv_path: Path) -> Dict:
        # read your csv to a dataframe

        df = pd.read_csv(input_csv_path)
        dict_bc = {}

        amount_of_bc = self.amount_of_bc + 1
        for i in range(amount_of_bc):
            dict_bc_i = {}
            for payload in self.cycles_array:
                dict_p = {
                    payload: {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0,
                              14: 0,
                              15: 0, 16: 0}}
                dict_bc_i.update(dict_p)
            dict_bc[i] = dict_bc_i
        df_bc = df.sort_values(self.bc_cycles_array, ascending=True, key=np.sin)
        df_bc.fillna(0)
        for row_idx, row in tqdm(df_bc.iterrows()):
            for location_payload, payload in row[1:].items():
                try:
                    dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

                except:
                    print(f"An exception occurred {row_idx}, row={row}")
                    continue

                dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

        utilities.write_dict_to_csv(dict_bc, self.foreach_bc_payload_count_file)

        return dict_bc

    def most_common_for_each_bc(self, dict_bc: Dict) -> None:
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

        utilities.write_dict_to_csv(dict_bc_most_common, self.results_most_common_file)

    def convert_most_common_to_letters_in_new_alphabet(self) -> Dict:
        df = pd.read_csv(self.results_most_common_file, index_col=0)
        dict_most_common = df.to_dict("list")
        result_payload = {}
        print(dict_most_common)

        for bc_i in range(1, (self.amount_of_bc + 1)):
            result_payload[bc_i] = []
            for payload in self.cycles_array:
                bc_i_str = str(bc_i)
                d = ast.literal_eval(dict_most_common[bc_i_str][0])
                d = d[payload]
                dict_convert_to_x = []
                for i in d:
                    dict_convert_to_x.append('X' + str(i))
                dict_convert_to_x = tuple(utilities.sorted_human(dict_convert_to_x))

                try:
                    z = self.k_mer_representative_to_z[dict_convert_to_x]
                except KeyError:
                    z = 'Z1'  # The X tuple is out of range
                result_payload[bc_i].append({z: dict_convert_to_x})

        return result_payload

    def compare_most_common_to_design(self, result_payload: Dict) -> None:
        with open(self.compare_design_to_experiment_results_output_file, "ab") as f:
            cols_names = [
                ['bc', 'design_z', 'design_x', 'experiment_results_z', 'experiment_results_x', 'is_retrieved_all_z',
                 'count_all_z_mm']]
            np.savetxt(f, cols_names, fmt='%s', delimiter=",")

        with open(self.design_simulation_file, 'r') as inf:
            for bc_i_idx, bc_i_design in enumerate(inf):
                count_all_z_mm = 0

                array_bc_i_design = bc_i_design.split("\n")
                array_bc_i_design = array_bc_i_design[0].split(",")
                for z_i in zip(array_bc_i_design[1:self.subset_size], result_payload[bc_i_idx + 1]):
                    if z_i[0] != list(z_i[1])[0]:
                        count_all_z_mm += 1
                with open(self.compare_design_to_experiment_results_output_file, "ab") as f:
                    bc = bc_i_idx + 1
                    design_z = " ".join(array_bc_i_design[1:self.subset_size])
                    design_x = str(
                        {z: self.z_to_k_mer_representative[z] for z in array_bc_i_design[1:self.subset_size]}).replace(
                        ",",
                        " ")
                    experiment_results_z = " ".join([list(z.keys())[0] for z in result_payload[bc_i_idx + 1]])
                    experiment_results_x = str(result_payload[bc_i_idx + 1]).replace(",", " ")
                    is_retrieved_all_z = (count_all_z_mm == 0)
                    count_all_z_mm = count_all_z_mm
                    cols = [[bc, design_z, design_x, experiment_results_z, experiment_results_x, is_retrieved_all_z,
                             count_all_z_mm]]
                    np.savetxt(f, cols, fmt='%s', delimiter=",")

    def create_heatmap_with_rectangles_on_most_common(self, dict_foreach_bc_and_x_count_all_cycles_matrix: Dict,
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

        plt.savefig(self.heatmap_foreach_bc_and_x_count_with_most_common_file, dpi=400)
        plt.close()

    def heatmap_foreach_bc_and_x_count_with_most_common(self) -> None:
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
                                                           dict_foreach_bc_and_x_count_all_cycles_matrix, ax=ax)

    def find_most_common(self, input_csv_path: Path) -> None:
        dict_bc = self.analyze_results_good_reads(input_csv_path=input_csv_path)
        self.most_common_for_each_bc(dict_bc=dict_bc)
        result_payload = self.convert_most_common_to_letters_in_new_alphabet()
        self.compare_most_common_to_design(result_payload=result_payload)

        self.heatmap_foreach_bc_and_x_count_with_most_common()

    def missing_bc_to_csv(self, dict_append_missing_bc):
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        ser_append_missing_bc.to_csv(self.missing_bcs_file, mode='a', header=False)

    def create_sampling_rate_from_good_reads_graph(self, input_csv_path:Path) -> None:
        """
        This function takes in a list of integers and creates a graph with the x-axis being the sampling rate
        and the y-axis being the count of different values in arr[0].
        The sampling rate will be in increments of 10% starting from 0% until 100%.
        """

        df_good_reads = pd.read_csv(self.input_csv_path)

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
        idx = 0
        for bc, count in count_sorted.iteritems():
            if idx != bc:
                for idx in range(idx, bc):
                    dict_append_missing_bc[idx] = 0
                    idx += 1
            idx += 1

        self.missing_bc_to_csv(dict_append_missing_bc)
        ser_append_missing_bc = pd.Series(dict_append_missing_bc)
        count_sorted = count_sorted.append(ser_append_missing_bc)
        count_sorted_with_missing_bc = count_sorted.sort_index()
        count_sorted_with_missing_bc.to_csv(output_file, mode='a', header=False)

        return count_sorted_with_missing_bc

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
        plt.hist(x, bins=5)
        plt.xlabel('errors in bc')
        plt.ylabel('count')
        plt.savefig(self.hist_foreach_error_count_of_bc_file)
        plt.close()

    def hist_foreach_read_count_count_bc(self, csv_output_file: Union[Path, str]) -> None:
        count_sorted_with_missing_bc = pd.read_csv(csv_output_file)
        x = count_sorted_with_missing_bc["count"]
        plt.hist(x, bins=20)
        plt.xlabel('# reads')
        plt.ylabel('bc')
        plt.savefig(self.hist_foreach_read_count_count_bc_file)
        plt.close()

    def for_each_bc_count_reads_read(self, csv_output_file: Union[Path, str], input_csv_path: Path) -> None:
        self.for_each_bc_count_reads_read_to_csv(output_file=csv_output_file, input_csv_path=input_csv_path)
        self.hist_foreach_bc_read_count(csv_output_file)
        self.hist_foreach_read_count_count_bc(csv_output_file=csv_output_file)
        self.hist_foreach_error_count_of_bc()

    def create_folders(self) -> None:
        utilities.is_dir_exists(self.output_hist_folder)
        utilities.is_dir_exists(self.output_folder)
        utilities.is_dir_exists(self.output_heatmap_folder)
        utilities.is_dir_exists(self.output_graphs_folder)
        utilities.is_dir_exists(self.output_csv_folder)
        utilities.is_dir_exists(self.output_line_graphs_folder)

    def run(self):
        self.create_folders()

        # upload design
        const_design_pd, payload_design_pd, barcodes_design_pd = self.upload_design()

        # reads
        reads = self.open_fastq()

        # # reads len showed in histogram
        # self.reads_len_hist(reads=reads)

        # extract the universal2 -> extract the rest of the seqs
        self.extract_start_position_and_reads_results_to_csv(reads=reads,
                                                             const_design=const_design_pd,
                                                             payload_design=payload_design_pd,
                                                             barcodes_design=barcodes_design_pd,
                                                             dist_option='levenshtein',
                                                             output_csv_path=self.results_good_reads_with_len_bigger_then_y)

        input_csv_path = self.results_good_reads_with_len_bigger_then_y
        # # good reads with len 575
        # good_reads = self.retrieve_reads_in_specific_len(reads=reads, length=self.len_reads_to_retrieve)
        #
        # # Write the good reads with len 575 to results_good_reads.csv
        # self.reads_results_to_csv(reads=good_reads,
        #                           const_design=const_design_pd,
        #                           payload_design=payload_design_pd,
        #                           barcodes_design=barcodes_design_pd,
        #                           output_csv_path=results_good_reads_file)

        # Find most common for each bc and for every cycle in that bc in results of good reads
        self.find_most_common(input_csv_path=input_csv_path)

        # For each bc count amount of reads sequenced
        self.for_each_bc_count_reads_read(csv_output_file=self.count_reads_for_each_bc_file, input_csv_path=input_csv_path)

        # Create graph with sampling rate
        self.create_sampling_rate_from_good_reads_graph(input_csv_path=input_csv_path)
