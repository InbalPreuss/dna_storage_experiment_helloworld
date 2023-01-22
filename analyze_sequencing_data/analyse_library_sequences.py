from pathlib import Path
from typing import Union, List, Dict

import numpy as np
import pandas as pd
from Levenshtein import distance as lev
from matplotlib import pyplot as plt

import utilities.utilities as uts


def matrix_to_csv(dist_matrix: List[List[int]], df: pd.DataFrame, csv_path: str) -> None:
    df_dist_matrix = pd.DataFrame(dist_matrix)
    idx = 0
    new_col = df.Seq.values
    df_dist_matrix.insert(loc=idx, column='BCIdx', value=new_col)
    df_name_and_seq = df.ID.astype(str) + '_' + df.Seq
    newfile = np.savetxt(csv_path, df_dist_matrix.to_numpy(), delimiter=',',
                         header='GuideName,' + ','.join(df_name_and_seq), fmt="%s")


class AnalyzeLibrarySeqs:
    def __init__(self, barcodes_design_file: Union[Path, str],
                 levenshtein_per_dist_path: str,
                 hamming_per_dist_path: str,
                 heatmap_hamming_dist_path: str,
                 heatmap_levenshtein_dist_path: str,
                 amount_of_bc: int,
                 hist_count_seq_per_dist_hamming: str,
                 hist_count_seq_per_dist_levenshtein: str,
                 ):
        self.barcodes_design_file = barcodes_design_file
        self.levenshtein_per_dist_path = levenshtein_per_dist_path
        self.hamming_per_dist_path = hamming_per_dist_path
        self.amount_of_bc = amount_of_bc
        self.heatmap_levenshtein_dist_path = heatmap_levenshtein_dist_path
        self.heatmap_hamming_dist_path = heatmap_hamming_dist_path
        self.hist_count_seq_per_dist_hamming = hist_count_seq_per_dist_hamming
        self.hist_count_seq_per_dist_levenshtein = hist_count_seq_per_dist_levenshtein

    def calculate_dist(self, df: pd.DataFrame, dist_calculation_name: str):
        dict_dist = {}
        dict_dist_count = {}
        for dist in range(0, self.amount_of_bc):
            dict_dist[dist] = {}
            dict_dist_count[dist] = 0

        rows, cols = (df.shape[0], df.shape[0])
        hamming_distance_seq1_seq2_matrix_temp = [[0] * cols] * rows
        distance_seq1_seq2_matrix = np.copy(hamming_distance_seq1_seq2_matrix_temp)

        for seq1_idx, seq1 in enumerate(df.Seq):
            for seq2_idx, seq2 in enumerate(df.Seq):
                if seq1_idx == seq2_idx:
                    continue
                if 'hamming' == dist_calculation_name:
                    distance_seq1_seq2 = uts.hamming_distance(seq1, seq2)
                elif 'levenshtein' == dist_calculation_name:
                    distance_seq1_seq2 = lev(seq1, seq2)
                distance_seq1_seq2_matrix[seq1_idx][seq2_idx] = distance_seq1_seq2

                seq1_name = seq1 + "," + str(df.ID[seq1_idx])
                seq2_name = seq2 + "," + str(df.ID[seq2_idx])

                if seq1 in dict_dist[distance_seq1_seq2]:
                    dict_dist[distance_seq1_seq2][seq1_name].add(seq2_name)
                else:
                    dict_dist[distance_seq1_seq2][seq1_name] = {seq2_name}

                dict_dist_count[distance_seq1_seq2] = dict_dist_count[distance_seq1_seq2] + 1

        dict_dist_count_each_pair = {k: int(v / 2) for k, v in dict_dist_count.items()}
        print(f'dict_dist_count_each_pair: {dict_dist_count_each_pair}')
        return distance_seq1_seq2_matrix, dict_dist, dict_dist_count_each_pair

    def matrix_to_heatmap(self, dist_calculation_name: str, heatmap_dist_path: str, csv_path: str) -> None:
        df = pd.read_csv(csv_path, index_col=0)
        plt.figure(figsize=(20, 15))
        plt.xlabel('BCs', fontsize=30)
        plt.ylabel('BCs', fontsize=30)
        plt.suptitle('BCs - ' + dist_calculation_name + ' distance between each pair', fontsize=40)
        plt.imshow(df, cmap='hot', interpolation='nearest')
        cbar = plt.colorbar()
        cbar.set_label(label="red - low dist. yellow - high dist.", size=20)
        plt.savefig(heatmap_dist_path)
        plt.close()

    def hist_per_dist(self, dict_dist_count: Dict, dist_calculation_name: str,
                      hist_count_seq_per_dist_path: str) -> None:

        plt.bar(dict_dist_count.keys(), dict_dist_count.values())
        plt.xlabel('distance')
        plt.ylabel('amount of sequences')
        plt.title(dist_calculation_name)
        plt.suptitle('Distance between each pair BC seq. ' + str(self.amount_of_bc) + ' seqs')
        plt.savefig(hist_count_seq_per_dist_path)
        plt.close()

    def run(self):
        df_bc = uts.open_csv_file_to_df(self.barcodes_design_file)

        dist_info = [
            {'dist_calculation_name': 'hamming',
             'csv_path': self.hamming_per_dist_path,
             'heatmap_dist_path': self.heatmap_hamming_dist_path,
             'hist_count_seq_per_dist_path': self.hist_count_seq_per_dist_hamming},
            {'dist_calculation_name': 'levenshtein',
             'csv_path': self.levenshtein_per_dist_path,
             'heatmap_dist_path': self.heatmap_levenshtein_dist_path,
             'hist_count_seq_per_dist_path': self.hist_count_seq_per_dist_levenshtein}]
        for dict_dist_info in dist_info:
            # dist between a pair of BC
            dist_matrix, dict_per_dist, dict_dist_count = self.calculate_dist(df=df_bc,
                                                                              dist_calculation_name=dict_dist_info[
                                                                                  'dist_calculation_name'])
            matrix_to_csv(dist_matrix=dist_matrix, df=df_bc, csv_path=dict_dist_info['csv_path'])
            self.matrix_to_heatmap(dist_calculation_name=dict_dist_info['dist_calculation_name'],
                                   heatmap_dist_path=dict_dist_info['heatmap_dist_path'],
                                   csv_path=dict_dist_info['csv_path'])
            uts.dict_to_csv(dict_dist=dict_per_dist, file_name=dict_dist_info['csv_path'])
            self.hist_per_dist(dict_dist_count=dict_dist_count,
                               dist_calculation_name=dict_dist_info['dist_calculation_name'],
                               hist_count_seq_per_dist_path=dict_dist_info['hist_count_seq_per_dist_path']
                               )



if __name__ == '__main__':
    pass
