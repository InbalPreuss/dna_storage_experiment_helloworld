from pathlib import Path
from typing import Union, Dict, List
from analyze_fastq import AnalyzeFastqData
import utilities.utilities as uts


class ConvertDesignToZ(AnalyzeFastqData):
    def __init__(self, design_before_conversion_file: Path, input_file: Union[Path, str],
                 design_results_as_z_file: Union[Path, str],
                 design_results_only_z_file: Union[Path, str],
                 design_foreach_bc_payload_count_file_dict_to_csv: Union[Path, str],
                 design_most_common_dict_to_csv_path: Union[Path, str],
                 design_results_as_x_file: Union[Path, str],
                 const_design_file: Union[Path, str], payload_design_file: Union[Path, str],
                 barcodes_design_file: Union[Path, str], len_reads_hist_output_file: Union[Path, str],
                 results_good_reads_file: Union[Path, str], results_good_reads_with_len_bigger_then_y: Union[Path, str],
                 results_most_common_file: Union[Path, str], design_simulation_file: Union[Path, str],
                 compare_design_to_experiment_results_output_file: Union[Path, str],
                 foreach_bc_payload_count_file: Union[Path, str],
                 heatmap_foreach_bc_and_x_count_with_most_common_file: Union[Path, str],
                 count_reads_for_each_bc_file: Union[Path, str], missing_bcs_file: Union[Path, str],
                 output_hist_folder: Union[Path, str], output_folder: Union[Path, str],
                 output_graphs_folder: Union[Path, str], output_csv_folder: Union[Path, str],
                 output_heatmap_folder: Union[Path, str], hist_per_bc_file: Union[Path, str],
                 hist_foreach_bc_read_count_file: Union[Path, str],
                 hist_foreach_error_count_of_bc_file: Union[Path, str],
                 hist_foreach_read_count_count_bc_file: Union[Path, str], len_reads_to_retrieve: int, amount_of_bc: int,
                 barcode_len: int, subset_size: int, amount_of_payloads: int, z_to_k_mer_representative: Dict,
                 k_mer_representative_to_z: Dict, payload_pos: List,
                 sampling_rate_from_good_reads_graph: Union[Path, str], output_line_graphs_folder: Union[Path, str],
                 cycles_array: List, bc_cycles_array: List, universal_len: int, payload_len: int, five_prime_len: int,
                 three_prime_len: int, th_minimum_len_reads_to_analyse: int):
        super().__init__(input_file, const_design_file, payload_design_file, barcodes_design_file,
                         len_reads_hist_output_file, results_good_reads_file, results_good_reads_with_len_bigger_then_y,
                         results_most_common_file, design_simulation_file,
                         compare_design_to_experiment_results_output_file, foreach_bc_payload_count_file,
                         heatmap_foreach_bc_and_x_count_with_most_common_file, count_reads_for_each_bc_file,
                         missing_bcs_file, output_hist_folder, output_folder, output_graphs_folder, output_csv_folder,
                         output_heatmap_folder, hist_per_bc_file, hist_foreach_bc_read_count_file,
                         hist_foreach_error_count_of_bc_file, hist_foreach_read_count_count_bc_file,
                         len_reads_to_retrieve, amount_of_bc, barcode_len, subset_size, amount_of_payloads,
                         z_to_k_mer_representative, k_mer_representative_to_z, payload_pos,
                         sampling_rate_from_good_reads_graph, output_line_graphs_folder, cycles_array, bc_cycles_array,
                         universal_len, payload_len, five_prime_len, three_prime_len, th_minimum_len_reads_to_analyse)
        self.design_most_common_dict_to_csv_path = design_most_common_dict_to_csv_path
        self.design_foreach_bc_payload_count_file_dict_to_csv = design_foreach_bc_payload_count_file_dict_to_csv
        self.design_before_conversion_file = design_before_conversion_file
        self.design_results_as_z_file = design_results_as_z_file
        self.design_results_only_z_file = design_results_only_z_file
        self.design_results_as_x_file = design_results_as_x_file

    def find_most_common_convert_to_z(self, input_csv_path: Union[Path, str],
                                      foreach_bc_payload_count_file_dict_to_csv: Union[Path, str],
                                      most_common_dict_to_csv_path: Union[Path, str]) -> None:
        # dict_bc = self.analyze_results_good_reads(input_csv_path=input_csv_path,
        #                                           dict_to_csv=foreach_bc_payload_count_file_dict_to_csv)
        # self.most_common_for_each_bc(dict_bc=dict_bc, dict_to_csv_path=most_common_dict_to_csv_path)
        result_payload, result_only_z = self.convert_most_common_to_letters_in_new_alphabet(
            results_most_common_file=most_common_dict_to_csv_path)
        # uts.dict_to_csv(dict=result_payload, file_name=self.design_results_as_z_file)
        uts.dict_to_csv(dict=result_only_z, file_name=self.design_results_only_z_file)

    def run(self):
        # # upload design
        # const_design_pd, payload_design_pd, barcodes_design_pd = AnalyzeFastqData.upload_design(self=self)
        #
        # AnalyzeFastqData.create_folders(self=self)
        #
        # # reads
        # reads = uts.open_fasta(input_file=self.design_before_conversion_file)
        #
        # AnalyzeFastqData.extract_start_position_and_reads_results_to_csv(self=self,
        #                                                                  reads=reads,
        #                                                                  const_design=const_design_pd,
        #                                                                  payload_design=payload_design_pd,
        #                                                                  barcodes_design=barcodes_design_pd,
        #                                                                  dist_option='levenshtein',
        #                                                                  output_csv_path=self.design_results_as_x_file)
        self.find_most_common_convert_to_z(input_csv_path=self.design_results_as_x_file,
                                           foreach_bc_payload_count_file_dict_to_csv=self.design_foreach_bc_payload_count_file_dict_to_csv,
                                           most_common_dict_to_csv_path=self.design_most_common_dict_to_csv_path)
