from pathlib import Path
from typing import Union, Dict, List
from analyze_fastq import AnalyzeFastqData
import utilities.utilities as uts


class ConvertDesignToZ(AnalyzeFastqData):
    def find_most_common_convert_to_z(self, input_csv_path: Union[Path, str],
                                      foreach_bc_payload_count_file_dict_to_csv: Union[Path, str],
                                      most_common_dict_to_csv_path: Union[Path, str]) -> None:
        dict_bc = self.analyze_results_good_reads(input_csv_path=input_csv_path,
                                                  dict_to_csv=foreach_bc_payload_count_file_dict_to_csv)
        self.most_common_for_each_bc(dict_bc=dict_bc, dict_to_csv_path=most_common_dict_to_csv_path)
        result_payload, result_only_z = self.convert_most_common_to_letters_in_new_alphabet(
            results_most_common_file=most_common_dict_to_csv_path)
        uts.dict_to_csv(dict=result_payload, file_name=self.design_results_as_z_file)
        uts.dict_to_csv(dict=result_only_z, file_name=self.design_results_only_z_file)

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
                 three_prime_len: int, th_minimum_len_reads_to_analyse: int, design_folder: Union[Path, str], uni_start_pos: List,
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
        super().__init__(
            input_file,
            const_design_file,
            payload_design_file,
            barcodes_design_file,
            len_reads_hist_output_file,
            results_good_reads_file,
            results_good_reads_with_len_bigger_then_y,
            results_most_common_file,
            design_simulation_file,
            design_results_only_z_file,
            compare_design_to_experiment_results_output_file,
            foreach_bc_payload_count_file,
            heatmap_foreach_bc_and_x_count_with_most_common_file,
            count_reads_for_each_bc_file,
            missing_bcs_file,
            output_hist_folder,
            output_folder,
            output_graphs_folder,
            output_csv_folder,
            output_heatmap_folder,
            hist_per_bc_file,
            hist_foreach_bc_read_count_file,
            hist_foreach_error_count_of_bc_file,
            hist_foreach_read_count_count_bc_file,
            len_reads_to_retrieve,
            amount_of_bc,
            barcode_len,
            subset_size,
            amount_of_payloads,
            z_to_k_mer_representative,
            k_mer_representative_to_z,
            payload_pos,
            sampling_rate_from_good_reads_graph,
            output_line_graphs_folder,
            cycles_array,
            bc_cycles_array,
            universal_len,
            payload_len,
            five_prime_len,
            three_prime_len,
            th_minimum_len_reads_to_analyse,
            design_folder,
            uni_start_pos,
            general_information_file,
            count_reads_len_file,
            algo_approach,
            output_fasta_folder,
            universals_fasta_format_file,
            reads_chunk_to_fasta_format_file,
            amount_of_universals,
            blast_database_folder,
            blast_db,
            is_sampling_rate,
            sampling_rate_array,
            output_all_sampling_rate_folder,
            bc_list,
            amount_of_cycles

        )

        self.design_most_common_dict_to_csv_path = design_most_common_dict_to_csv_path
        self.design_foreach_bc_payload_count_file_dict_to_csv = design_foreach_bc_payload_count_file_dict_to_csv
        self.design_before_conversion_file = design_before_conversion_file
        self.design_results_as_z_file = design_results_as_z_file
        self.design_results_only_z_file = design_results_only_z_file
        self.design_results_as_x_file = design_results_as_x_file
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

    def run(self):
        # upload design
        const_design_pd, payload_design_pd, barcodes_design_pd = AnalyzeFastqData.upload_design(self=self)

        AnalyzeFastqData.create_folders(self=self)

        # reads
        reads = uts.open_fasta(input_file=self.design_before_conversion_file)

        AnalyzeFastqData.extract_start_position_and_reads_results_to_csv(self=self,
                                                                         reads=reads,
                                                                         const_design=const_design_pd,
                                                                         payload_design=payload_design_pd,
                                                                         barcodes_design=barcodes_design_pd,
                                                                         dist_option='levenshtein',
                                                                         output_csv_path=self.design_results_as_x_file)
        # number_of_reads, is_run_loop, reads_length_counts, blast_database_name = AnalyzeFastqData.prepare_run_extract_all_pos_and_reads_results_to_csv(self=self,
        #         output_csv_path=self.design_before_conversion_file)
        # AnalyzeFastqData.analyze_data_with_the_chosen_approach(self=self,
        #                                                        algo_approach=self.algo_approach,
        #                                                        const_design_pd=const_design_pd,
        #                                                        payload_design_pd=payload_design_pd,
        #                                                        barcodes_design_pd=barcodes_design_pd)
        self.find_most_common_convert_to_z(input_csv_path=self.design_results_as_x_file,
                                           foreach_bc_payload_count_file_dict_to_csv=self.design_foreach_bc_payload_count_file_dict_to_csv,
                                           most_common_dict_to_csv_path=self.design_most_common_dict_to_csv_path)
