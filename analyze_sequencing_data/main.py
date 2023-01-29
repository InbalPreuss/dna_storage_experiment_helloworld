from analyze_sequencing_data.config import config
from analyze_sequencing_data.analyze_fastq import AnalyzeFastqData
from analyze_sequencing_data.analyse_library_sequences import AnalyzeLibrarySeqs
from analyze_sequencing_data.convert_design_file_to_z import ConvertDesignToZ

if __name__ == '__main__':
    # Analyze fastq results from sequencing
    config = config.build_config()

    # convert_design_file_to_z = ConvertDesignToZ(design_before_conversion_file=config['design_before_conversion_file'],
    #                                             design_results_as_z_file=config['design_results_as_z_file'],
    #                                             design_results_only_z_file=config['design_results_only_z_file'],
    #                                             design_results_as_x_file=config['design_results_as_x_file'],
    #                                             design_foreach_bc_payload_count_file_dict_to_csv=config['design_foreach_bc_payload_count_file_dict_to_csv'],
    #                                             design_most_common_dict_to_csv_path=config['design_most_common_dict_to_csv_path'],
    #                                             input_file=config['input_file'],
    #                                             const_design_file=config['const_design_file'],
    #                                             payload_design_file=config['payload_design_file'],
    #                                             barcodes_design_file=config['barcodes_design_file'],
    #                                             len_reads_hist_output_file=config['len_reads_hist_output_file'],
    #                                             results_good_reads_file=config['results_good_reads_file'],
    #                                             results_good_reads_with_len_bigger_then_y=config[
    #                                                 'results_good_reads_with_len_bigger_then_y'],
    #                                             results_most_common_file=config['results_most_common_file'],
    #                                             design_simulation_file=config['design_simulation_file'],
    #                                             compare_design_to_experiment_results_output_file=config[
    #                                                 'compare_design_to_experiment_results_output_file'],
    #                                             foreach_bc_payload_count_file=config['foreach_bc_payload_count_file'],
    #                                             heatmap_foreach_bc_and_x_count_with_most_common_file=config[
    #                                                 'heatmap_foreach_bc_and_x_count_with_most_common_file'],
    #                                             count_reads_for_each_bc_file=config['count_reads_for_each_bc_file'],
    #                                             missing_bcs_file=config['missing_bcs_file'],
    #                                             output_hist_folder=config['output_hist_folder'],
    #                                             output_folder=config['output_folder'],
    #                                             output_graphs_folder=config['output_graphs_folder'],
    #                                             output_csv_folder=config['output_csv_folder'],
    #                                             output_heatmap_folder=config['output_heatmap_folder'],
    #                                             hist_per_bc_file=config['hist_per_bc_file'],
    #                                             hist_foreach_bc_read_count_file=config[
    #                                                 'hist_foreach_bc_read_count_file'],
    #                                             hist_foreach_error_count_of_bc_file=config[
    #                                                 'hist_foreach_error_count_of_bc_file'],
    #                                             hist_foreach_read_count_count_bc_file=config[
    #                                                 'hist_foreach_read_count_count_bc_file'],
    #                                             len_reads_to_retrieve=config['design_len'],
    #                                             amount_of_bc=config['amount_of_bc'],
    #                                             barcode_len=config['barcode_len'],
    #                                             subset_size=config['subset_size'],
    #                                             amount_of_payloads=config['amount_of_payloads'],
    #                                             z_to_k_mer_representative=config['z_to_k_mer_representative'],
    #                                             k_mer_representative_to_z=config['k_mer_representative_to_z'],
    #                                             payload_pos=config['payload_pos'],
    #                                             sampling_rate_from_good_reads_graph=config[
    #                                                 'sampling_rate_from_good_reads_graph'],
    #                                             output_line_graphs_folder=config['output_line_graphs_folder'],
    #                                             cycles_array=config['cycles_array'],
    #                                             bc_cycles_array=config['bc_cycles_array'],
    #                                             payload_len=config['payload_len'],
    #                                             universal_len=config['universal_len'],
    #                                             five_prime_len=config['five_prime_len'],
    #                                             three_prime_len=config['three_prime_len'],
    #                                             th_minimum_len_reads_to_analyse=config[
    #                                                 'th_minimum_len_reads_to_analyse']
    #                                             )
    # convert_design_file_to_z.run()

    # Fastq analysis
    analyze_fastq_data = AnalyzeFastqData(input_file=config['input_file'],
                                          const_design_file=config['const_design_file'],
                                          payload_design_file=config['payload_design_file'],
                                          barcodes_design_file=config['barcodes_design_file'],
                                          len_reads_hist_output_file=config['len_reads_hist_output_file'],
                                          results_good_reads_file=config['results_good_reads_file'],
                                          results_good_reads_with_len_bigger_then_y=config[
                                              'results_good_reads_with_len_bigger_then_y'],
                                          results_most_common_file=config['results_most_common_file'],
                                          design_simulation_file=config['design_results_only_z_file'],
                                          design_results_only_z_file=config['design_results_only_z_file'],
                                          compare_design_to_experiment_results_output_file=config[
                                              'compare_design_to_experiment_results_output_file'],
                                          foreach_bc_payload_count_file=config['foreach_bc_payload_count_file'],
                                          heatmap_foreach_bc_and_x_count_with_most_common_file=config[
                                              'heatmap_foreach_bc_and_x_count_with_most_common_file'],
                                          count_reads_for_each_bc_file=config['count_reads_for_each_bc_file'],
                                          missing_bcs_file=config['missing_bcs_file'],
                                          output_hist_folder=config['output_hist_folder'],
                                          output_folder=config['output_folder'],
                                          output_graphs_folder=config['output_graphs_folder'],
                                          output_csv_folder=config['output_csv_folder'],
                                          output_heatmap_folder=config['output_heatmap_folder'],
                                          hist_per_bc_file=config['hist_per_bc_file'],
                                          hist_foreach_bc_read_count_file=config['hist_foreach_bc_read_count_file'],
                                          hist_foreach_error_count_of_bc_file=config[
                                              'hist_foreach_error_count_of_bc_file'],
                                          hist_foreach_read_count_count_bc_file=config[
                                              'hist_foreach_read_count_count_bc_file'],
                                          len_reads_to_retrieve=config['design_len'],
                                          amount_of_bc=config['amount_of_bc'],
                                          barcode_len=config['barcode_len'],
                                          subset_size=config['subset_size'],
                                          amount_of_payloads=config['amount_of_payloads'],
                                          z_to_k_mer_representative=config['z_to_k_mer_representative'],
                                          k_mer_representative_to_z=config['k_mer_representative_to_z'],
                                          payload_pos=config['payload_pos'],
                                          sampling_rate_from_good_reads_graph=config[
                                              'sampling_rate_from_good_reads_graph'],
                                          output_line_graphs_folder=config['output_line_graphs_folder'],
                                          cycles_array=config['cycles_array'],
                                          bc_cycles_array=config['bc_cycles_array'],
                                          payload_len=config['payload_len'],
                                          universal_len=config['universal_len'],
                                          five_prime_len=config['five_prime_len'],
                                          three_prime_len=config['three_prime_len'],
                                          th_minimum_len_reads_to_analyse=config['th_minimum_len_reads_to_analyse']
                                          )
    analyze_fastq_data.run()

    # analyze_library_seqs = AnalyzeLibrarySeqs(barcodes_design_file=config['barcodes_design_file'],
    #                                           levenshtein_per_dist_path=config['levenshtein_per_dist_path'],
    #                                           hamming_per_dist_path=config['hamming_per_dist_path'],
    #                                           heatmap_hamming_dist_path=config['heatmap_hamming_dist_path'],
    #                                           heatmap_levenshtein_dist_path=config['heatmap_levenshtein_dist_path'],
    #                                           amount_of_bc=config['amount_of_bc'],
    #                                           hist_count_seq_per_dist_hamming=config['hist_count_seq_per_dist_hamming'],
    #                                           hist_count_seq_per_dist_levenshtein=config['hist_count_seq_per_dist_levenshtein']
    #                                           )
    # analyze_library_seqs.run()
