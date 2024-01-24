import itertools

# Need to choose the algorithm approach you want to use to analyze the data
algo_approachs = {
    0: 'analyze_only_good_reads',
    1: 'use_reads_with_len_bigger_then_y_and_u2_as_start_pos',
    2: 'use_reads_with_len_bigger_then_y_and_use_all_u_for_pos'
}


def build_config():
    amount_of_payloads = 32
    subset_size = 32
    bits_per_z = 84

    shrink_dict_size = amount_of_payloads
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    x_combinations = [set(k) for k in k_mer_representative]
    z = itertools.combinations(['Z' + str(i) for i in range(1, len(x_combinations) + 1)], 1)
    z = [i[0] for i in z]

    z = z[:2 ** bits_per_z]
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    k_mer_representative = list(k_mer_representative)[:2 ** bits_per_z]
    k_mer_representative_to_z = dict(zip(k_mer_representative, z))
    z_to_k_mer_representative = dict(zip(z, k_mer_representative))  # Z1: {X1,X2,X3,X4}

    cycles_array = ["c1"]

    config = {
        # Need to choose the algorithm approach you want to use to analyze the data
        'algo_approach': algo_approachs[2],
        'payload_pos': [50],
        'uni_start_pos': [25],
        'amount_of_bc': 8,
        'amount_of_payloads': amount_of_payloads,
        'amount_of_universals': 1,
        'design_len': 75,
        'payload_len': 25,
        'universal_len': 25,
        'barcode_len': 24,
        'five_prime_len': 0,
        'three_prime_len': 0,
        'subset_size': subset_size,
        'th_minimum_len_reads_to_analyse': 50,
        'is_sampling_rate': True,
        'sampling_rate_array': [100, 50, 30, 20, 10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001],
        'bc_list': [1, 2,3,4,5,6,7,8],
        'k_mer_representative_to_z': k_mer_representative_to_z,
        'z_to_k_mer_representative': z_to_k_mer_representative,
        "cycles_array": cycles_array,
        "bc_cycles_array": ["bc"] + cycles_array,
        'amount_of_cycles': len(cycles_array),

        ### Config files ###
        'const_design_file': "config/design.csv",
        'barcodes_design_file': "config/barcodes_design.csv",
        'payload_design_file': "config/payload_design.csv",

        ### Data for analysis files ###
        # 'input_file_or_folder': "data/sequencing_output_fastq/all_sequencing_output.fastq",
        'input_file_or_folder': "data/sequencing_output_fastq/",
        'design_before_conversion_file': "data/design/all_design_before_parse.fa",

        ### Analysis files ###
        ## output ##
        'output_folder': "output/",
        'output_all_sampling_rate_folder': "output_all_sampling_rate/",
        ## Blast database ##
        'blast_database_folder': "output/blast_database/",
        'blast_db': "output/blast_database/blast_db",
        # FASTA
        'output_fasta_folder': "output/fasta/",
        'universals_fasta_format_file': "output/fasta/universals_fasta_format.fasta",
        'reads_chunk_to_fasta_format_file': "output/fasta/reads_chunk_to_fasta_format",
        # csv
        'output_csv_folder': "output/csv/",
        'general_information_file': "output/csv/general_information.csv",
        'count_reads_len_file': "output/csv/count_reads_len.csv",
        'results_most_common_file': "output/csv/results_most_common.csv",
        'results_good_reads_file': "output/csv/results_good_reads.csv",
        'results_good_reads_with_len_bigger_then_y': "output/csv/results_good_reads_with_len_bigger_then_y.csv",
        'count_reads_for_each_bc_file': "output/csv/count_reads_for_each_bc.csv",
        'missing_bcs_file': "output/csv/missing_bc.csv",
        'foreach_bc_payload_count_file': "output/csv/foreach_bc_payload_count.csv",
        'compare_design_to_experiment_results_output_file': "output/csv/compare_design_to_experiment_results.csv",
        # graphs
        'output_graphs_folder': 'output/graphs/',
        'output_hist_folder': "output/graphs/hist/",
        'len_reads_hist_output_file': "output/graphs/hist/hist_len_reads.png",
        'output_line_graphs_folder': 'output/graphs/line_graphs/',
        'sampling_rate_from_good_reads_graph': 'output/graphs/line_graphs/sampling_rate_from_good_reads_graph',
        'output_heatmap_folder': 'output/graphs/heatmap/',
        'heatmap_foreach_bc_and_x_count_with_most_common_file':
            "output/graphs/heatmap/heatmap_foreach_bc_and_x_count_with_most_common.png",
        'hist_per_bc_file': "output/graphs/hist/hist_per_bc",
        'hist_foreach_bc_read_count_file': "output/graphs/hist/hist_foreach_bc_read_count",
        'hist_foreach_read_count_count_bc_file': "output/graphs/hist/hist_foreach_read_count_count_bc",
        'hist_foreach_error_count_of_bc_file': "output/graphs/hist/hist_foreach_error_count_of_bc",

        ### Design parsed output files ###
        'design_folder': 'output/csv/design/',
        'design_results_as_z_file': 'output/csv/design/design_results_as_z.csv',
        'design_results_only_z_file': 'output/csv/design/design_results_only_z.csv',
        'design_results_as_x_file': 'output/csv/design/design_results_as_x.csv',
        'design_most_common_dict_to_csv_path': 'output/csv/design/design_most_common_dict_to_csv.csv',
        'design_foreach_bc_payload_count_file_dict_to_csv': 'output/csv/design/design_foreach_bc_payload_count_file_dict_to_csv.csv',

        ### motif analysis files ###
        'motif_analysis_folder': 'output/motif_analysis/',
        'motif_analysis_csv_folder': 'output/motif_analysis/csv/',
        'motif_analysis_graphs_folder': 'output/motif_analysis/graphs/',
        'motif_analysis_hist_folder': 'output/motif_analysis/graphs/hist/',
        'motif_analysis_heatmap_folder': 'output/motif_analysis/graphs/heatmap/',
        'levenshtein_per_dist_path': 'output/motif_analysis/csv/levenshtein_per_distance_',
        'hamming_per_dist_path': 'output/motif_analysis/csv/hamming_per_distance_',
        'heatmap_levenshtein_dist_path': 'output/motif_analysis/graphs/heatmap/heatmap_levenshtein_dist_',
        'heatmap_hamming_dist_path': 'output/motif_analysis/graphs/heatmap/heatmap_hamming_dist_',
        'hist_count_seq_per_dist_hamming': 'output/motif_analysis/graphs/hist/hist_count_seq_per_dist_hamming_',
        'hist_count_seq_per_dist_levenshtein': 'output/motif_analysis/graphs/hist/hist_count_seq_per_dist_levenshtein_',

    }

    return config
