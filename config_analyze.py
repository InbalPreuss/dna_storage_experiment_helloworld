import itertools


def build_config():
    shrink_dict_3_mer = {'AAT': 'X1',
                         'ACA': 'X2',
                         'ATG': 'X3',
                         'AGC': 'X4',
                         'TAA': 'X5',
                         'TCT': 'X6',
                         'TTC': 'X7',
                         'TGG': 'X8',
                         'GAG': 'X9',
                         'GCC': 'X10',
                         'GTT': 'X11',
                         'GGA': 'X12',
                         'CAC': 'X13',
                         'CCG': 'X14',
                         'CTA': 'X15',
                         'CGT': 'X16'}

    shrink_dict_size = len(shrink_dict_3_mer)

    subset_size = 5
    bits_per_z = 12
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    x_combinations = [set(k) for k in k_mer_representative]
    z = itertools.combinations(['Z' + str(i) for i in range(1, len(x_combinations) + 1)], 1)
    z = [i[0] for i in z]

    z = z[:2 ** bits_per_z]
    k_mer_representative = itertools.combinations(['X' + str(i) for i in range(1, shrink_dict_size + 1)], subset_size)
    k_mer_representative = list(k_mer_representative)[:2 ** bits_per_z]
    k_mer_representative_to_z = dict(zip(k_mer_representative, z))
    z_to_k_mer_representative = dict(zip(z, k_mer_representative))

    config = {
        'payload_pos': [60, 100, 140, 180],
        'amount_of_bc': 167,
        'design_len': 220,
        'payload_len': 20,
        'universal_len': 20,
        'barcode_len': 20,
        'subset_size': 5,
        'amount_of_payloads': 16,
        'k_mer_representative_to_z': k_mer_representative_to_z,
        'z_to_k_mer_representative': z_to_k_mer_representative,

        # 167 BC analysis
        'input_file': "data/output_prefix.assembled.fastq",
        'results_most_common_file': "output/csv/results_most_common.csv",
        'const_design_file': "config/design.csv",
        'barcodes_design_file': "config/barcodes_design.csv",
        'payload_design_file': "config/payload_design.csv",
        'results_good_reads_file': "output/csv/results_good_reads.csv",
        'count_reads_for_each_bc_file': "output/csv/count_reads_for_each_bc.csv",
        'missing_bcs_file': "output/csv/missing_bc.csv",
        'output_csv_folder': "output/csv/",
        'foreach_bc_payload_count_file': "output/csv/foreach_bc_payload_count.csv",
        'compare_design_to_experiment_results_output_file': "output/csv/compare_design_to_experiment_results.csv",
        'output_hist_folder': "output/graphs/hist/",
        'len_reads_hist_output_file': "output/graphs/hist/len_reads_hist.png",
        'output_graphs_folder': 'output/graphs/',
        'output_heatmap_folder': 'output/graphs/heatmap/',
        'heatmap_foreach_bc_and_x_count_with_most_common_file':
            "output/graphs/heatmap/heatmap_foreach_bc_and_x_count_with_most_common.png",
        'hist_per_bc_file': "output/graphs/hist/hist_per_bc",
        'hist_foreach_bc_read_count_file': "output/graphs/hist/hist_foreach_bc_read_count",
        'hist_foreach_read_count_count_bc_file': "output/graphs/hist/hist_foreach_read_count_count_bc",
        'hist_foreach_error_count_of_bc_file': "output/graphs/hist/hist_foreach_error_count_of_bc",
        'design_simulation_file': 'data/simulation_data.3.encoder_results_file.dna'
    }

    return config