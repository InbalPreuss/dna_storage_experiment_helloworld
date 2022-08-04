config = {
    'payload_pos': [60, 100, 140, 180],
    'amount_of_bc': 167,
    'design_len': 220,
    'payload_len': 20,
    'universal_len': 20,
    'barcode_len': 20,
    'subset_size': 5,
    'amount_of_payloads': 16,

    # # 2 BC analysis
    # 'input_file': "D:\\Downloads\\merged_robot_two_BCs\\merged_robot_two_BCs.fastq",
    # 'const_design_file': "D:\\Downloads\\merged_robot_two_BCs\\design.csv",
    # 'barcodes_design_file': "D:\\Downloads\\merged_robot_two_BCs\\barcodes_design.csv",
    # 'payload_design_file': "D:\\Downloads\\merged_robot_two_BCs\\payload_design.csv",
    # 'data_csv': "D:\\Downloads\\merged_robot_two_BCs\\data\\csv\\",
    # 'data_hist': "D:\\Downloads\\merged_robot_two_BCs\\data\\histograms\\",

    # 167 BC analysis
    'input_file': "data/output_full_167bc_.assembled.fastq",
    'const_design_file': "config/design.csv",
    'barcodes_design_file': "config/barcodes_design.csv",
    'payload_design_file': "config/payload_design.csv",
    'data_csv': "output/csv/",
    'data_hist': "output/hist/"
}


def get_config():
    return config
