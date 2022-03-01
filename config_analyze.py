config = {
    'payload_pos': [60, 100, 140, 180],
    'design_len': 220,
    'payload_len': 20,
    'universal_len': 20,
    'barcode_len': 20,
    'input_file': "D:\\Downloads\\merged_robot_two_BCs\\merged_robot_two_BCs.fastq",
    'const_design_file': "D:\\Downloads\\merged_robot_two_BCs\\design.csv",
    'barcodes_design_file': "D:\\Downloads\\merged_robot_two_BCs\\barcodes_design.csv",
    'payload_design_file': "D:\\Downloads\\merged_robot_two_BCs\\payload_design.csv",
    'data_csv': "D:\\Downloads\\merged_robot_two_BCs\\data\\csv\\",
    'data_hist': "D:\\Downloads\\merged_robot_two_BCs\\data\\histograms\\"
}


def get_config():
    return config
