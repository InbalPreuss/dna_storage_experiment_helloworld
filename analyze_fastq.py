import csv
import heapq
import math

from Bio import SeqIO, Seq
from matplotlib import pyplot as plt
import numpy as np, pandas as pd
from tqdm import tqdm

from config_analyze import get_config

config = get_config()


def compare_with_errors(read, seq, max_dist=3):
    return sum(r != s for r, s in zip(read, seq)) <= max_dist


def verify_const_universal_and_reverse_complement(read, const_design):
    if verify_const_universal(read, const_design):
        return read
    # elif verify_const_universal(read.reverse_complement(), const_design):
    #     return read.reverse_complement()
    return SeqIO.SeqRecord("None")


# Verify universal
def verify_const_universal(read, const_design):
    for s_idx, s in const_design.iterrows():
        pos = int(s['Pos'])
        if not compare_with_errors(read[pos:pos + 20], s['Seq']):
            return False
    return True

def identify_oligo(read, payload_design, payload_pos, barcodes_design):
    bc = identify_barcode(read, barcodes_design)
    payload = identify_payload(read, payload_design, payload_pos)
    payload.insert(0,bc)
    return payload

# Identify which payload in each position
def identify_payload(read, payload_design, payload_pos):
    res = [0, ] * len(payload_pos)
    for p_idx, pos in enumerate(payload_pos):
        for s_idx, s in payload_design.iterrows():
            if compare_with_errors(read[pos:pos + 20], s['Seq']):
                res[p_idx] = s_idx
                break
    return res


# Identify which barcode in each position
def identify_barcode(read, barcodes_design):
    pos = config['barcode_len']
    for s_idx, s in barcodes_design.iterrows():
        if compare_with_errors(read[pos:pos + 20], s['Seq']):
            return s_idx
    return 0


def open_fastq(input_file):
    with open(input_file, 'r') as inf:
        reads = list(SeqIO.parse(inf, 'fastq'))

    return reads


def reads_len_his(reads):
    len_reads = len(reads)
    print(len_reads)
    lens = [len(r) for r in reads]
    plt.hist(lens, bins=50)
    plt.savefig(config['data_hist'] + '/len_reads_hist.png')
    # plt.show()


def upload_design(const_design_file, payload_design_file, barcodes_design_file):
    const_design_pd = pd.read_csv(const_design_file, index_col=0, header=0)
    payload_design_pd = pd.read_csv(payload_design_file, index_col=0, header=0)
    barcodes_design_pd = pd.read_csv(barcodes_design_file, index_col=0, header=0)

    return const_design_pd, payload_design_pd, barcodes_design_pd


def retrieve_good_reads(reads):
    good_reads = [r for r in reads if len(r) == config['design_len']]
    return good_reads


def good_reads_results_to_csv(reads, const_design, payload_design, barcodes_design):
    res = list()
    failed = 0
    none = SeqIO.SeqRecord("None")
    results_path = config['data_csv'] + '/results.csv'

    for read_idx, read in enumerate(reads):
        if read_idx % 1000 == 999:
            print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
            with open(results_path, "ab") as f:
                f.write(b"\n")
                np.savetxt(f, res, fmt='%i', delimiter=",")
            res = list()
        read = verify_const_universal_and_reverse_complement(read, const_design)
        if read.seq == none.seq:
            failed += 1
            continue

        res.append(identify_oligo(read, payload_design, config['payload_pos'], barcodes_design))
        # res.append(identify_payload(read, payload_design, config['payload_pos']))
        if res[-1].__contains__(0):
            failed += 1
    print(f'processed {read_idx + 1} reads, {failed} ({100 * failed / (read_idx + 1) : .2f}%) failed')
    np.savetxt(config['data_csv'] + '/results.csv', res, fmt='%i', delimiter=",")


def analyze_results_good_reads(results_path, barcodes_design):
    import pandas as pd
    # read your csv to a dataframe
    df = pd.read_csv(results_path)
    dict_bc = {}
    count = 0

    for i in range(168):
        dict_bc_i={}
        for payload in['p1', 'p2', 'p3', 'p4']:
            dict_p = {payload: {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0}}    # for bc, bc_len in barcodes_design.iterrows():
            dict_bc_i.update(dict_p)
        dict_bc[i]=dict_bc_i
    df_bc = df.sort_values(["bc", "p1","p2","p3", "p4"], ascending=True, key=np.sin)
    df_bc.fillna(0)
    # df_bc = pd.to_numeric(df_bc, errors='coerce')
    for row_indx, row in tqdm(df_bc.iterrows()):
        for location_payload, payload in row[1:].items():
            # print(f'payload={payload}, location_payload={location_payload}, int(row[bc])={r}, dict_bc[int(row[bc])][location_payload][int(payload)]={dbc}' )

            try:
                # if math.isnan(int(payload)) or math.isnan(int(row['bc'])):
                #     continue
                # elif math.isnan(dict_bc[int(row['bc'])][location_payload][int(payload)]):
                #     continue
                dict_bc[int(row['bc'])][location_payload][int(payload)] += 1

            except:
                print(f"An exception occurred {row_indx}, row={row}")
                continue

            dict_bc[int(row['bc'])][location_payload][int(payload)] += 1
        # df_bc['total_'+p] = df_bc.p.apply(lambda x: df_bc.p.value_counts()[x])
        # print('start')
        # df[['bc', 'total_'+p]].to_csv(config['data_csv']+"total_"+p+".csv")
        # print(df_bc)

    with open(config['data_csv'] + 'results_count.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
        w = csv.DictWriter(f, dict_bc.keys())
        w.writeheader()
        w.writerow(dict_bc)

    return dict_bc


def most_common_for_each_bc(dict_bc):
    for bc_i in range(168):
        for payload in ['p1', 'p2', 'p3', 'p4']:
            most_common = heapq.nlargest(5, dict_bc[bc_i][payload], key=dict_bc[bc_i][payload].get)
            print(f'bc = {bc_i}, payload = {payload}, 5 most common = {most_common}')


def run():
    # upload design
    const_design_pd, payload_design_pd, barcodes_design_pd = upload_design(
        const_design_file=config['const_design_file'],
        payload_design_file=config['payload_design_file'],
        barcodes_design_file=config['barcodes_design_file'])
    #
    # # reads
    # reads = open_fastq(input_file=config['input_file'])
    #
    # # reads len showed in histogram
    # reads_len_his(reads=reads)
    #
    # # good reads with len  220
    # good_reads = retrieve_good_reads(reads=reads)
    #
    # # Write the good reads with len 220 to results.csv
    # good_reads_results_to_csv(reads=good_reads,
    #                           const_design=const_design_pd,
    #                           payload_design=payload_design_pd,
    #                           barcodes_design=barcodes_design_pd)

    results_path = config['data_csv'] + 'results_temp.csv'
    # Analyze results of good reads
    dict_bc=analyze_results_good_reads(results_path=results_path,barcodes_design=barcodes_design_pd)
    most_common_for_each_bc(dict_bc)
