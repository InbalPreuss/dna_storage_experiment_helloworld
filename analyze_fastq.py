from Bio import SeqIO, Seq
from matplotlib import pyplot as plt
import numpy as np, pandas as pd
import csv

from config_analyze import get_config
config = get_config()

# TODO: if verify_const is false, then try verify_const with the reverse complement of the read and add it to the list of reads.
def compare_with_errors(read, seq, max_dist=3):
    return sum(r != s for r, s in zip(read, seq)) <= max_dist


def verify_const_universal_and_reverse_complement(read, const_design):
    if verify_const_universal(read, const_design):
        return read
    elif verify_const_universal(read.reverse_complement(), const_design):
        return read.reverse_complement()
    return SeqIO.SeqRecord("None")


# Verify universal
def verify_const_universal(read, const_design):
    for s_idx, s in const_design.iterrows():
        pos = int(s['Pos'])
        if not compare_with_errors(read[pos:pos + 20], s['Seq']):
            return False
    return True

# Identify which payload in each position
def identify_payload(read, payload_design, payload_pos):
    res = [0, ] * len(payload_pos)
    for p_idx, pos in enumerate(payload_pos):
        for s_idx, s in payload_design.iterrows():
            if compare_with_errors(read[pos:pos + 20], s['Seq']):
                res[p_idx] = s_idx
                break
    return res


def open_fastq(input_file):
    with open(input_file, 'r') as inf:
        reads = list(SeqIO.parse(inf, 'fastq'))

    return reads


def reads_len_his(reads):
    len_reads = len(reads)
    print(len_reads)
    lens = [len(r) for r in reads]
    plt.hist(lens, bins=50)
    plt.savefig(config['data_hist']+'/len_reads_hist.png')
    # plt.show()


def upload_design(const_design_file, payload_design_file):
    const_design_pd = pd.read_csv(const_design_file, index_col=0, header=0)
    payload_design_pd = pd.read_csv(payload_design_file, index_col=0, header=0)

    return const_design_pd, payload_design_pd


def retrieve_good_reads(reads):
    good_reads = [r for r in reads if len(r) == config['design_len']]
    return good_reads


def good_reads_results_to_csv(reads, const_design, payload_design):
    res=list()
    failed = 0
    none = SeqIO.SeqRecord("None")

    for r_idx, r in enumerate(reads):
        if r_idx % 1000 == 999:
            print(f'processed {r_idx + 1} reads, {failed} ({100 * failed / (r_idx + 1) : .2f}%) failed')
            np.savetxt(config['data_csv']+'/results.csv', res, fmt='%i', delimiter=",")
            res = list()
        r = verify_const_universal_and_reverse_complement(r, const_design)
        if r.seq == none.seq:
            failed += 1
            continue
        res.append(identify_payload(r, payload_design, config['payload_pos']))
        if res[-1].__contains__(0):
            failed += 1
    print(f'processed {r_idx + 1} reads, {failed} ({100 * failed / (r_idx + 1) : .2f}%) failed')
    np.savetxt(config['data_csv']+'/results.csv', res, fmt='%i', delimiter=",")


def run():

    # upload design
    const_design_pd, payload_design_pd = upload_design(const_design_file=config['const_design_file'], payload_design_file=config['payload_design_file'])

    # reads
    reads = open_fastq(input_file=config['input_file'])

    # reads len showed in histogram
    reads_len_his(reads=reads)

    # good reads with len  220
    good_reads = retrieve_good_reads(reads=reads)

    # Write the good reads with len 220 to results.csv
    good_reads_results_to_csv(reads=good_reads, const_design=const_design_pd, payload_design=payload_design_pd)
