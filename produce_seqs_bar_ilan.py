import collections
import itertools
import math
import random
from typing import Tuple
from Bio.Seq import Seq

from config_seqs import get_configs

def get_amount_of_sequences(number_of_cycles, data_in_kb, bits_per_synthesis_cycle):
    data_in_bit = data_in_kb * 8000
    return math.ceil(data_in_bit / bits_per_synthesis_cycle / number_of_cycles), config['number_of_cycles'] + 2


def dna_sequence_generator(sequence_len=12, symbols=('A', 'C', 'G', 'T')) -> Tuple[str]:
    barcodes = itertools.product(symbols, repeat=sequence_len)
    while True:
        try:
            yield next(barcodes)
        except StopIteration:
            return


def is_forbidden_x_letters_in_a_row(sequence, forbidden_x_letters_in_a_row):
    forbidden_x_letters = [k for k, g in itertools.groupby(sequence) if len(list(g)) >= forbidden_x_letters_in_a_row]
    if len(forbidden_x_letters) > 0:
        return False
    else:
        return True


def hamming_dist(str1, str2):
    i = 0
    count = 0

    while (i < len(str1)):
        if (str1[i] != str2[i]):
            count += 1
        i += 1
    return count


def is_hamming_distance_valid(sequence_list, new_sequence, num_hamming_distance):
    for seq in sequence_list:
        if hamming_dist(new_sequence, seq) < num_hamming_distance:
            return False
    return True


def is_each_letter_has_to_be_x_percent_occurences(new_sequence,
                                                  each_letter_has_to_be_x_percent_occurences_from,
                                                  each_letter_has_to_be_x_percent_occurences_to,
                                                  sequence_len):
    each_letter_has_to_be_x_occurences_from = math.ceil(each_letter_has_to_be_x_percent_occurences_from * sequence_len)
    each_letter_has_to_be_x_occurences_to = math.floor(each_letter_has_to_be_x_percent_occurences_to * sequence_len)
    collection = collections.Counter(new_sequence)
    for key, value in collection.items():
        if value > each_letter_has_to_be_x_occurences_to or value < each_letter_has_to_be_x_occurences_from:
            return False
    return True


def get_dna_experiment_sequence(sequence_len, amount_of_barcodes, amount_of_universal, num_of_building_blocks,
                                num_hamming_distance,
                                forbidden_x_letters_in_a_row, each_letter_has_to_be_x_percent_occurrences_from,
                                each_letter_has_to_be_x_percent_occurrences_to,
                                number_of_cycles):
    f = open(f"sequences/cycles_{number_of_cycles}/num_cycles_{number_of_cycles}_sequences.txt", "w")
    amount_of_sequence = amount_of_barcodes + amount_of_universal + num_of_building_blocks
    sequence_list = []
    # sequence_generator = dna_sequence_generator(sequence_len=sequence_len)

    while amount_of_sequence > 0:
        new_sequence = ''.join(random.choice('ACGT') for _ in range(sequence_len))
        # new_sequence = ''.join(next(sequence_generator))
        if is_forbidden_x_letters_in_a_row(new_sequence, forbidden_x_letters_in_a_row) \
                and is_hamming_distance_valid(sequence_list, new_sequence, num_hamming_distance) \
                and is_each_letter_has_to_be_x_percent_occurences(new_sequence,
                                                                  each_letter_has_to_be_x_percent_occurrences_from,
                                                                  each_letter_has_to_be_x_percent_occurrences_to,
                                                                  sequence_len):
            amount_of_sequence = amount_of_sequence - 1
            sequence_list.append(new_sequence)

            f.write(f'{new_sequence}\n')
            print(f'amount_of_sequence = {amount_of_sequence}')
    sequence_list_barcode = sequence_list[:amount_of_barcodes]
    sequence_list_universal = sequence_list[amount_of_barcodes:amount_of_barcodes + amount_of_universal]
    sequence_list_building_blocks = sequence_list[amount_of_barcodes + amount_of_universal:]
    f.close()
    return sequence_list_barcode, sequence_list_universal, sequence_list_building_blocks


def make_gibson_sequences(sequence_list_barcode,
                          sequence_list_universal,
                          sequence_list_building_blocks,
                          number_of_cycles,
                          ):
    sequence_list_universal1 = [sequence_list_universal[0]] * len(sequence_list_barcode)
    sequence_list_with_universal1_barcode = [','.join(unbc) for unbc in
                                             zip(sequence_list_universal1, sequence_list_barcode)]
    sequence_list_universal2 = [sequence_list_universal[1]] * len(sequence_list_barcode)
    sequence_list_universal1_barcode_universal2 = [','.join(unbc) for unbc in
                                                   zip(sequence_list_with_universal1_barcode, sequence_list_universal2)]

    with open(f'sequences/cycles_{number_of_cycles}/gibson_num_cycles_{number_of_cycles}_barcodes', 'w') as filehandle:
        for bc in sequence_list_universal1_barcode_universal2:
            filehandle.write('%s\n' % bc)

    for universal_index in range(1, len(sequence_list_universal) - 1):
        sequence_list_universal1 = [sequence_list_universal[universal_index]] * len(sequence_list_building_blocks)
        sequence_list_with_universal_payload = [','.join(unbc) for unbc in
                                                zip(sequence_list_universal1, sequence_list_building_blocks)]
        sequence_list_universal2 = [sequence_list_universal[universal_index + 1]] * len(sequence_list_building_blocks)
        sequence_list_with_universal_payload = [','.join(unbc) for unbc in
                                                zip(sequence_list_with_universal_payload, sequence_list_universal2)]
        with open(f'sequences/cycles_{number_of_cycles}/gibson_num_cycles_{number_of_cycles}_payloads{universal_index}',
                  'w') as filehandle:
            for payload in sequence_list_with_universal_payload:
                filehandle.write('%s\n' % payload)

    return


def make_sticky_ends_ligation_sequences(sequence_list_barcode,
                                        sequence_list_universal,
                                        sequence_list_building_blocks,
                                        number_of_cycles,
                                        stick_ends_overhang,
                                        sequence_len):
    sequence_list_universal1 = [sequence_list_universal[0]] * len(sequence_list_barcode)
    sequence_list_with_universal1_barcode = [','.join(unbc) for unbc in
                                             zip(sequence_list_universal1, sequence_list_barcode)]
    sequence_list_universal2 = [sequence_list_universal[1]] * len(sequence_list_barcode)
    sequence_list_universal1_barcode_universal2 = [','.join(unbc) for unbc in
                                                   zip(sequence_list_with_universal1_barcode, sequence_list_universal2)]

    with open(f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_barcodes',
              'w') as filehandle:
        for bc in sequence_list_universal1_barcode_universal2:
            filehandle.write('%s\n' % bc)
    with open(f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_barcodes_rc',
              'w') as filehandle:
        for bc in sequence_list_universal1_barcode_universal2:
            bc_rc = 'P,' + Seq(bc[:(-stick_ends_overhang)]).reverse_complement()
            filehandle.write('%s\n' % bc_rc)

    for universal_index in range(1, len(sequence_list_universal) - 2):
        sequence_list_universal1 = [sequence_list_universal[universal_index]] * len(sequence_list_building_blocks)
        sequence_list_with_universal_payload = [','.join(unbc) for unbc in
                                                zip(sequence_list_universal1, sequence_list_building_blocks)]
        sequence_list_universal2 = [sequence_list_universal[universal_index + 1]] * len(sequence_list_building_blocks)
        sequence_list_with_universal_payload = [','.join(unbc) for unbc in
                                                zip(sequence_list_with_universal_payload, sequence_list_universal2)]
        with open(f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_payloads{universal_index}',
                  'w') as filehandle:
            for payload in sequence_list_with_universal_payload:
                payload = 'P,' + payload[sequence_len:]
                filehandle.write('%s\n' % payload)
        with open(f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_payloads{universal_index}_rc',
                  'w') as filehandle:
            for payloads in sequence_list_with_universal_payload:
                payload_rc = 'P,' + Seq(payloads[(sequence_len-stick_ends_overhang):(-stick_ends_overhang)]).reverse_complement()
                filehandle.write('%s\n' % payload_rc)


    sequence_list_universal1 = [sequence_list_universal[-2]] * len(sequence_list_building_blocks)
    sequence_list_with_universal1_last_payload = [','.join(unbc) for unbc in
                                             zip(sequence_list_universal1, sequence_list_building_blocks)]
    sequence_list_universal2 = [sequence_list_universal[-1]] * len(sequence_list_building_blocks)
    sequence_list_universal1_last_payload_universal2 = [','.join(unbc) for unbc in
                                                   zip(sequence_list_with_universal1_last_payload, sequence_list_universal2)]

    with open(
            f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_payloads{sequence_list_universal[-1]}',
            'w') as filehandle:
        for payload in sequence_list_universal1_last_payload_universal2:
            payload = 'P,' + payload[sequence_len:]
            filehandle.write('%s\n' % payload)
    with open(
            f'sequences/cycles_{number_of_cycles}/sticky_ends_num_cycles_{number_of_cycles}_payloads{sequence_list_universal[-1]}_rc',
            'w') as filehandle:
        for payloads in sequence_list_universal1_last_payload_universal2:
            payload_rc = 'P,' + Seq(
                payloads[(sequence_len - stick_ends_overhang):]).reverse_complement()
            filehandle.write('%s\n' % payload_rc)

    return


def run():
    configs = get_configs()
    for config in configs:
        if(configs.__sizeof__() > 4):
            config = configs
        amount_of_barcodes, amount_of_universal = get_amount_of_sequences(number_of_cycles=config['number_of_cycles'],
                                                                          data_in_kb=config['data_in_kb'],
                                                                          bits_per_synthesis_cycle=config[
                                                                              'bits_per_synthesis_cycle'])
        sequence_list_barcode, sequence_list_universal, sequence_list_building_blocks = get_dna_experiment_sequence(
            sequence_len=config['sequence_len'],
            amount_of_barcodes=amount_of_barcodes,
            amount_of_universal=amount_of_universal,
            num_of_building_blocks=config['num_of_building_blocks'],
            num_hamming_distance=config['num_hamming_distance'],
            forbidden_x_letters_in_a_row=config['forbidden_x_letters_in_a_row'],
            each_letter_has_to_be_x_percent_occurrences_from=config['each_letter_has_to_be_x_percent_occurrences_from'],
            each_letter_has_to_be_x_percent_occurrences_to=config['each_letter_has_to_be_x_percent_occurrences_to'],
            number_of_cycles=config['number_of_cycles'])

        make_gibson_sequences(sequence_list_barcode=sequence_list_barcode,
                              sequence_list_universal=sequence_list_universal,
                              sequence_list_building_blocks=sequence_list_building_blocks,
                              number_of_cycles=config['number_of_cycles'])
        make_sticky_ends_ligation_sequences(sequence_list_barcode=sequence_list_barcode,
                                            sequence_list_universal=sequence_list_universal,
                                            sequence_list_building_blocks=sequence_list_building_blocks,
                                            number_of_cycles=config['number_of_cycles'],
                                            stick_ends_overhang=config['stick_ends_overhang'],
                                            sequence_len=config['sequence_len'])
