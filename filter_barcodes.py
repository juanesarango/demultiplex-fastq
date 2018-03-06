# %%
#!/usr/bin/env python

from os.path import join
from Bio import SeqIO

from barcodes import BARCODES_A, BARCODES_B, BARCODES_C, BARCODES_D


INPUT_DIR = '/ifs/work/leukgen/home/arangooj/single_cell_data'
POOL_NUMBER = {
    'A': 1,
    'B': 2,
    'C': 3,
    'D': 4
}

# Utils
def hamming_distance(pattern1, pattern2):
    if len(pattern1) == len(pattern2):
        return sum([pattern1[index] != pattern2[index] for index in range(len(pattern1))])
    raise Exception('Length of both reads do not match')


# Main Function
def filter_files(pool, barcodes, stats={}):

    number = POOL_NUMBER[pool]
    folder = f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}'

    for cell, barcode in barcodes.items():

        stats = {index: 0 for index in range(7)}

        print(f'Analyzing Cell {cell} - Pool {pool}')

        input = join(
            INPUT_DIR, folder, f'Sample_pool{pool}_R1_{cell}.fastq'
        )

        with open(input, 'rt') as fr1:
            records_gen = SeqIO.parse(fr1, 'fastq')

            for record in records_gen:
                score1 = hamming_distance(barcode, record.seq[:6])
                stats[score1] += 1

        # Print Stats
        print(f'Finish Stats of Cell {cell} - Pool {pool}')
        print(stats)
        stats_pool[cell] = stats

    print(stats_pool)
    return stats_pool


if __name__=="__main__":
    stats_pool = {}
    # stats_pool = filter_files('A', BARCODES_A, stats_pool)
    # stats_pool = filter_files('B', BARCODES_B, stats_pool)
    stats_pool = filter_files('C', BARCODES_C, stats_pool)
    stats_pool = filter_files('D', BARCODES_D, stats_pool)
    print(stats_pool)

