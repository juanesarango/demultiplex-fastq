# %%
#!/usr/bin/env python

from os.path import join

from Bio import SeqIO

from barcodes import BARCODES_A, BARCODES_B, BARCODES_C, BARCODES_D


INPUT_DIR = '/home/arangooj/fastq-trim'
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
def filter_files(pool, barcodes):

    number = POOL_NUMBER[pool]
    folder = f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}'
    records_by_cell = {key: 0 for key in barcodes.keys()}

    barcodes = {'A1': 'ACAACC'}
    stats_pool = {}

    for cell, barcode in barcodes.items():

        stats_r1 = {index: 0 for index in range(7)}
        stats_r2 = {index: 0 for index in range(7)}

        print(f'Analyzing Cell {cell} - Pool {pool}')

        input_r1 = join(
            INPUT_DIR, folder, f'Sample_pool{pool}_R1_{cell}.fastq'
        )
        input_r2 = join(
            INPUT_DIR, folder, f'Sample_pool{pool}_R2_{cell}.fastq'
        )

        # print(input_r1)
        # print(input_r1)
        # input_r1 = 'sampleA.fastq'
        # input_r2 = 'sampleA.fastq'

        with open(input_r1, 'rt') as fr1, open(input_r2, 'rt') as fr2:

            records_r1_gen = SeqIO.parse(fr1, 'fastq')
            records_r2_gen = SeqIO.parse(fr2, 'fastq')

            for record_r1 in records_r1_gen:
                record_r2 = next(records_r2_gen)

                score1 = hamming_distance(barcode, record_r1.seq[:6])
                score2 = hamming_distance(barcode, record_r1.seq[:6])

                stats_r1[score1] += 1
                stats_r2[score2] += 1

        # Print Stats
        print(f'Finish Stats of Cell {cell} - Pool {pool}')
        print(stats_r1)
        print(stats_r2)


if __name__=="__main__":
    filter_files('A', BARCODES_A)


