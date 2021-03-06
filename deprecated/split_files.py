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
def split_files(pool, barcodes):

    number = POOL_NUMBER[pool]
    folder_in = f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}'
    folder_out = f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}_trimmed'

    for cell, barcode in barcodes.items():

        print(f'Spliting Cell {cell} - Pool {pool}')

        # Define input/output files
        input_file_r1 = join(INPUT_DIR, folder_in, f'Sample_pool{pool}_R1_{cell}.fastq')
        input_file_r2 = join(INPUT_DIR, folder_in, f'Sample_pool{pool}_R2_{cell}.fastq')

        output_file_r1 = join(INPUT_DIR, folder_out, f'Sample_pool{pool}_R1_{cell}_trimmed.fastq')
        output_file_r2 = join(INPUT_DIR, folder_out, f'Sample_pool{pool}_R2_{cell}_trimmed.fastq')

        with open(input_file_r1, 'rt') as fr1, open(input_file_r2, 'rt') as fr2:

            # Filter and trim sequences
            filtered_seq_r1 = [
                record[6:] for record in SeqIO.parse(fr1, 'fastq')
                if hamming_distance(record[:6], barcode) <= 1
            ]
            filtered_seq_r2 = [
                record[6:] for record in SeqIO.parse(fr2, 'fastq')
                if hamming_distance(record[:6], barcode) <= 1
            ]

            # Write to files
            SeqIO.write(filtered_seq_r1, output_file_r1, 'fastq')
            SeqIO.write(filtered_seq_r2, output_file_r2, 'fastq')


if __name__=="__main__":
    split_files('A', BARCODES_A)
    split_files('B', BARCODES_B)
    split_files('C', BARCODES_C)
    split_files('D', BARCODES_D)

