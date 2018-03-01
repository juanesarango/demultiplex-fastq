#!/usr/bin/env python

from os.path import join
from os import listdir

from Bio import SeqIO
import gzip

from barcodes import BARCODES_A, BARCODES_B, BARCODES_C, BARCODES_D


ROOT_DIR_DATA = '/ifs/archive/BIC/share/leukgen/papaemme/LIZ_0827_AH5MH5BCX2/Project_08099_C/'
OUTPUT_DIR = '/home/arangooj/fastq-trim'
POOL_NUMBER = {
    'A': 1,
    'B': 2,
    'C': 3,
    'D': 4
}

# Utils
def similarity_score(pattern, read, quality):
    if len(pattern) == len(read) == len(quality):
        return sum([
            abs((1 if read[index] == pattern[index] else 0) - quality[index]/40) * 100
            for index in range(len(pattern))
        ])
    raise Exception('Length of both reads do not match')


def find_array_positions(num, array):
    return [index for index, item in enumerate(array) if item == num][0]


def find_most_similar_paired(record1, record2, barcodes):
    if record1.id == record2.id:
        sequence1 = record1.seq[:6]
        sequence2 = record2.seq[:6]
        quality1 = record1.letter_annotations['phred_quality'][:6]
        quality2 = record2.letter_annotations['phred_quality'][:6]
        distance = [
            similarity_score(barcode, sequence1, quality1) +
            similarity_score(barcode, sequence2, quality2)
            for barcode in barcodes.values()
        ]
        matched_barcode = [
            barcode
            for barcode in barcodes.items()
        ][find_array_positions(min(distance), distance)]
        return matched_barcode
    raise Exception("R1 read sequence daoesn't match the R2 read sequence.")


# Main Function
def demultiplex(pool, barcodes):

    number = POOL_NUMBER[pool]
    folder = f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}'
    records_by_cell = {key: 0 for key in barcodes.keys()}

    # Get R1 and R2 file paths
    file_dir = join(ROOT_DIR_DATA, folder)

    fastq_r1_path = [
        join(file_dir, file_name)
        for file_name in listdir(file_dir)
        if file_name.endswith('fastq.gz') and 'R1' in file_name
    ][0]

    fastq_r2_path = [
        join(file_dir, file_name)
        for file_name in listdir(file_dir)
        if file_name.endswith('fastq.gz') and 'R2' in file_name
    ][0]

    # Open record one by one of each file, classify it and output to file.
    with gzip.open(fastq_r1_path, 'rt') as fr1, gzip.open(fastq_r2_path, 'rt') as fr2:

        records_r1_gen = SeqIO.parse(fr1, 'fastq')
        records_r2_gen = SeqIO.parse(fr2, 'fastq')

        n = 0

        for record_r1 in records_r1_gen:
            record_r2 = next(records_r2_gen)

            n += 1
            if n> 50: break

            # Classify records according to barcode
            cell, _ = find_most_similar_paired(record_r1, record_r2, barcodes)

            # Write in respective files
            output_r1 = join(
                OUTPUT_DIR, folder, f'Sample_pool{pool}_R1_{cell}.fastq'
                )
            output_r2 = join(
                OUTPUT_DIR, folder, f'Sample_pool{pool}_R1_{cell}.fastq'
                )

            SeqIO.write(record_r1, output_r1, 'fastq')
            SeqIO.write(record_r2, output_r2, 'fastq')

            # Store some stats
            records_by_cell[cell] += 1

    print(f'Finish Demultiplexing of Pool {pool} Sample {number}')
    print(records_by_cell)


if __name__=="__main__":
    demultiplex('A', BARCODES_A)
    demultiplex('B', BARCODES_B)
    demultiplex('C', BARCODES_C)
    demultiplex('D', BARCODES_D)
