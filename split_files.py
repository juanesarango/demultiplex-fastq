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
        return sum([
            pattern1[index] != pattern2[index]
            for index in range(len(pattern1))
        ])
    raise Exception('Length of both reads do not match')


# Main Function
def split_files(pool, barcodes, stats_pool):

    number = POOL_NUMBER[pool]
    folder_in = (
        f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}'
    )
    folder_out = (
        f'Sample_171205-HBM-0006-pool{pool}_IGO_08099_C_{number}_trimmed'
    )

    for cell, barcode in barcodes.items():

        print(f'Spliting Cell {cell} - Pool {pool}')

        # Define input/output files
        input_file_r1 = join(
            INPUT_DIR, folder_in, f'Sample_pool{pool}_R1_{cell}.fastq'
        )
        input_file_r2 = join(
            INPUT_DIR, folder_in, f'Sample_pool{pool}_R2_{cell}.fastq'
        )
        output_file_r1 = join(
            INPUT_DIR, folder_out, f'Sample_pool{pool}_R1_{cell}_trimmed.fastq'
        )
        output_file_r2 = join(
            INPUT_DIR, folder_out, f'Sample_pool{pool}_R2_{cell}_trimmed.fastq'
        )

        with open(input_file_r1, 'rt') as fr1, open(input_file_r2, 'rt') as fr2:

            filtered_seq_r1 = []
            filtered_seq_r2 = []

            records_r1_gen = SeqIO.parse(fr1, 'fastq')
            records_r2_gen = SeqIO.parse(fr2, 'fastq')

            number_of_sequences = 0

            # Evaluate R1 and R2 sequences one by one at the same time
            for record_r1 in records_r1_gen:
                record_r2 = next(records_r2_gen)

                number_of_sequences += 1

                mismatches_r1 = hamming_distance(record_r1[:6], barcode)
                mismatches_r2 = hamming_distance(record_r2[:6], barcode)

                # Include only the trimmed sequences with 0 or 1 mismatch
                if mismatches_r1 <= 1 or mismatches_r2 <= 1:
                    filtered_seq_r1.append(record_r1[6:])
                    filtered_seq_r2.append(record_r2[6:])

            # Show some stats
            stats_pool[cell] = {
                'total': number_of_sequences,
                'filtered': len(filtered_seq_r1),
                'percentage':round(
                    len(filtered_seq_r1)/number_of_sequences * 100, 2
                )
            }
            print(stats_pool[cell])

            # Write to files
            SeqIO.write(filtered_seq_r1, output_file_r1, 'fastq')
            SeqIO.write(filtered_seq_r2, output_file_r2, 'fastq')

        print(f'Finished Pool {pool}!')
        print(stats_pool)


if __name__ == "__main__":
    stats_pool = {}
    split_files('A', BARCODES_A, stats_pool)
    split_files('B', BARCODES_B, stats_pool)
    split_files('C', BARCODES_C, stats_pool)
    split_files('D', BARCODES_D, stats_pool)
    print(stats_pool)


"""
Filtering Stats:

stats = {
    'A1': {'total': 979515, 'filtered': 921027, 'percentage': 94.03},
    'A2': {'total': 1458954, 'filtered': 1410721, 'percentage': 96.69},
    'A3': {'total': 2356877, 'filtered': 2292536, 'percentage': 97.27},
    'A4': {'total': 1192324, 'filtered': 1128004, 'percentage': 94.61},
    'A5': {'total': 864906, 'filtered': 806501, 'percentage': 93.25},
    'A6': {'total': 1862990, 'filtered': 1821791, 'percentage': 97.79},
    'A7': {'total': 833566, 'filtered': 784912, 'percentage': 94.16},
    'A8': {'total': 1336569, 'filtered': 1296804, 'percentage': 97.02},
    'A9': {'total': 5188707, 'filtered': 5124591, 'percentage': 98.76},
    'A10': {'total': 869596, 'filtered': 827066, 'percentage': 95.11},
    'A11': {'total': 1885650, 'filtered': 1836159, 'percentage': 97.38},
    'A12': {'total': 2640814, 'filtered': 2581124, 'percentage': 97.74},
    'B1': {'total': 1461334, 'filtered': 1415157, 'percentage': 96.84},
    'B2': {'total': 1430574, 'filtered': 1379040, 'percentage': 96.4},
    'B3': {'total': 1461151, 'filtered': 1425033, 'percentage': 97.53},
    'B4': {'total': 1387233, 'filtered': 1341448, 'percentage': 96.7},
    'B5': {'total': 787276, 'filtered': 745782, 'percentage': 94.73},
    'B6': {'total': 2516497, 'filtered': 2489294, 'percentage': 98.92},
    'B7': {'total': 3051455, 'filtered': 2973889, 'percentage': 97.46},
    'B8': {'total': 2534009, 'filtered': 2489580, 'percentage': 98.25},
    'B9': {'total': 2196957, 'filtered': 2169815, 'percentage': 98.76},
    'B10': {'total': 1672517, 'filtered': 1628737, 'percentage': 97.38},
    'B11': {'total': 2238618, 'filtered': 2205269, 'percentage': 98.51},
    'B12': {'total': 2239936, 'filtered': 2204463, 'percentage': 98.42},
    'C1': {'total': 1512722, 'filtered': 1469875, 'percentage': 97.17},
    'C2': {'total': 1097656, 'filtered': 1062505, 'percentage': 96.8},
    'C3': {'total': 1199330, 'filtered': 1154391, 'percentage': 96.25},
    'C4': {'total': 1345924, 'filtered': 1312012, 'percentage': 97.48},
    'C5': {'total': 721761, 'filtered': 676816, 'percentage': 93.77},
    'C6': {'total': 1427646, 'filtered': 1396440, 'percentage': 97.81},
    'C7': {'total': 855286, 'filtered': 820043, 'percentage': 95.88},
    'C8': {'total': 1639615, 'filtered': 1609686, 'percentage': 98.17},
    'C9': {'total': 2516984, 'filtered': 2471237, 'percentage': 98.18},
    'C10': {'total': 803838, 'filtered': 770980, 'percentage': 95.91},
    'C11': {'total': 1201528, 'filtered': 1167304, 'percentage': 97.15},
    'C12': {'total': 1826806, 'filtered': 1786706, 'percentage': 97.8},
    'D1': {'total': 827358, 'filtered': 794256, 'percentage': 96.0},
    'D2': {'total': 862868, 'filtered': 826832, 'percentage': 95.82},
    'D3': {'total': 759212, 'filtered': 732260, 'percentage': 96.45},
    'D4': {'total': 1961655, 'filtered': 1926536, 'percentage': 98.21},
    'D5': {'total': 453363, 'filtered': 422042, 'percentage': 93.09},
    'D6': {'total': 1516236, 'filtered': 1496135, 'percentage': 98.67},
    'D7': {'total': 2341031, 'filtered': 2287950, 'percentage': 97.73},
    'D8': {'total': 1589464, 'filtered': 1557651, 'percentage': 98.0},
    'D9': {'total': 1831497, 'filtered': 1811784, 'percentage': 98.92},
    'D10': {'total': 1376601, 'filtered': 1345684, 'percentage': 97.75},
    'D11': {'total': 2287288, 'filtered': 2260895, 'percentage': 98.85},
    'D12': {'total': 4153978, 'filtered': 4126127, 'percentage': 99.33},
    'E1': {'total': 4087374, 'filtered': 3931511, 'percentage': 96.19},
    'E2': {'total': 2535006, 'filtered': 2416961, 'percentage': 95.34},
    'E3': {'total': 5127710, 'filtered': 4954498, 'percentage': 96.62},
    'E4': {'total': 3666963, 'filtered': 3552579, 'percentage': 96.88},
    'E5': {'total': 2171522, 'filtered': 2043669, 'percentage': 94.11},
    'E6': {'total': 3260881, 'filtered': 3155742, 'percentage': 96.78},
    'E7': {'total': 1008226, 'filtered': 889862, 'percentage': 88.26},
    'E8': {'total': 2994327, 'filtered': 2891047, 'percentage': 96.55},
    'E9': {'total': 6435796, 'filtered': 6281949, 'percentage': 97.61},
    'E10': {'total': 1843911, 'filtered': 1716700, 'percentage': 93.1},
    'E11': {'total': 4000474, 'filtered': 3871461, 'percentage': 96.78},
    'E12': {'total': 20144524, 'filtered': 19972071, 'percentage': 99.14},
    'F1': {'total': 1348877, 'filtered': 1220088, 'percentage': 90.45},
    'F2': {'total': 1635533, 'filtered': 1499894, 'percentage': 91.71},
    'F3': {'total': 3590634, 'filtered': 3490541, 'percentage': 97.21},
    'F4': {'total': 4438806, 'filtered': 4323177, 'percentage': 97.4},
    'F5': {'total': 2223385, 'filtered': 2105221, 'percentage': 94.69},
    'F6': {'total': 4319206, 'filtered': 4254285, 'percentage': 98.5},
    'F7': {'total': 6315481, 'filtered': 6149325, 'percentage': 97.37},
    'F8': {'total': 4785108, 'filtered': 4667870, 'percentage': 97.55},
    'F9': {'total': 1973297, 'filtered': 1911229, 'percentage': 96.85},
    'F10': {'total': 3783472, 'filtered': 3678118, 'percentage': 97.22},
    'F11': {'total': 6770083, 'filtered': 6660286, 'percentage': 98.38},
    'F12': {'total': 18269344, 'filtered': 18150198, 'percentage': 99.35},
    'G1': {'total': 973631, 'filtered': 931625, 'percentage': 95.69},
    'G2': {'total': 389498, 'filtered': 354837, 'percentage': 91.1},
    'G3': {'total': 629533, 'filtered': 585878, 'percentage': 93.07},
    'G4': {'total': 923959, 'filtered': 888843, 'percentage': 96.2},
    'G5': {'total': 715666, 'filtered': 675968, 'percentage': 94.45},
    'G6': {'total': 892856, 'filtered': 860633, 'percentage': 96.39},
    'G7': {'total': 399522, 'filtered': 363651, 'percentage': 91.02},
    'G8': {'total': 693722, 'filtered': 664478, 'percentage': 95.78},
    'G9': {'total': 1533084, 'filtered': 1487683, 'percentage': 97.04},
    'G10': {'total': 701222, 'filtered': 668796, 'percentage': 95.38},
    'G11': {'total': 1147998, 'filtered': 1105771, 'percentage': 96.32},
    'G12': {'total': 3510720, 'filtered': 3469098, 'percentage': 98.81},
    'H1': {'total': 766380, 'filtered': 729979, 'percentage': 95.25},
    'H2': {'total': 701858, 'filtered': 662047, 'percentage': 94.33},
    'H3': {'total': 1259075, 'filtered': 1232856, 'percentage': 97.92},
    'H4': {'total': 479193, 'filtered': 447566, 'percentage': 93.4},
    'H5': {'total': 747185, 'filtered': 713966, 'percentage': 95.55},
    'H6': {'total': 652355, 'filtered': 633651, 'percentage': 97.13},
    'H7': {'total': 1706330, 'filtered': 1662641, 'percentage': 97.44},
    'H8': {'total': 1361278, 'filtered': 1331245, 'percentage': 97.79},
    'H9': {'total': 2340887, 'filtered': 2324405, 'percentage': 99.3},
    'H10': {'total': 3405120, 'filtered': 3370477, 'percentage': 98.98},
    'H11': {'total': 5032458, 'filtered': 5000336, 'percentage': 99.36},
    'H12': {'total': 1881020, 'filtered': 1855736, 'percentage': 98.66}
}
"""
