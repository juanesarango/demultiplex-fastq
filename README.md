# Demultiplex-fastq

These are the scripts used for demultiplexing the pool of samples from fastq files.

## How they work

### demultiplex_files.py

Each sample has R1 and R2 files as sequences are pair-ended.
Each read has 3 lines for:

- **header**: @SN835:827:H5MH5BCX2:1:1107:3300:2125 1:N:0:TAAGGCGA
- **sequence**: NTAAGGCCTCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGG
- **quality score**: #<<DBFGHIIIIHHIIHIIEHIIHIIHIIIHHIHIHHIHIHHHHHIIIICE

So a R1 sequence example is:

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 1:N:0:TAAGGCGA
    NTAAGGCCTCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGG
    +
    #<<DBFGHIIIIHHIIHIIEHIIHIIHIIIHHIHIHHIHIHHHHHIIIICE

And its corresponding R2 sequence will have the same header:

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 2:N:0:TAAGGCGA
    GTGAGGCCTTACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT
    +
    DDDDDHHIHHIIIIHIHIIIIEHEHDHHHIEEHHEHHIIHEHHIGHHIIHI

Now if those files are multiplexed by the following barcodes:

| Cell | Barcode |
| --- | ------ |
| A01 | ACAACC |
| A02 | ACTCAC |
| A03 | AGGATG |
| A04 | ATCGAC |
| A05 | CAAGAG |
| A06 | CATGAC |
| A07 | CCTTCG |
| A08 | CGGTAG |
| A09 | CTATTG |
| A10 | CTCAGC |
| A11 | GCATTC |
| A12 | GTGAGG |

Then the script takes the first 6 base pairs of each reads, **score** them against each of the barcodes taking the quality of each read into account, and classify both reads to the best match.

    barcode_r1 = NTAAGG
    barcode_r2 = GTGAGG

Scoring both, `r1` will have a bad score for A12, but `r2` will have a perfect match with a good score, so both reads will be demultiplexed into the output files:

`R1_A12.fastq`

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 1:N:0:TAAGGCGA
    NTAAGGCCTCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGG
    +
    #<<DBFGHIIIIHHIIHIIEHIIHIIHIIIHHIHIHHIHIHHHHHIIIICE

`R2_A12.fastq`

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 2:N:0:TAAGGCGA
    GTGAGGCCTTACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT
    +
    DDDDDHHIHHIIIIHIHIIIIEHEHDHHHIEEHHEHHIIHEHHIGHHIIHI


### split_files.py

This script does 2 main things:
- It **filters out** all sequences whose barcode differs by 2 or more base pairs with the barcode.
- It **trimms** the barcode of each sequence.

In the previous example as both barcodes differ from the pattern by <2 base pairs. Then the sequences are kept and stored without the barcode:

`R1_A12.fastq`

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 1:N:0:TAAGGCGA
    CCTCACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAAGG
    +
    GHIIIIHHIIHIIEHIIHIIHIIIHHIHIHHIHIHHHHHIIIICE

`R2_A12.fastq`

    @SN835:827:H5MH5BCX2:1:1107:3300:2125 2:N:0:TAAGGCGA
    CCTTACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT
    +
    HIHHIIIIHIHIIIIEHEHDHHHIEEHHEHHIIHEHHIGHHIIHI

