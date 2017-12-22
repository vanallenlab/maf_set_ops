from mafSetOps.MafSetOps import MafSetOps
import argparse

"""
Example usage: 
    python maf_set_intersect.py GCT001-TP-NT-SM-DPBZL-SM-DPBZK.snp.wgs.maf.annotated GCT001-TP-NT-SM-DPC5P-SM-DPBZK.snp.wgs.maf.annotated /output/path
"""


def main():
    parser = argparse.ArgumentParser(description='Find the set intersection between variants in two maf files')
    parser.add_argument('maf_file_1', metavar='maf_file_1', type=str)
    parser.add_argument('maf_file_2', metavar='maf_file_2', type=str)
    parser.add_argument('output_path', metavar='output_path', type=str)

    args = parser.parse_args()

    maf_1 = args.maf_file_1
    maf_2 = args.maf_file_2
    output_path = args.output_path

    MafSetOps.intersection(maf_1, maf_2, output_path=output_path)


if __name__ == '__main__':
    main()
