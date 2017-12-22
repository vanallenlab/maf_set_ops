import pandas as pd
import sys
import os


class MafSetOps:
    """A class that can be used to compare MAF files using set operations"""
    @staticmethod
    def __get_variant_hash(maf_row):
        return '{}:{}:{}-{}:{}:{}'.format(
            maf_row.Hugo_Symbol,
            maf_row.Chromosome,
            maf_row.Start_position,
            maf_row.End_position,
            maf_row.Reference_Allele,
            maf_row.Tumor_Seq_Allele2
        )

    @staticmethod
    def intersection(maf_file_1, maf_file_2, output_path=None):
        """Determine which variants are found in both MAF files"""
        intersection_variants = []

        maf_variants = {}

        maf_files = [maf_file_1, maf_file_2]
        maf_file_lengths = []

        sys.stdout.write('Reading from MAF {}\n'.format(maf_files[0]))
        df = pd.read_csv(maf_files[0], sep='\t', engine='python')
        maf_file_lengths.append(len(df))
        for (idx, maf_row) in df.iterrows():
            variant_hash = MafSetOps.__get_variant_hash(maf_row)
            maf_variants[variant_hash] = {'count': 1, 'row': maf_row}

        for other_maf in maf_files[1:]:
            sys.stdout.write('Reading from MAF {}\n'.format(other_maf))
            df = pd.read_csv(other_maf, sep='\t', engine='python')
            maf_file_lengths.append(len(df))
            for (idx, maf_row) in df.iterrows():
                variant_hash = MafSetOps.__get_variant_hash(maf_row)
                # If the variant hash is found in
                if variant_hash in maf_variants.keys():
                    maf_variants[variant_hash]['count'] += 1

        for maf_file_name, variant_count in zip(maf_files, maf_file_lengths):
            sys.stdout.write('{} variants in {}\n'.format(variant_count, maf_file_name))

        num_mafs = len(maf_files)
        for variant_hash, variant_info in maf_variants.items():
            count = variant_info['count']
            row = variant_info['row']
            if count == num_mafs:
                intersection_variants.append(row)

        sys.stdout.write('{} variants shared across all MAFs\n'.format(len(intersection_variants)))
        output_filename = '{}.{}'\
            .format(('.'.join([os.path.basename(n).split('.')[0] for n in maf_files])), 'intersection.maf')

        output_location = os.path.join(output_path, output_filename)
        intersection_df = pd.concat(intersection_variants, axis=1)
        intersection_df.transpose().to_csv(output_location, sep='\t', index=False)

