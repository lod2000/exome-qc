import os
import sys
import glob

import pandas

import sample_parser
import analysis

# Get path to this file
src_path = sys.path[0]
# Set path to data directory
data_path = os.path.join(src_path, '..', 'data')
gt_file = glob.glob(os.path.join(data_path, 'ground_truth', '*.xlsx'))[0]
bed_file = glob.glob(os.path.join(data_path, 'panel', '*.bed'))[0]
samples_dir = os.path.join(data_path, 'samples')
# Combined data tab file path
df_file = os.path.join(samples_dir, 'combined.tab')

# Generate combined data file
if not os.path.isfile(df_file):
    sample_parser.combine(gt_file, bed_file, samples_dir)

df = pandas.read_csv(df_file, sep='\t', low_memory=False)
# Convert SAMPLE_IDs to Strings
df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
gt_parsed = sample_parser.parse_gt(gt_file)
match_list = sample_parser.find_matches(
        next(os.walk(samples_dir))[1], gt_parsed)
gt = pandas.concat(
        sample_parser.split_gt(gt_parsed, match_list)
).reset_index(drop=True)
panel = sample_parser.parse_bed(bed_file)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 14)
pandas.set_option('display.max_rows', 13)
pandas.set_option('display.width', 300)

# analysis.generate_combined_caller(df)

analysis_df = analysis.analyze_callers(df, panel, gt)
analysis_file = os.path.join(data_path, 'analysis.csv')
try:
    analysis_df.to_csv(
            analysis_file, sep='\t', encoding='utf-8', index=False
    )
except PermissionError:
    print('Another program is using the file analysis.csv.\
           Please close and try again.')
