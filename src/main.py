import os
import sys
import glob

import pandas

import parser
import analysis

# Get path to this file
src_path = sys.path[0]
# Set path to data directory
data_path = os.path.join(src_path, '..', 'data')
gt_dir = os.path.join(data_path, 'ground_truth')
gt_file = os.path.join(gt_dir, 'simple_ground_truth.csv')
bed_file = glob.glob(os.path.join(data_path, 'panel', '*.bed'))[0]
samples_dir = os.path.join(data_path, 'samples')
# Combined data tab file path
df_file = os.path.join(samples_dir, 'parsed.tab')
# Combined caller weights file path
weights_file = os.path.join(data_path, 'weights.tab')

# Generate combined data file
if not os.path.isfile(df_file):
    parser.combine('2017-06-30_NgsReviewer_master', bed_file, samples_dir)

# Read parsed CSV
df = pandas.read_csv(df_file, sep='\t', low_memory=False)
# Read ground truth CSV
gt = pandas.read_csv(gt_file, sep='\t', low_memory=False)
# Convert SAMPLE_IDs to Strings
df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
gt['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
panel = parser.parse_bed(bed_file)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 7)
pandas.set_option('display.max_rows', 13)
pandas.set_option('display.width', 300)

analysis_df = analysis.analyze_callers(df, panel, gt)
# analysis.plot_callers(analysis_df)
analysis_file = os.path.join(data_path, 'analysis.csv')
try:
    df.to_csv(
            os.path.join(samples_dir, 'combined.tab'), sep='\t',
            encoding='utf-8', index=False
    )
    print('\nOutput to data/samples/combined.tab')
    analysis_df.to_csv(
            analysis_file, sep='\t', encoding='utf-8', index=False
    )
    print('\nOutput to file ' + analysis_file)
except PermissionError:
    print('Another program is using the file analysis.csv or combined.tab.\
           Please close and try again.')
