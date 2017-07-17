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
df_file = os.path.join(samples_dir, 'parsed.tab')
# Combined caller weights file path
weights_file = os.path.join(data_path, 'weights.tab')

# Generate combined data file
if not os.path.isfile(df_file):
    sample_parser.combine(gt_file, bed_file, samples_dir)

# Read parsed CSV
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

# Generate combined callers
callers = sample_parser.get_og_caller_names(df)
if os.path.isfile(weights_file):
    weights_df = pandas.read_csv(weights_file, sep='\t')
    weights = weights_df['WEIGHT']
else:
    weights = analysis.generate_combined_caller_weights(df)
    weights_df = pandas.DataFrame({'CALLER': callers, 'WEIGHT': weights})
    weights_df.to_csv(weights_file, sep='\t', encoding='utf-8', index=False)

print(weights)
analysis.add_combined_caller(df, weights)
analysis.add_2_or_more(df)

# analysis.add_differences(df)

analysis_df = analysis.analyze_callers(df, panel, gt)
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
