import os
import sys
import glob

import pandas

import parser
import analysis

# Directory structure
main_dir = os.path.abspath(os.path.join(sys.path[0], '..'))
output_dir = os.path.join(main_dir, 'output')
data_path = os.path.join(main_dir, 'data')
# Panel bed file
bed_file = glob.glob(os.path.join(data_path, 'panel', '*.bed'))[0]
# Samples directory
samples_dir = os.path.join(data_path, 'samples')
# Combined data tab file path
df_file = os.path.join(output_dir, 'parsed.tab')
# Combined caller weights file path
weights_file = os.path.join(output_dir, 'weights.tab')

# Generate combined data file
if os.path.isfile(df_file):
    # Read parsed CSV
    df = pandas.read_csv(df_file, sep='\t', low_memory=False)
    # Convert SAMPLE_IDs to Strings
    df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
else:
    df = parser.combine('2017-06-30_NgsReviewer_master', bed_file, samples_dir)

# Get small panel coverage file
panel = parser.parse_bed(bed_file)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 7)
pandas.set_option('display.max_rows', 13)
pandas.set_option('display.width', 300)

analysis_df = analysis.analyze_callers(df, panel)
# analysis.plot_callers(analysis_df)
analysis_file = os.path.join(output_dir, 'analysis.csv')
try:
#    df.to_csv(
#            os.path.join(samples_dir, 'combined.tab'), sep='\t',
#            encoding='utf-8', index=False
#    )
#    print('\nOutput to data/samples/combined.tab')
    analysis_df.to_csv(
            analysis_file, sep='\t', encoding='utf-8', index=False
    )
    print('Output to file ' + analysis_file)
except PermissionError:
    print('Another program is using the file analysis.csv or combined.tab.\
           Please close and try again.')
