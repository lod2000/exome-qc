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
# Parsed sample data file path
df_file = os.path.join(output_dir, 'parsed.tab')
# Parsed sample file including combined caller
combined_file = os.path.join(output_dir, 'parsed_combined.tab')

# Get small panel coverage file
panel = parser.parse_bed(bed_file)
    
# Get combined callers
if os.path.isfile(combined_file):
    print('Reading combined callers data...')
    df = pandas.read_csv(combined_file, sep='\t', low_memory=False)
    # Convert SAMPLE_IDs to Strings
    df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
else:
    # Generate combined data file
    if os.path.isfile(df_file):
        print('Reading parsed samples file...')
        # Read parsed CSV
        df = pandas.read_csv(df_file, sep='\t', low_memory=False)
        # Convert SAMPLE_IDs to Strings
        df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
    else:
        df = parser.combine('2017-06-30_NgsReviewer_master', bed_file, samples_dir)

    # Add combined callers to df
    print('Adding combined callers...')
    analysis.add_x_or_more(df)
    analysis.add_prob_caller(df)
    analysis.add_weight_caller(df)
    print('Classifying combined calls...')
    parser.classify(df, parser.get_new_caller_names(df))
    print('Creating new tab file...')
    df.to_csv(
            os.path.join(output_dir, 'parsed_combined.tab'), sep='\t', 
            encoding='utf-8', index=False
    )

analysis_df = analysis.analyze_callers(df, panel)
analysis.plot_callers(df, analysis_df, output_dir)
analysis_file = os.path.join(output_dir, 'analysis.csv')
try:
    analysis_df.to_csv(
            analysis_file, sep='\t', encoding='utf-8', index=False
    )
    print('Output to file ' + analysis_file)
except PermissionError:
    print('Another program is using the file analysis.csv or combined.tab.\
           Please close and try again.')
