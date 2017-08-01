import os
import sys
import glob

import pandas
import numpy

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
    # Get DataFrame of covered variants
    covered = df.iloc[[
            i for i in range(0, df.shape[0])
            if df['COVERED'][i] or df['REPORTABLE'][i]
    ]].reset_index(drop=True)
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

    # Get DataFrame of covered variants
    covered = df.iloc[[
            i for i in range(0, df.shape[0])
            if df['COVERED'][i] or df['REPORTABLE'][i]
    ]].reset_index(drop=True)
    # Add combined callers to df
    print('Adding combined callers...')
    analysis.add_x_or_more(df)
    training = covered[int(covered.shape[0] / 2):].reset_index(drop=True)
    print(covered.shape[0])
    print(training.shape[0])

    # Generate probability combined caller
    prob_cutoffs = numpy.arange(0, 1, 0.01)
    prob_vals = analysis.prob_fn(prob_cutoffs, training)
    prob_ratios = [
            analysis.get_mcc(
                    prob_vals['TP'][i], prob_vals['TN'][i],
                    prob_vals['FP'][i], prob_vals['FN'][i]
            )
            for i in range(0, len(prob_vals['FP']))
    ]
    prob_cutoff = prob_cutoffs[[
            i for i, ratio in enumerate(prob_ratios) 
            if ratio == max(prob_ratios)
    ][-1]]
    analysis.add_prob_caller(prob_cutoff, df)

    # Generate weights combined caller by calculating optimum cutoff
    weights = analysis.get_caller_weights(df)
    weight_cutoffs = numpy.arange(0, sum(weights), 0.01)
    weight_vals = analysis.weight_fn(weight_cutoffs, training, weights)
    weight_ratios = [
            analysis.get_mcc(
                    weight_vals['TP'][i], weight_vals['TN'][i],
                    weight_vals['FP'][i], weight_vals['FN'][i]
            )
            for i in range(0, len(weight_vals['FP']))
    ]
    weight_cutoff = weight_cutoffs[[
            i for i, ratio in enumerate(weight_ratios)
            if ratio == max(weight_ratios)
    ][-1]]
    analysis.add_weight_caller(weight_cutoff, df, weights)
    print('Classifying combined calls...')
    parser.classify(df, parser.get_new_caller_names(df))
    print('Creating new tab file...')
    df.to_csv(
            os.path.join(output_dir, 'parsed_combined.tab'), sep='\t', 
            encoding='utf-8', index=False
    )

training = covered[:int(covered.shape[0] / 2)].reset_index(drop=True)
analysis_df = analysis.analyze_callers(training, panel)
analysis.plot_callers(training, analysis_df, output_dir)
analysis_file = os.path.join(output_dir, 'analysis.csv')
try:
    analysis_df.to_csv(
            analysis_file, sep='\t', encoding='utf-8', index=False
    )
    print('Output to file ' + analysis_file)
except PermissionError:
    print('Another program is using the file analysis.csv or combined.tab.\
           Please close and try again.')
