import os
import sys
import glob
import importlib
import re

import pandas
import numpy

import utils
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

# Import joint caller modules
callers_dir = os.path.join(main_dir, 'callers')
sys.path.append(callers_dir)
caller_names = [
        caller.split('.py')[0] for caller in os.listdir(callers_dir)
        if re.search('^[A-Za-z]', caller)
]
if len(caller_names) > 0:
    print('Found mutation caller modules:')
    print(caller_names)
    callers = []
    for i, caller in enumerate(caller_names):
        callers.append( importlib.import_module(caller))

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
    training = covered[int(covered.shape[0] / 2):].reset_index(drop=True)
    # Add joint callers to df
    print('Adding joint callers...')
    for caller in callers:
        caller.add_caller(df, training)

    # Generate probability combined caller
    prob_cutoffs = numpy.arange(0, 1, 0.01)
    prob_vals = analysis.prob_fn(prob_cutoffs, training)
    prob_ratios = [
            utils.get_mcc(
                    prob_vals['TP'][i], prob_vals['TN'][i],
                    prob_vals['FP'][i], prob_vals['FN'][i]
            )
            for i in range(0, len(prob_vals['FP']))
    ]
    prob_cutoff = prob_cutoffs[[
            i for i, ratio in enumerate(prob_ratios) 
            if ratio == max(prob_ratios)
    ][-1]]
    print('Combined cutoff: ' + str(prob_cutoff))
    analysis.add_prob_caller(prob_cutoff, df)
    print('Classifying combined calls...')
    parser.classify(df, utils.get_joint_callers(df))
    print('Creating new tab file...')
    df.to_csv(
            os.path.join(output_dir, 'parsed_combined.tab'), sep='\t', 
            encoding='utf-8', index=False
    )

# Get DataFrame of covered variants
covered = df.iloc[[
        i for i in range(0, df.shape[0])
        if df['COVERED'][i] or df['REPORTABLE'][i]
]].reset_index(drop=True)
training = covered[:int(covered.shape[0] / 2)].reset_index(drop=True)
analysis_df = analysis.analyze_callers(training, panel)
analysis.plot_callers(training, analysis_df, output_dir)
analysis_file = os.path.join(output_dir, 'analysis.csv')
analysis_df.to_csv(
        analysis_file, sep='\t', encoding='utf-8', index=False
)
print('Output to file ' + analysis_file)
