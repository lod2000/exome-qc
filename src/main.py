import os
import sys
import glob
import importlib
import re

import pandas
import numpy
import matplotlib.pyplot as pyplot
import matplotlib
from adjustText import adjust_text

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
    caller_modules = []
    for i, caller in enumerate(caller_names):
        caller_modules.append(importlib.import_module(caller))

# Get small panel coverage file
panel = parser.parse_bed(bed_file)
    
# Get combined callers
#if os.path.isfile(combined_file):
#    print('Reading combined callers data...')
#    df = pandas.read_csv(combined_file, sep='\t', low_memory=False)
#    # Convert SAMPLE_IDs to Strings
#    df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
#    # Get DataFrame of covered variants
#    covered = df.iloc[[
#            i for i in range(0, df.shape[0])
#            if df['COVERED'][i] or df['REPORTABLE'][i]
#    ]].reset_index(drop=True)
#else:

# Read / generate combined data file
if os.path.isfile(df_file):
    print('Reading parsed samples file...')
    # Read parsed CSV
    df = pandas.read_csv(df_file, sep='\t', low_memory=False)
    # Convert SAMPLE_IDs to Strings
    df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
else:
    df = parser.combine('2017-06-30_NgsReviewer_master', bed_file, samples_dir)

# Get indices of covered variants
covered_indices = [
        i for i in range(0, df.shape[0])
        if df['COVERED'][i] or df['REPORTABLE'][i]
]
covered_split = int(len(covered_indices) / 2)
# Get training set
training_indices = covered_indices[covered_split:]
training = df.iloc[training_indices].reset_index(drop=True)
# Add joint callers to df
print('Adding joint callers...')
for caller in caller_modules:
    caller.add_caller(df, training)

# Analyze combined callers
print('Classifying combined calls...')
parser.classify(df, utils.get_joint_callers(df))
print('Creating new tab file...')
df.to_csv(
        os.path.join(output_dir, 'parsed_combined.tab'), sep='\t', 
        encoding='utf-8', index=False
)

# Get testing set
testing_indices = covered_indices[:covered_split]
testing = df.iloc[testing_indices].reset_index(drop=True)

analysis_df = analysis.analyze_callers(testing, panel)
#analysis.plot_callers(testing, analysis_df, output_dir)
analysis_file = os.path.join(output_dir, 'analysis.csv')
analysis_df.to_csv(
        analysis_file, sep='\t', encoding='utf-8', index=False
)
print('Output to file ' + analysis_file)

# Plot callers
matplotlib.rcParams.update({'font.size': 16})
figure, axes = pyplot.subplots()
callers = [c.split('_')[-1] for c in utils.get_all_callers(df)]
at = analysis_df.transpose().reset_index()
at.rename(columns = at.iloc[0], inplace=True)
at = at[1:].reset_index(drop=True)
at_original = at.iloc[[
        i for i, caller in enumerate(at['ANALYSIS']) 
        if re.search('GT_', caller)
]] 
at_joint = at.iloc[[
        i for i, caller in enumerate(at['ANALYSIS']) 
        if re.search('JOINT_', caller)
]]

# Plot labels
pyplot.title('Mutation Caller ROC')
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.ylim(ymin=0)
pyplot.xlim(xmin=0)

for caller in caller_modules:
    if hasattr(caller, 'get_roc'):
        vals = caller.get_roc(testing) 
        axes.plot(
                # False Positive Rate
                [vals['FP'][i] / (vals['FP'][i] + vals['TN'][i])
                        for i in range(0, len(vals['TN']))
                ],
                # True Positive Rate
                [vals['TP'][i] / (vals['TP'][i] + vals['FN'][i])
                        for i in range(0, len(vals['TN']))
                ],
                '-', 
                label=caller.NAME
        )

# Original callers scatter plot
axes.scatter(
        at_original['False Positive Rate'],
        at_original['True Positive Rate'],
        marker='o', 
        color='b',
        s=50,
        label='Original callers'
)
# Joint callers scatter plot
axes.scatter(
        at_joint['False Positive Rate'],
        at_joint['True Positive Rate'],
        marker='^',
        color='g',
        s=100,
        label='Joint callers'
)

pyplot.legend(loc=4)

# Point labels
r = re.compile('^(GT_|JOINT_)')
callers = [c.split('_')[-1] for c in filter(r.match, list(analysis_df))]
annotations = []
for x, y, caller in zip(
        at['False Positive Rate'], at['True Positive Rate'], callers
):
    annotations.append(pyplot.text(x, y, caller))
adjust_text(annotations, arrowprops=dict(arrowstyle="-", color='r', lw=0.5))

pyplot.show()
