import os
import sys
import glob
import importlib
import re
import math
import cProfile

import pandas
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
from adjustText import adjust_text

import utils
import parser

#pr = cProfile.Profile()
#pr.enable()

# Change these if needed
#db_name = '2017-06-30_NgsReviewer_master'
db_name = 'NgsReviewer'
hostname = '172.17.0.9'

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

# Import joint caller modules
print('Importing joint caller modules...')
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
    
# Read / generate combined data file
if os.path.isfile(df_file):
    print('Reading parsed samples file...')
    # Read parsed CSV
    df = pandas.read_csv(df_file, sep='\t', low_memory=False)
    # Convert SAMPLE_IDs to Strings
    df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
else:
    df = parser.combine(db_name, hostname, bed_file, samples_dir)

# Get indices of covered variants
print('Splitting training and testing sets...')
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
    print('Adding ' + caller.NAME + '...')
    caller.add_caller(df, training)

# Analyze combined callers
print('Classifying combined calls...')
parser.classify(df, utils.get_joint_callers(df))

# Get testing set
testing_indices = covered_indices[:covered_split]
testing = df.iloc[testing_indices].reset_index(drop=True)

# Analyze callers
print('Analyzing callers...')
# Initialize analysis DataFrame
analysis_df = pandas.DataFrame({
        'ANALYSIS' : ['True Positives', 'True Negatives',
        'False Positives', 'False Negatives', 'True Positive Rate',
        'True Negative Rate', 'Positive Predictive Value',
        'Negative Predictive Value', 'False Negative Rate',
        'False Positive Rate', 'False Discovery Rate',
        'False Omission Rate',
        'Accuracy', 'Matthews Correlation Coefficient']
})
callers = utils.get_all_callers(df)
for caller in callers:
    # Get true and false positives and negatives
    tp_df = utils.get_tp(testing, caller)
    fp_df = utils.get_fp(testing, caller)
    tn_df = utils.get_tn(testing, caller)
    fn_df = utils.get_fn(testing, caller)

    # Count rows in DataFrames
    tp = tp_df.shape[0]
    tn = tn_df.shape[0]
    fp = fp_df.shape[0]
    fn = fn_df.shape[0]

    # Analysis
    tpr = tp / (tp + fn) # True positive rate (sensitivity)
    tnr = tn / (tn + fp) # True negative rate (specificity)
    ppv = tp / (tp + fp) # Positive predictive value (precision)
    npv = tn / (tn + fn) # Negative predictive value
    fnr = fn / (fn + tp) # False negative rate (miss rate)
    fpr = fp / (fp + tn) # False positive rate (fall-out)
    fdr = fp / (fp + tp) # False discovery rate
    fom = fn / (fn + tn) # False omission rate
    acc = (tp + tn) / (tp + tn + fp + fn) # Accuracy
    # Matthews correlation coefficient
    mcc = ((tp * tn - fp * fn)
            / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))

    # Add column to analysis dataframe
    analysis_df[caller] = [
            tp, tn, fp, fn, tpr, tnr, ppv, npv, fnr, fpr, fdr, fom, acc, mcc
    ]
analysis_file = os.path.join(output_dir, 'analysis.csv')
analysis_df.to_csv(
        analysis_file, sep='\t', encoding='utf-8', index=False
)
print('Output to file ' + analysis_file)

# Plot callers
print('Plotting callers...')
# Transpose analysis DataFrame
at = analysis_df.transpose().reset_index()
at.rename(columns = at.iloc[0], inplace=True)
at = at[1:].reset_index(drop=True)
# DataFrame of original callers
at_original = at.iloc[[
        i for i, caller in enumerate(at['ANALYSIS']) 
        if re.search('GT_', caller)
]] 
# Joint callers
at_joint = at.iloc[[
        i for i, caller in enumerate(at['ANALYSIS']) 
        if re.search('JOINT_', caller)
]]

# Plot params
matplotlib.rcParams.update({'font.size': 10})
figure, axes = pyplot.subplots()
pyplot.title('Mutation Caller ROC')
pyplot.xlabel('False Positive Rate')
pyplot.ylabel('True Positive Rate')
pyplot.ylim(ymin=0)
pyplot.xlim(xmin=0)

# Plot ROC curves if applicable
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
# Plot original callers
axes.scatter(
        at_original['False Positive Rate'],
        at_original['True Positive Rate'],
        marker='o', 
        color='b',
        s=50,
        label='Original callers'
)
# Plot joint callers
axes.scatter(
        at_joint['False Positive Rate'],
        at_joint['True Positive Rate'],
        marker='^',
        color='g',
        s=100,
        label='Joint callers'
)

# Create legend
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

#pyplot.show()
figure.set_size_inches(10, 10)
pyplot.savefig(os.path.join(output_dir, 'plot.png'), bbox_inches='tight')

#pr.disable()
#pr.print_stats(sort='time')
