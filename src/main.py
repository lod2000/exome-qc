import os
import math
import sys
import glob

import pandas

import parse_samples
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
    parse_samples.combine(gt_file, bed_file, samples_dir)

df = pandas.read_csv(df_file, sep='\t', low_memory=False)
# Convert SAMPLE_IDs to Strings
df['SAMPLE_ID'] = df['SAMPLE_ID'].astype(str)
gt_parsed = parse_samples.parse_gt(gt_file)
match_list = parse_samples.find_matches(
        next(os.walk(samples_dir))[1], gt_parsed)
gt = pandas.concat(
        parse_samples.split_gt(gt_parsed, match_list)
).reset_index(drop=True)
panel = parse_samples.parse_bed(bed_file)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 14)
pandas.set_option('display.max_rows', 13)
pandas.set_option('display.width', 300)

# Initialize analysis DataFrame
analysis_df = pandas.DataFrame({
        'ANALYSIS' : ['True Positives', 'True Negatives',
        'False Positives', 'False Negatives', 'True Positive Rate',
        'True Negative Rate', 'Positive Predictive Value',
        'Negative Predictive Value', 'False Negative Rate',
        'False Positive Rate', 'False Discovery Rate', 'False Omission Rate',
        'Accuracy', 'Matthews Correlation Coefficient']
})

for caller in parse_samples.get_caller_names(df):
    print('\n')
    print(caller)

    # True positives DataFrame
    tp_df = analysis.get_true_positives(df, caller)
    # False positives DataFrame
    fp_df = analysis.get_false_positives(df, caller)
    # True negatives in DataFrame
    tn_df = analysis.get_true_negatives(fp_df, caller, panel)
    # False negatives DataFrame
    fn_df = analysis.get_false_negatives(tp_df, gt)

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

    analysis_df[caller] = [
            tp, tn, fp, fn, tpr, tnr, ppv, npv, fnr, fpr, fdr, fom, acc, mcc
    ]

analysis_file = os.path.join(data_path, 'analysis.csv')
try:
    analysis_df.to_csv(analysis_file, sep='\t', encoding='utf-8', index=False)
except PermissionError:
    print('Another program is using the file analysis.csv.\
           Please close and try again.')
