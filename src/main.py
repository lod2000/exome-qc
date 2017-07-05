import os, pandas, sys, glob
import sample_parser, analysis

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
match_list = sample_parser.find_matches(next(os.walk(samples_dir))[1], gt_parsed)
gt = pandas.concat(sample_parser.split_gt(gt_parsed, match_list)).reset_index(drop=True)
panel = sample_parser.parse_bed(bed_file)
positions = analysis.get_positions(panel)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 14)
pandas.set_option('display.max_rows', 13)
pandas.set_option('display.width', 300)

print(positions.shape[0])

# Initialize analysis DataFrame
analysis_df = pandas.DataFrame({'ANALYSIS':['True Positives','False Positives',\
        'False Negatives','True Positive Rate','Positive Predictive Value',\
        'False Negative Rate','False Discovery Rate']})

for caller in sample_parser.get_caller_names(df):
    print('\n')
    print(caller)

    # True positives DataFrame
    tp = analysis.get_true_positives(df, caller)
    # False positives DataFrame
    fp = analysis.get_false_positives(df, caller)
    # False negatives DataFrame
    fn = analysis.get_false_negatives(tp, gt)
    # True negatives in DataFrame
    tn = analysis.get_true_negatives(fp, caller, positions)
    # tn.to_csv(os.path.join(samples_dir, caller + '_tn.tab'), sep='\t', encoding='utf-8', index=False)

    # Count rows in DataFrames
    tps = tp.shape[0]
    tns = tn.shape[0]
    fps = fp.shape[0]
    fns = fn.shape[0]

    # Analysis
    # tpr = tps / (tps + fns)
    # ppv = tps / (tps + fps)
    # fnr = fns / (fns + tps)
    # fdr = fps / (fps + tps)

    # print(tps+fns)
    print(fps+tns)
    # print(fp)
    # print(tns)
    # print(fp)
    # print(fns)
    # print(fp)

    # analysis_df[caller] = [tps, fps, fns, tpr, ppv, fnr, fdr]

# analysis_file = os.path.join(data_path, 'analysis.csv')
# try:
    # analysis_df.to_csv(analysis_file, sep='\t', encoding='utf-8', index=False)
# except PermissionError:
#     print('Another program is using the file analysis.csv. Please close and try again.')
