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

if not os.path.isfile(df_file):
    sample_parser.combine(gt_file, bed_file, samples_dir)

df = pandas.read_csv(df_file, sep='\t', low_memory=False)#, dtype=object)
        # dtype={'CHROMOSOME': str, 'POSITION': int, 'TOTAL_CALLERS': int,\
        # 'GENE': str, 'DNA_CHANGE': str, 'PROTEIN_CHANGE': str,\
        # 'EXAC_FREQ_ESTIMATE': str, '1KG_FREQ_ESTIMATE': str,\
        # 'GT_VARSCAN': str, 'GT_GATK': str, 'GT_SAMTOOLS': str,\
        # 'GT_FREEBAYES': str, 'GT_VARDICT': str, 'SMPLE_ID': str,\
        # 'COVERED': bool, 'REPORTABLE': bool})
gt_parsed = sample_parser.parse_gt(gt_file)
match_list = sample_parser.find_matches(next(os.walk(samples_dir))[1], gt_parsed)
gt = pandas.concat(sample_parser.split_gt(gt_parsed, match_list)).reset_index(drop=True)

# DataFrame terminal display options
pandas.set_option('display.max_columns', 14)
pandas.set_option('display.max_rows', 80)
pandas.set_option('display.width', 100)

for caller in sample_parser.get_caller_names(df):
    print(caller)
    # print('True positives:')
    tp = analysis.get_true_positives(df, caller)
    fn = analysis.get_false_negatives(tp, gt)
    # print(tp)
    # print('False negatives:')
    # print(analysis.get_false_negatives(tp, gt))
    print('TP+FN: ' + str(tp.shape[0]+fn.shape[0]))
