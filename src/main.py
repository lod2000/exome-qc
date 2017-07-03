import argparse, os, pandas
import sample_parser, analysis

# Parser / command line instructions
# Description of script
arg_parser = argparse.ArgumentParser(description='Compares sample variants to\
        ground truth file and produces false negatives, false positives,\
        and true positives.')
# First argument: ground truth Excel file
arg_parser.add_argument('file1', help='ground truth .xlsx', action='store')
# Second argument: small panel gene list, bed file
arg_parser.add_argument('file2', help='.bed', action='store')
# Third argument: directory containing sample directories
arg_parser.add_argument('directory', help='directory containing all sample\
        directories', action='store')

# Parse arguments
args = arg_parser.parse_args()
gt_file = args.file1
bed_file = args.file2
samples_dir = args.directory

df_file = os.path.join(samples_dir, 'combined.tab')

if not os.path.isfile(df_file):
    sample_parser.combine(gt_file, bed_file, samples_dir)

df = pandas.read_csv(df_file, sep='\t')
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
    print('TP+FP: ' + str(tp.shape[0]+fn.shape[0]))
