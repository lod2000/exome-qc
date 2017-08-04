import pandas

import parser

# Add callers which detect variants called by at least x original callers
def add_caller(df):
    lower = 2
    upper = len(parser.get_og_caller_names(df)) + 1
    for cutoff in range(lower, upper):
        name = 'JOINT_' + str(cutoff) + 'ORMORE'
        print('Adding ' + name + '...')
        df[name] = [
                True if callers >= cutoff else './.'
                for callers in df['TOTAL_CALLERS']
        ]
