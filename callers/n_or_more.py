import pandas

import utils

# Add callers which detect variants called by at least x original callers
def add_caller(df, training):
    lower = 2
    upper = len(utils.get_og_caller(df)) + 1
    for cutoff in range(lower, upper):
        name = 'JOINT_' + str(cutoff) + 'ORMORE'
        print('Adding ' + name + '...')
        df[name] = [
                True if callers >= cutoff else './.'
                for callers in df['TOTAL_CALLERS']
        ]
