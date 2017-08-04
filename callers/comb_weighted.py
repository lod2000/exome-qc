import numpy
import pandas

import utils

# Counts only variants marked 'positive' by specific callers and 'negative' by
# all other callers.
def combine_callers(df, select_callers):
    all_callers = utils.get_og_callers(df)
    tps = df
    fps = df
    for caller in all_callers:
        if caller in select_callers:
            tps = utils.get_tp(tps, caller)
            fps = utils.get_fp(fps, caller)
        else:
            tps = utils.get_fn(tps, caller)
            fps = utils.get_tn(fps, caller)
    return {'TP': tps.shape[0], 'FP': fps.shape[0]}

# Returns a dict of probabilities that a variant is real for all possible
# combinations of callers. Gives a value 0 to any combination not found in the
# data. Example:
# GT_VARDICT            : 0.2
# GT_VARDICT&GT_VARSCAN : 0.8
# GT_GATK&GT_VARDICT    : 0.0
# (That the callers not listed didn't call the variant is implied)
def get_comb_weights(df):
    reportables = len([rep for rep in df['REPORTABLE'] if rep])
    covered = len([
            i for i in range(0, df.shape[0])
            if df['REPORTABLE'][i] or df['COVERED'][i]
    ])
    callers = utils.get_og_callers(df)
    weights = {}
    for k in range(1, len(callers) + 1):
        for subset in itertools.combinations(callers, k):
            sublist = list(subset)
            comb = combine_callers(df, sublist)
            called = comb['TP'] + comb['FP']
            if called > 0:
                weights['&'.join(sublist)] = utils.p_real_given_called(
                        comb['TP'], comb['FP'], reportables, covered
                )
            else:
                weights['&'.join(sublist)] = 0.0
    weights[''] = 0.0
    return weights

def get_calls(df, weights, cutoff, true_str, false_str):
    return [true_str if probs[
                    '&'.join([c for c in callers if df[c][i][1] == 'P'])
            ] > cutoff else false_str for i in range(0, df.shape[0])
    ]

# Function to return the number of true positives and false positives based on
# the probability of a particular combination of callers, with the true/false
# cutoff as the dependent variable
def get_roc(cutoffs, df):
    callers = utils.get_og_callers(df)
    probs = get_comb_weights(df)
    calls = [get_calls(df, weights, cutoff, 'P', 'N') for cutoff in cutoffs]
    #nps = [[
    #        'P'
    #        if probs['&'.join([c for c in callers if df[c][i][1] == 'P'])] > cutoff
    #        else 'N' for i in range(0, df.shape[0])
    #] for cutoff in cutoffs]
    statuses = [[
            str((np[i] == 'P') == df['REPORTABLE'][i])[0] + np[i]
            for i in range(0, df.shape[0])
    ] for np in nps]
    tp = [len([s for s in status if s == 'TP']) for status in statuses]
    fp = [len([s for s in status if s == 'FP']) for status in statuses]
    tn = [len([s for s in status if s == 'TN']) for status in statuses]
    fn = [len([s for s in status if s == 'FN']) for status in statuses]
    return {'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn}

# Add a caller based on the probability of a particular combination of callers
def add_prob_caller(df, training):
    name = 'JOINT_COMB'
    print('Adding ' + name + '...')
    callers = utils.get_og_callers(df)
    weights = get_comb_weights(training)
    cutoffs = numpy.arange(0, 1, 0.01)
    vals = get_roc(training, weights, cutoffs)
    df['JOINT_COMB'] = [
            True 
            if probs['&'.join([c for c in callers if df[c][i][1] == 'P'])] > cutoff
            else './.' for i in range(0, df.shape[0])
    ]
