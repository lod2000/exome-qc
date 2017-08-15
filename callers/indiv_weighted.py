import numpy
import pandas

import utils

NAME = 'JOINT_INDIV'

# Returns the probabilities of individual callers predicting a reportable variant
def get_indiv_weights(df):
    callers = utils.get_og_callers(df)
    weights = {}
    for caller in callers:
        tp = utils.get_tp(df, caller).shape[0]
        fp = utils.get_fp(df, caller).shape[0]
        tn = utils.get_tn(df, caller).shape[0]
        fn = utils.get_fn(df, caller).shape[0]
        all_calls = tp + fp + tn + fn
        weights[caller] = utils.p_real_given_called(tp, fp, tp + fn, all_calls)
    return weights

def get_calls(df, weights, cutoff, true_str, false_str):
    callers = utils.get_og_callers(df)
    return [true_str if sum([
                    weights[caller] for caller in callers
                    if df[caller][i][1] == 'P'
            ]) > cutoff else false_str for i in range(0, df.shape[0])
    ]

def get_cutoffs(df):
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(df)
    cutoffs = numpy.arange(0, sum([weights[c] for c in callers]), 0.01)
    return cutoffs

# Returns true positive and false positive variants based on individual caller
# probabilities, with the cutoff as the dependent variable
def get_roc(df):
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(df)
    cutoffs = get_cutoffs(df)
    calls = [get_calls(df, weights, cutoff, 'P', 'N') for cutoff in cutoffs]
    statuses = [[
            str((call[i] == 'P') == df['REPORTABLE'][i])[0] + call[i]
            for i in range(0, df.shape[0])
    ] for call in calls]
    tp = [len([s for s in status if s == 'TP']) for status in statuses]
    fp = [len([s for s in status if s == 'FP']) for status in statuses]
    tn = [len([s for s in status if s == 'TN']) for status in statuses]
    fn = [len([s for s in status if s == 'FN']) for status in statuses]
    return {'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn}

# Add a caller based on the probability of an individual caller correctly
# finding a reportable variant
def add_caller(df, training):
    #name = 'JOINT_INDIV'
    print('Adding ' + NAME + '...')
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(training)
    cutoffs = get_cutoffs(df)
    vals = get_roc(training)
    mccs = [
            utils.get_mcc(
                    vals['TP'][i], vals['TN'][i], vals['FP'][i], vals['FN'][i]
            ) for i in range(0, len(vals['TP']))
    ]
    cutoff = cutoffs[[i for i, mcc in enumerate(mccs) if mcc == max(mccs)][-1]]
    print('Individual cutoff: ' + str(cutoff))
    df[NAME] = get_calls(df, weights, cutoff, True, './.')
