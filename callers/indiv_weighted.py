import numpy
import pandas

import utils

# All joint callers must have a NAME
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

# Returns a list of weights for each variant based on the number of callers
def get_vals(df, weights):
    callers = utils.get_og_callers(df)
    return [sum([weights[caller] for caller in callers 
                    if df[caller][i][1] == 'P'
            ]) for i in range(0, df.shape[0])
    ]

# Run the caller on all variants
def get_calls(df, vals, cutoff, true_str, false_str):
    callers = utils.get_og_callers(df)
    return [true_str if vals[i] > cutoff else false_str 
            for i in range(0, df.shape[0])
    ]

# List of cutoffs from 0 to the maximum possible, in increments of 0.01
def get_cutoffs(df):
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(df)
    cutoffs = numpy.arange(0, sum([weights[c] for c in callers]), 0.01)
    return cutoffs

# Joint callers which can be tuned by adjusting a cutoff value should have
# a get_roc function used for plotting an ROC curve
#######
# Returns true positive and false positive variants based on individual caller
# probabilities, with the cutoff as the dependent variable
def get_roc(df):
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(df)
    cutoffs = get_cutoffs(df)
    vals = get_vals(df, weights)
    calls = [get_calls(df, vals, cutoff, 'P', 'N') for cutoff in cutoffs]
    # Determine TP / TN / FP / FN
    statuses = [[
            str((call[i] == 'P') == df['REPORTABLE'][i])[0] + call[i]
            for i in range(0, df.shape[0])
    ] for call in calls]
    tp = [len([s for s in status if s == 'TP']) for status in statuses]
    fp = [len([s for s in status if s == 'FP']) for status in statuses]
    tn = [len([s for s in status if s == 'TN']) for status in statuses]
    fn = [len([s for s in status if s == 'FN']) for status in statuses]
    return {'TP': tp, 'FP': fp, 'TN': tn, 'FN': fn}

# All joint callers MUST have an add_caller function
######
# Add a caller based on the probability of an individual caller correctly
# finding a reportable variant
def add_caller(df, training):
    callers = utils.get_og_callers(df)
    weights = get_indiv_weights(training)
    cutoffs = get_cutoffs(df)
    vals = get_roc(training)
    # Get list of all possible Matthews Correlation Coefficients
    mccs = [
            utils.get_mcc(
                    vals['TP'][i], vals['TN'][i], vals['FP'][i], vals['FN'][i]
            ) for i in range(0, len(vals['TP']))
    ]
    # Get cutoff with the highest MCC
    cutoff = cutoffs[[i for i, mcc in enumerate(mccs) if mcc == max(mccs)][-1]]
    print('Individual cutoff: ' + str(cutoff))
    # Add caller
    weight_vals = get_vals(df, weights)
    df[NAME] = get_calls(df, weight_vals, cutoff, True, './.')
