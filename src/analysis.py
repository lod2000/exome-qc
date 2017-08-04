import math
import itertools
import os
import re
import sys

import pandas
import numpy
import matplotlib.pyplot as pyplot
from adjustText import adjust_text
import matplotlib

import utils
sys.path.append(os.path.join(sys.path[0], '..', 'callers'))
import indiv_weighted

# Counts only variants marked 'positive' by specific callers and 'negative' by
# all other callers.
def comb_callers(df, select_callers):
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
def get_caller_comb_probs(df):
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
            comb = comb_callers(df, sublist)
            called = comb['TP'] + comb['FP']
            if called > 0:
                weights['&'.join(sublist)] = utils.p_real_given_called(
                        comb['TP'], comb['FP'], reportables, covered
                )
            else:
                weights['&'.join(sublist)] = 0.0
    weights[''] = 0.0
    return weights

# Function to return the number of true positives and false positives based on
# the probability of a particular combination of callers, with the true/false
# cutoff as the dependent variable
def prob_fn(cutoffs, df):
    callers = utils.get_og_callers(df)
    probs = get_caller_comb_probs(df)
    nps = [[
            'P'
            if probs['&'.join([c for c in callers if df[c][i][1] == 'P'])] > cutoff
            else 'N' for i in range(0, df.shape[0])
    ] for cutoff in cutoffs]
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
def add_prob_caller(cutoff, df):
    callers = utils.get_og_callers(df)
    probs = get_caller_comb_probs(df)
    df['JOINT_COMB'] = [
            True 
            if probs['&'.join([c for c in callers if df[c][i][1] == 'P'])] > cutoff
            else './.' for i in range(0, df.shape[0])
    ]

# Plot false positives vs true positives for all callers
def plot_callers(df, analysis_df, plots_dir, combined=True):
    matplotlib.rcParams.update({'font.size': 16})
    figure, axes = pyplot.subplots()
    caller_pref = '^GT_'
    if combined:
        caller_pref = '^(GT_|JOINT_)' 
    r = re.compile(caller_pref)
    callers = [c.split('_')[-1] for c in filter(r.match, list(analysis_df))]
    at = analysis_df.transpose().reset_index()
    at.rename(columns = at.iloc[0], inplace=True)
    at = at[1:].reset_index(drop=True)
    at_gt = at.iloc[[
            i for i, caller in enumerate(at['ANALYSIS']) 
            if re.search('GT_', caller)
    ]] 
    at_ormore = at.iloc[[
            i for i, caller in enumerate(at['ANALYSIS']) 
            if re.search('ORMORE', caller)
    ]]
    at_other = at.iloc[[
            i for i, caller in enumerate(at['ANALYSIS']) 
            if re.search('JOINT_[A-Z]', caller)
    ]]

    covered = df.iloc[[
            i for i in range(0, df.shape[0])
            if df['COVERED'][i] or df['REPORTABLE'][i]
    ]].reset_index(drop=True)

    # Plot labels
    pyplot.title('Mutation Caller ROC')
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')

    # Individual probability caller curve
    callers = utils.get_og_callers(df)
    indiv_weights = indiv_weighted.get_indiv_weights(df)
    indiv_cutoffs = numpy.arange(0, sum([indiv_weights[c] for c in callers]), 0.1)
    indiv_vals = indiv_weighted.get_roc(covered, indiv_weights, indiv_cutoffs)
    indiv_line = axes.plot(
            # False Positive Rate
            [indiv_vals['FP'][i] / (indiv_vals['FP'][i] + indiv_vals['TN'][i])
                    for i in range(0, len(indiv_vals['TN']))
            ],
            # True Positive Rate
            [indiv_vals['TP'][i] / (indiv_vals['TP'][i] + indiv_vals['FN'][i])
                    for i in range(0, len(indiv_vals['TN']))
            ],
            '--', 
            color='y', 
            label='Individual weighted caller'
    )

    # Combined probability caller curve
    comb_cutoffs = numpy.arange(0, 1, 0.01)
    comb_vals = prob_fn(comb_cutoffs, covered)
    axes.plot(
            # False Positive Rate
            [comb_vals['FP'][i] / (comb_vals['TN'][i] + comb_vals['FP'][i])
                    for i in range(0, len(comb_vals['TN']))
            ],
            # True Positive Rate
            [comb_vals['TP'][i] / (comb_vals['TP'][i] + comb_vals['FN'][i])
                    for i in range(0, len(comb_vals['TN']))
            ],
            '-', 
            color='m', 
            label='Combined weighted caller'
    )

    # Original callers scatter plot
    axes.scatter(
            at_gt['False Positive Rate'],
            at_gt['True Positive Rate'],
            marker='o', 
            color='b',
            s=50,
            label='Original callers'
    )
    # n or more callers scatter plot
    axes.scatter(
            at_ormore['False Positive Rate'],
            at_ormore['True Positive Rate'],
            marker='^',
            color='g',
            s=100,
            label='N or more'
    )
    # Probability callers scatter plot
    axes.scatter(
            at_other['False Positive Rate'],
            at_other['True Positive Rate'],
            marker='*',
            color='r',
            s=100,
            label='Weighted caller cutoff'
    )

    pyplot.ylim(ymin=0)
    pyplot.xlim(xmin=0)
    pyplot.legend(loc=4)
    # Point labels
    callers = [c.split('_')[-1] for c in filter(r.match, list(analysis_df))]
    annotations = []
    for x, y, caller in zip(
            at['False Positive Rate'], at['True Positive Rate'], callers
    ):
        annotations.append(pyplot.text(x, y, caller))
    adjust_text(annotations, arrowprops=dict(arrowstyle="-", color='r', lw=0.5))

    pyplot.show()

# Returns DataFrame of analysis
def analyze_callers(df, panel):
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

    print('Analyzing callers...')
    callers = utils.get_all_callers(df)
    for caller in callers:
        # True positives DataFrame
        tp_df = utils.get_tp(df, caller)
        # False positives DataFrame
        fp_df = utils.get_fp(df, caller)
        # True negatives in DataFrame
        tn_df = utils.get_tn(df, caller)
        # False negatives DataFrame
        fn_df = utils.get_fn(df, caller)

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
    
    return analysis_df

