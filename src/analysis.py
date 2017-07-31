import math
import itertools
import os
import re

import pandas
import numpy
from scipy.optimize import fmin
from scipy.optimize import minimize
import matplotlib.pyplot as pyplot
from adjustText import adjust_text

import parser

def get_true_positives(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'TP'
    ]].reset_index(drop=True)

def get_true_negatives(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'TN'
    ]].reset_index(drop=True)

def get_false_positives(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'FP'
    ]].reset_index(drop=True)

def get_false_negatives(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'FN'
    ]].reset_index(drop=True)

def get_unclassified(df, caller_name):
    return df.iloc[[
            i for i, call in enumerate(df[caller])
            if call == 'UP' or call == 'UN'
    ]].reset_index(drop=True)

# Counts only variants marked 'positive' by specific callers and 'negative' by
# all other callers.
def comb_callers(df, select_callers):
    all_callers = parser.get_og_caller_names(df)
    tps = df
    fps = df
    for caller in all_callers:
        if caller in select_callers:
            tps = get_true_positives(tps, caller)
            fps = get_false_positives(fps, caller)
        else:
            tps = get_false_negatives(tps, caller)
            fps = get_true_negatives(fps, caller)
    return {'TP': tps.shape[0], 'FP': fps.shape[0]}

# The probability that a variant is real given that it has been called.
# Bayes' rule: P(A|B) = P(B|A) * P(A) / P(B)
def p_real_given_called(tp, fp, reportables, covered):
    p_called_given_real = tp / reportables
    p_real = reportables / covered
    p_called = (tp + fp) / covered
    p_real_given_called = p_called_given_real * p_real / p_called
    return p_real_given_called

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
    callers = parser.get_og_caller_names(df)
    weights = {}
    for k in range(1, len(callers) + 1):
        for subset in itertools.combinations(callers, k):
            sublist = list(subset)
            comb = comb_callers(df, sublist)
            called = comb['TP'] + comb['FP']
            if called > 0:
                weights['&'.join(sublist)] = p_real_given_called(
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
    callers = parser.get_og_caller_names(df)
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
    return {'TP': tp, 'FP': fp}

# Add a caller based on the probability of a particular combination of callers
def add_prob_caller(df):
    print('adding prob caller')
    cutoff = 0.3
    callers = parser.get_og_caller_names(df)
    print('getting probs')
    probs = get_caller_comb_probs(df)
    print('adding column')
    df['COMB_PROB'] = [
            True 
            if probs['&'.join([c for c in callers if df[c][i][1] == 'P'])] > cutoff
            else './.' for i in range(0, df.shape[0])
    ]

# Returns the probabilities of individual callers predicting a reportable variant
def get_caller_weights(df):
    callers = parser.get_og_caller_names(df)
    weights = []
    for caller in callers:
        tp = get_true_positives(df, caller).shape[0]
        fp = get_false_positives(df, caller).shape[0]
        tn = get_true_negatives(df, caller).shape[0]
        fn = get_false_negatives(df, caller).shape[0]
        all_calls = tp + fp + tn + fn
        weights.append(p_real_given_called(tp, fp, tp + fn, all_calls))
    return weights

# Returns true positive and false positive variants based on individual caller
# probabilities, with the cutoff as the dependent variable
def weight_fn(cutoffs, cov, weights):
    callers = parser.get_og_caller_names(cov)
    nps = [[
            'P' if sum([
                    weight for k, weight in enumerate(weights) 
                    if cov[callers[k]][i][1] == 'P'
            ]) > cutoff else 'N' for i in range(0, cov.shape[0])
    ] for cutoff in cutoffs]
    statuses = [[
            str((np[i] == 'P') == cov['REPORTABLE'][i])[0] + np[i]
            for i in range(0, cov.shape[0])
    ] for np in nps]
    tp = [len([s for s in status if s == 'TP']) for status in statuses]
    fp = [len([s for s in status if s == 'FP']) for status in statuses]
    return {'TP': tp, 'FP': fp}

# Add a caller based on the probability of an individual caller correctly
# finding a reportable variant
def add_weight_caller(df):
    weights = get_caller_weights(df)
    callers = parser.get_og_caller_names(df)
    cutoff = 0.3
    df['COMB_WEIGHT' + str(cutoff)[-1]] = [
            True if sum([
                    weight for k, weight in enumerate(weights) 
                    if df[callers[k]][i][1] == 'P'
            ]) > cutoff else './.' for i in range(0, df.shape[0])
    ]

# Add callers which detect variants called by at least x original callers
def add_x_or_more(df):
    for cutoff in range(2, len(parser.get_og_caller_names(df))):
        name = 'COMB_' + str(cutoff) + 'ORMORE'
        df[name] = [
                True if callers >= cutoff else './.'
                for callers in df['TOTAL_CALLERS']
        ]

# Plot false positives vs true positives for all callers
def plot_callers(df, analysis_df, plots_dir, combined=True):
    caller_pref = '^GT_'
    if combined:
        caller_pref = '^(GT_|COMB_)' 
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
            if re.search('COMB_[A-Z]', caller)
    ]]

    weights = get_caller_weights(df)
    cov = df.iloc[[
            i for i in range(0, df.shape[0])
            if df['COVERED'][i] or df['REPORTABLE'][i]
    ]].reset_index(drop=True)
    c = numpy.arange(0, sum(weights), 0.1)
    c2 = numpy.arange(0, 1, 0.01)
    f = weight_fn(c, cov, weights)
    f2 = prob_fn(c2, cov)
    pyplot.plot(f['FP'], f['TP'], '--', color='y')
    pyplot.plot(f2['FP'], f2['TP'], '-', color='m')

    # Create scatter plot
    pyplot.scatter(
            at_gt['False Positives'], 
            at_gt['True Positives'], 
            marker='o', 
            color='b'
    )
    pyplot.scatter(
            at_ormore['False Positives'],
            at_ormore['True Positives'],
            marker='^',
            color='g'
    )
    pyplot.scatter(
            at_other['False Positives'],
            at_other['True Positives'],
            marker='*',
            color='r'
    )
    # Plot labels
    pyplot.title('Mutation caller positive hits')
    pyplot.ylabel('True positives')
    pyplot.xlabel('False positives')
    # Point labels
    annotations = []
    for x, y, caller in zip(
            at['False Positives'], at['True Positives'], callers
    ):
        annotations.append(pyplot.text(x, y, caller))
    adjust_text(annotations)

    #pyplot.axes([0, 0, 1, 1])
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmin=0)

    #pyplot.axes([0.6, 0.2, 0.8, 0.6])
    #pyplot.xlim(xmin=300, xmax=420)
    #pyplot.ylim(ymin=100, ymax=300)

    pyplot.savefig(os.path.join(plots_dir, 'callers.pdf'), format='pdf')
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
    callers = parser.get_caller_names(df)
    for caller in callers:
        # True positives DataFrame
        tp_df = get_true_positives(df, caller)
        # False positives DataFrame
        fp_df = get_false_positives(df, caller)
        # True negatives in DataFrame
        tn_df = get_true_negatives(df, caller)
        # False negatives DataFrame
        fn_df = get_false_negatives(df, caller)

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

