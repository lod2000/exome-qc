import math

import pandas
import numpy
from scipy.optimize import fmin
import matplotlib.pyplot as pyplot

import parser

# Returns DataFrame of true positives (reportables)
def get_true_positives(df, caller):
#    return df.iloc[[
#            i for i, reportable in enumerate(df['REPORTABLE'])
#            if reportable and not df[caller_name][i] == './.'
#    ]].reset_index(drop=True)
    return df.iloc[[i for i, call in enumerate(df[caller]) if call == 'TP']].reset_index(drop=True)

# Covered, not reportable
# Sometimes a caller will find multiple mutations in the same location. These
# are counted as separate false positives
def get_false_positives(df, caller):
#    return df.iloc[[
#            i for i, covered in enumerate(df['COVERED'])
#            if covered and not df['REPORTABLE'][i]
#            and not df[caller][i] == './.'
#    ]].reset_index(drop=True)
    return df.iloc[[i for i, call in enumerate(df[caller]) if call == 'FP']].reset_index(drop=True)

# Returns DataFrame from ground truth of variants not detected by caller
def get_false_negatives(df, caller):
#    return gt.iloc[[
#            i for i in range(0, gt.shape[0])
#            if not (gt['GENE'][i] in tp['GENE'].tolist()
#            and gt['DNA_CHANGE'][i] in tp['DNA_CHANGE'].tolist()
#            and gt['PROTEIN_CHANGE'][i] in tp['PROTEIN_CHANGE'].tolist()
#            and gt['SAMPLE_ID'][i] in tp['SAMPLE_ID'].tolist())
#    ]].reset_index(drop=True)
    return df.iloc[[i for i, call in enumerate(df[caller]) if call == 'FN']].reset_index(drop=True)

# Not covered, not reportable
def get_unclassified(df, caller_name):
    return df.iloc[[
            i for i, cov in enumerate(df['COVERED'])
            if not cov and not df['REPORTABLE'][i]
            and not df[caller_name][i] == './.'
    ]].reset_index(drop=True)

# Returns list of positions which were not called (correctly) and were covered
# by the small panel
# TODO speed up
def get_true_negatives(fp, caller_name, panel):
    positions = parser.split_panel(panel)
    all_covered = fp.iloc[[
            i for i, covered in enumerate(fp['COVERED'])
            if covered and not fp[caller_name][i] == './.'
    ]].reset_index(drop=True)
    samples = list(set(fp['SAMPLE_ID']))
    # Each sample has the same covered positions
    all_positions = pandas.DataFrame(
            columns=['SAMPLE_ID', 'POSITION', 'CHROMOSOME', 'GENE']
    )
    for sample in samples:
        indices = []
        positions['SAMPLE_ID'] = [sample] * len(positions)
        covered = all_covered.iloc[[
                i for i, sample_id in enumerate(all_covered['SAMPLE_ID'])
                if sample_id == sample
        ]].reset_index(drop=True)
        called_positions = [
                i for i, pos in enumerate(positions['POSITION'])
                if pos in covered['POSITION'].tolist()
        ]
        for p in called_positions:
            for c in range(0, covered.shape[0]):
                if (positions['POSITION'][p] == covered['POSITION'][c]
                        and positions['GENE'][p] == covered['GENE'][c]
                        and positions['CHROMOSOME'][p] == covered['CHROMOSOME'][c]):
                    indices.append(p)
        all_positions = pandas.concat(
                [all_positions, positions.drop(positions.index[indices])],
                ignore_index=True
        )
    return all_positions

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

    for caller in parser.get_caller_names(df):
        print('\n')
        print(caller)

        # True positives DataFrame
        tp_df = get_true_positives(df, caller)
        # False positives DataFrame
        fp_df = get_false_positives(df, caller)
        # True negatives in DataFrame
        tn_df = get_true_negatives(fp_df, caller, panel)
        # False negatives DataFrame
        fn_df = get_false_negatives(df, caller)

        # Count rows in DataFrames
        tp = tp_df.shape[0]
        print(tp)
        tn = tn_df.shape[0]
        print(tn)
        fp = fp_df.shape[0]
        print(fp)
        fn = fn_df.shape[0]
        print(fn)
        print(str(tp + fn))
        print(str(tn+fp))

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

def f(weight, caller, df):
    weighted_sum = 0
    for i in range(0, df.shape[0]):
        true_status = int(df['REPORTABLE'][i] == True)
        call = int(not df[caller][i] == './.')
        weighted_sum += (true_status - weight * call) ** 2
    return weighted_sum

# WIP
def generate_combined_caller_weights(df):
    callers = parser.get_og_caller_names(df)
    weights = []
    for k, caller in enumerate(callers):
        print(caller)
        weight = fmin(lambda weight: f(weight, caller, df), 0, maxfun=20)[0]
        weights.append(weight)
        # Test
        print(weight)
        print(f(weight, caller, df))
        print(f(1, caller, df))
    return weights

def add_combined_caller(df, weights):
    callers = parser.get_og_caller_names(df)
    combined_status = []
    for i in range(0, df.shape[0]):
        weighted_calls = []
        for k, caller in enumerate(callers):
            call = int(not df[caller][i] == './.')
            weight = weights['WEIGHT'][k]
            weighted_calls.append(weight * call)
        combined_status.append(sum(weighted_calls))
    df['COMBINED_STATUS'] = combined_status
    df['COMB_TF'] = [
            True if status >= 0.009 else './.' for status in df['COMBINED_STATUS']
    ]

# def add_x_or_more(df):
#     for cutoff in range(2, len(parser.get_og_caller_names(df)) + 1):
#         name = 'COMB_' + str(cutoff) + '_OR_MORE'
#         df[name] = [
#                 True if callers >= cutoff else './.'
#                 for callers in df['TOTAL_CALLERS']
#         ]

def add_2_or_more(df):
    name = 'COMB_2_OR_MORE'
    df[name] = [
            True if callers >= 2 else './.' for callers in df['TOTAL_CALLERS']
    ]

# WIP
def add_differences(df):
    callers = parser.get_og_caller_names(df)
    all_except = list(callers)
    for caller1 in callers:
        all_except.remove(caller1)
        print(caller1)
        print(all_except)
        # callers.remove(caller1)
        for caller2 in all_except:
            name = 'COMB_' + caller1.split('_')[-1] + '_U_' + caller2.split('_')[-1]
            df[name] = [
                    True if not df[caller1][i] == './.'
                    and not df[caller2][i] == './.'
                    else './.' for i in range(0, df.shape[0])
            ]

def plot_callers(analysis_df):
    callers = [c.split('_')[-1] for c in list(analysis_df.columns)[1:]]
    at = analysis_df.transpose().reset_index()
    at.rename(columns = at.iloc[0], inplace=True)
    at = at[1:].reset_index(drop=True)

    # Create scatter plot
    pyplot.scatter(
            at['False Positives'], 
            at['True Positives'], 
            marker='o', 
    )
    # Plot labels
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmin=0)
    pyplot.title('Mutation caller positive hits')
    pyplot.ylabel('True positives')
    pyplot.xlabel('False positives')
    # Point labels
    for caller, x, y in zip(
            callers, at['False Positives'], at['True Positives']
    ):
        pyplot.annotate(
                caller, xy=(x, y), 
                xytext=(-10, -10), 
                textcoords='offset points', 
                ha='right', 
                va='bottom'
        )
    pyplot.show()
