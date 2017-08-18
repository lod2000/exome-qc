import math
import re

import pandas

# This file contains functions which are used in multiple other files

# True positives
def get_tp(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'TP'
    ]].reset_index(drop=True)

# True negatives
def get_tn(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'TN'
    ]].reset_index(drop=True)

# False positives
def get_fp(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'FP'
    ]].reset_index(drop=True)

# False negatives
def get_fn(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call == 'FN'
    ]].reset_index(drop=True)

def get_unclassified(df, caller):
    return df.iloc[[
            i for i, call in enumerate(df[caller]) if call[0] == 'U'
    ]].reset_index(drop=True)

def get_accuracy(tp, tn, fp, fn):
    return (tp + tn) / (tp + tn + fp + fn)

# Matthews Correlation Coefficient
def get_mcc(tp, tn, fp, fn):
    return ((tp * tn - fp * fn)
                / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))

# The probability that a variant is real given that it has been called.
# Bayes' rule: P(A|B) = P(B|A) * P(A) / P(B)
def p_real_given_called(tp, fp, reportables, covered):
    p_called_given_real = tp / reportables
    p_real = reportables / covered
    p_called = (tp + fp) / covered
    p_real_given_called = p_called_given_real * p_real / p_called
    return p_real_given_called

# Returns the names of the variant callers
def get_og_callers(df):
    headers = list(df)
    return [h for h in headers if re.search('^GT_', h)]

def get_joint_callers(df):
    headers = list(df)
    return [h for h in headers if re.search('^JOINT_', h)]

def get_all_callers(sample):
    headers = list(sample)
    return [h for h in headers if re.search('^(GT_|JOINT_)', h)]

def query_yes_no(question, default='yes'):
    """
    Ask a yes/no question and return their answer as a boolean
    'question' is a string presented to the user
    'default' determines what to do if the user just hits <Enter>
    """
    valid = {
            'yes': True, 'ye': True, 'y': True, 
            'no': False, 'n': False, 
    }
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('Invalid default answer: "%s"' % default)

    while True:
        choice = input(question + prompt).lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            print('Please respond with yes/y or no/n')
