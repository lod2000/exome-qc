import pandas
import matplotlib.pyplot as pyplot

import sample_parser

# Reportable
def get_true_positives(df, caller_name):
    return df.iloc[[
            i for i, reportable in enumerate(df['REPORTABLE'])
            if reportable and not df[caller_name][i] == './.'
    ]].reset_index(drop=True)

# Covered, not reportable
# Sometimes a caller will find multiple mutations in the same location. These
# are counted as separate false positives
def get_false_positives(df, caller_name):
    return df.iloc[[
            i for i, covered in enumerate(df['COVERED'])
            if covered and not df['REPORTABLE'][i]
            and not df[caller_name][i] == './.'
    ]].reset_index(drop=True)

def get_false_negatives(tp, gt):
    return gt.iloc[[
            i for i in range(0, gt.shape[0])
            if not (gt['GENE'][i] in tp['GENE'].tolist()
            and gt['DNA_CHANGE'][i] in tp['DNA_CHANGE'].tolist()
            and gt['PROTEIN_CHANGE'][i] in tp['PROTEIN_CHANGE'].tolist()
            and gt['SAMPLE_ID'][i] in tp['SAMPLE_ID'].tolist())
    ]].reset_index(drop=True)

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
    positions = sample_parser.split_panel(panel)
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
