import pandas

import parse_samples

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
    positions = parse_samples.split_panel(panel)
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
    # return positions.iloc[[i for i, pos in enumerate(positions['POSITION']) if not (pos in covered['POSITION'].tolist() and positions['CHROMOSOME'][i] in covered['CHROMOSOME'] and positions['GENE'][i] in covered['GENE'])]].reset_index(drop=True)
