import pandas

# Reportable
def get_true_positives(df, caller_name):
    return df.iloc[[i for i, reportable in enumerate(df['REPORTABLE'])\
            if reportable and not df[caller_name][i] == './.']].reset_index(drop=True)

# Covered, not reportable
# Sometimes a caller will find multiple mutations in the same location. These
# are counted as separate false positives
def get_false_positives(df, caller_name):
    return df.iloc[[i for i, covered in enumerate(df['COVERED'])\
            if covered and not df['REPORTABLE'][i]\
            and not df[caller_name][i] == './.']].reset_index(drop=True)

def get_false_negatives(tp, gt):
    return gt.iloc[[i for i in range(0, gt.shape[0])
            if not (gt['GENE'][i] in tp['GENE'].tolist()
            and gt['DNA_CHANGE'][i] in tp['DNA_CHANGE'].tolist()
            and gt['PROTEIN_CHANGE'][i] in tp['PROTEIN_CHANGE'].tolist()
            and gt['SAMPLE_ID'][i] in tp['SAMPLE_ID'].tolist())]].reset_index(drop=True)

# Not covered, not reportable
def get_unclassified(df, caller_name):
    return df.iloc[[i for i, cov in enumerate(df['COVERED'])\
            if not cov and not df['REPORTABLE'][i]\
            and not df[caller_name][i] == './.']].reset_index(drop=True)

# TODO does the 'end' position include the last covered position?
# Returns a DataFrame of individual positions, genes, and chromosomes covered
def get_positions(panel):
    positions = []
    genes = []
    chromosomes = []
    for i in range(0, panel.shape[0]):
        pos_list = list(range(panel['START'][i], panel['END'][i]))
        positions += pos_list
        genes += len(pos_list) * [panel['GENE'][i]]
        chromosomes += len(pos_list) * [panel['CHROMOSOME'][i]]
    return pandas.DataFrame({ 'POSITION' : positions,\
                              'GENE' : genes,\
                              'CHROMOSOME' : chromosomes }).reset_index(drop=True)

# Returns list of positions which were not called (correctly) and were covered
# by the small panel
def get_true_negatives(fp, caller_name, positions):
    covered = fp.iloc[[i for i, covered in enumerate(fp['COVERED']) if covered\
            and not fp[caller_name][i] == './.']].reset_index(drop=True)
    indices = []
    for sample in list(set(covered['SAMPLE_ID'])):
        sample_covered = covered.iloc[[i for i, sample_id in enumerate(covered['SAMPLE_ID']) if sample_id == sample]].reset_index(drop=True)
        for p in range(0, positions.shape[0]):
            if int(positions['POSITION'][p]) in covered['POSITION'].astype(int).tolist():
                for c in range(0, covered.shape[0]):
                    if int(positions['POSITION'][p]) == int(covered['POSITION'][c]) and str(positions['GENE'][p]) == str(covered['GENE'][c]) and str(positions['CHROMOSOME'][p]) == str(covered['CHROMOSOME'][c]):
                        indices.append(p)
    return positions.drop(positions.index[indices])
    # return positions.iloc[[i for i, pos in enumerate(positions['POSITION']) if not (pos in covered['POSITION'].tolist() and positions['CHROMOSOME'][i] in covered['CHROMOSOME'] and positions['GENE'][i] in covered['GENE'])]].reset_index(drop=True)
