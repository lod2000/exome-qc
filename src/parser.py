import re
import os
import time
import sys

import pandas
import pymongo
import numpy
from bson import ObjectId

import analysis

def replace_key(dictionary, new_key, old_key):
    dictionary[new_key] = dictionary[old_key]
    del dictionary[old_key]

# Get DataFrame of variants deemed "reportable" in final mongo database
# Also outputs to ground_truth.csv
def get_reportables(db_name):
    # Get mongo database
    # Requires mongorestore to have been run
    client = pymongo.MongoClient()
    db = client[db_name]
    # Mongo collections
    final = db.mongo_final.final_results
    torrent = db.mongo_torrent.torrent_results
    server = db.mongo_server.server_downloads
    # Generate cursor for reportable final documents with a locus
    reportable_cursor = final.find({
            'reportable': True, 
            'germline': False, 
            'ion_junk': False,
            '_locus': {'$ne': ':'}
    })
    combined = []
    print('Combining final and torrent databases...')
#    # Set up progress indicator
#    sys.stdout.write('Parsing mongo database:  0%')
#    sys.stdout.flush()
    for i, final_doc in enumerate(reportable_cursor):
#        # Progress indicator
#        progress = int(100 * (i / reportable_cursor.count()))
#        if progress < 10:
#            sys.stdout.write('\b\b' + str(progress) + '%')
#        else:
#            sys.stdout.write('\b\b\b' + str(progress) + '%')
#        sys.stdout.flush()
        torrent_doc = torrent.find_one(final_doc['torrent_result_id'])
        # Generate combined dict of torrent and final information
        try:
            # Replace potential duplicate key names
            replace_key(final_doc, 'final_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_result_id', '_id')
            replace_key(final_doc, 'final_result_id', '_id')
            # Remove unnecessary fields
            if 'snapshot' in torrent_doc:
                del torrent_doc['snapshot']
            del final_doc['ion_junk']
            del final_doc['reportable']
            del final_doc['germline']
            # Get information from server_downloads
            server_doc = server.find_one(torrent_doc['server_download_id'])
            if not server_doc == None:
                # Sample ID
                if 'sample_id' in server_doc:
                    torrent_doc['server_sample_id'] = server_doc['sample_id']
                # Sample annotation tab file path
                if 'ppmp_annotation_path' in server_doc:
                    torrent_doc['ppmp_annotation_path'] = server_doc['ppmp_annotation_path']
            # Combine torrent and final documents
            combined.append({**torrent_doc,  **final_doc})
        # In case there isn't a matching torrent result
        except TypeError:
            print('No torrent result for ' + str(final_doc['torrent_result_id']))
    # Generate combined DataFrame
    df = pandas.DataFrame(combined)
    return df

# Simplify reportable ground truth data to contain only relevant information
# Also outputs to ground_truth.csv
def get_simplified_gt(db_name, output_dir):
    df = get_reportables(db_name)
    print('Generating simplified data frame...')
    # Remove MiSeq entries
    df = df.iloc[[i for i, miseq in enumerate(df['miseq_run']) if not miseq]]    
    df.reset_index(drop=True, inplace=True)
    # Initialize lists
    sample_names = []
    chromosomes = []
    positions = []
    genes = []
    dna_changes = []
    protein_changes = []
    ppmp_annotation_paths = []
    target_allele_counts = []
    exac_freq_estimates = []
    nkg_freq_estimates = []
    headers = list(df)
    # Cycle through variants
    for i in range(0, df.shape[0]):
        # Get paths to sample annotation tab files
        if df['ppmp_annotation_path'].notnull()[i]:
            ppmp_annotation_paths.append(df['ppmp_annotation_path'][i].split('Results/')[-1])
        else:
            ppmp_annotation_paths.append('')
        sample_names.append(df['sample_name'][i])
        chromosomes.append(df['_locus'][i].split(':')[0])
        positions.append(df['_locus'][i].split(':')[1])
        target_allele_counts.append(df['target_allele_count'][i])
        exac_freq_estimates.append(df['exac_freq_estimate'][i])
        nkg_freq_estimates.append(df['nkg_freq_estimate'][i])
        # Try to get gene, dna, and protein data from the same column
        variant = str(df['reported_variant'][i])
        if re.search('^[A-Z0-9]*\sc.[^\s]*:', variant):
            genes.append(variant.split()[0])
            dna_changes.append(variant.split()[1].split(':')[0])
            try:
                protein_changes.append(re.search('p\..*$', variant).group(0))
            except AttributeError:
                protein_changes.append('')
        # Pull gene, dna, and protein data from different locations
        else:
            # Gene
            genes.append('')
            gene_columns = [h for h in headers if re.search('gene', h)]
            gene_columns += ['biomarker_name', 'protein_id']
            for column in gene_columns:
                if re.search('^[A-Z]', str(df[column][i])):
                    genes[i] = df[column][i]
                    break
            # DNA
            dna_changes.append('')
            dna_columns = [
                    'hgvs_coding_change', 'snp_eff_coding',
                    'snpeff_annotated_dna_change'
            ]
            for column in dna_columns:
                dna_search = re.search('c.[^\s]*$', str(df[column][i]))
                if dna_search:
                    dna_changes[i] = dna_search.group(0)
                    break
            # Protein
            protein_changes.append('')
            protein_columns = [
                    'hgvs_protein_change', 'snp_eff_protein',
                    'snpeff_annotated_protein_change'
            ]
            for column in protein_columns:
                protein_search = re.search('p.[^\s]*$', str(df[column][i]))
                if protein_search:
                    protein_changes[i] = protein_search.group(0)
                    break
    # Create DataFrame with only necessary information
    simple_df = pandas.DataFrame({
            'SAMPLE_NAME': sample_names, 'CHROMOSOME': chromosomes,
            'POSITION': positions, 'GENE': genes, 'DNA_CHANGE': dna_changes,
            'PROTEIN_CHANGE': protein_changes,
            'SAMPLE_PATH': ppmp_annotation_paths,
            'TARGET_ALLELE_COUNT': target_allele_counts,
            'EXAC_FREQ_ESTIMATE': exac_freq_estimates,
#            '1KG_FREQ_ESTIMATE': nkg_freq_estimates
    })
    # Strip out variants without DNA or protein information
    simple_df = simple_df.iloc[[
            i for i, dna in enumerate(simple_df['DNA_CHANGE']) 
            if not dna == '' and not simple_df['PROTEIN_CHANGE'][i] == ''
    ]].reset_index(drop=True)
    simple_df['SAMPLE_ID'] = [
            path.split(os.sep)[0] for path in simple_df['SAMPLE_PATH']
    ]
    # Output CSV
    simple_df.to_csv(
            os.path.join(output_dir, 'ground_truth.csv'),
            sep='\t', encoding='utf8', index=False
    )
    return simple_df

def parse_bed(bed):
    parsed = pandas.read_csv(
            bed, names=['CHROMOSOME', 'START', 'END', 'GENE'], sep='\t'
    )
    genes = [
            re.search('_[A-Z0-9]*_', gene_string).group(0).split('_')[1]
            for gene_string in parsed['GENE']
    ]
    parsed['GENE'] = genes
    return parsed

# TODO does the 'end' position include the last covered position?
# Returns a DataFrame of individual positions, genes, and chromosomes covered
def split_panel(panel):
    positions = []
    genes = []
    chromosomes = []
    for i in range(0, panel.shape[0]):
        pos_list = list(range(panel['START'][i], panel['END'][i]))
        positions += pos_list
        genes += len(pos_list) * [panel['GENE'][i]]
        chromosomes += len(pos_list) * [panel['CHROMOSOME'][i]]
    return pandas.DataFrame({
            'POSITION' : positions,
            'GENE' : genes,
            'CHROMOSOME' : chromosomes
    }).reset_index(drop=True)

# Returns the names of the variant callers
def get_og_caller_names(sample):
    headers = list(sample)
    return [h for h in headers if re.search('^GT_', h)]

def get_caller_names(sample):
    headers = list(sample)
    return [h for h in headers if re.search('^GT_', h) or re.search('^COMB_', h)]

# Returns list of sample paths from ground truth that are in data/samples
def find_samples(gt, samples_dir):
    samples = next(os.walk(samples_dir))[1]
    potentials = set(gt['SAMPLE_PATH'])
    finals = []
    for potential in potentials:
        if os.path.isfile(os.path.join(samples_dir, str(potential))):
            finals.append(os.path.join(*potential.split('/')))            
    return finals

# Returns a DataFrames of sample for which matches exist in the ground truth
def combine_samples(samples_dir, sample_paths):
    # List of sample DataFrames
    df_list = []
    for sample_path in sample_paths:
        sample_id = sample_path.split(os.sep)[0]
        tab_path = os.path.join(samples_dir, sample_path)
        # Import sample CSV
        sample_df = pandas.read_csv(tab_path, sep='\t')
        # Get variant caller names
        caller_names = get_caller_names(sample_df)
        # List of columns to keep
        col_list = [
                'CHROMOSOME', 'POSITION', 'TOTAL_CALLERS',
                'SNPEFF_ANNOTATED_GENE', 'SNPEFF_ANNOTATED_DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE', 'EXAC_FREQ_ESTIMATE',
                # '1KG_FREQ_ESTIMATE', 
                'TARGET_ALLELE_COUNT',
        ] + caller_names
        # New sample DataFrame with relevant columns only
        altered_df = sample_df[col_list].rename(columns={
                'SNPEFF_ANNOTATED_GENE': 'GENE',
                'SNPEFF_ANNOTATED_DNA_CHANGE': 'DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE': 'PROTEIN_CHANGE'
        })
        # Add SAMPLE_ID column (for when this DataFrame will be merged with others)
        altered_df['SAMPLE_ID'] = sample_df.shape[0] * [sample_id]
        altered_df['SAMPLE_PATH'] = sample_df.shape[0] * [sample_path]
        # Add sample DataFrame to list
        df_list.append(altered_df.reset_index(drop=True))
    df = pandas.concat(df_list).reset_index(drop=True)
    return df

def combine(db_name, bed_file, samples_dir):
    # File system
    main_dir = os.path.abspath(os.path.join(sys.path[0], '..'))
    output_dir = os.path.join(main_dir, 'output')
    gt_file = os.path.join(output_dir, 'ground_truth.csv')
    # Get ground truth
    if os.path.isfile(gt_file):
        print('Reading ground truth file...')
        gt = pandas.read_csv(gt_file, sep='\t')
        gt['SAMPLE_ID'] = gt['SAMPLE_ID'].astype(str)
    else:
        # Parse ground truth DataFrame from mongo database
        gt = get_simplified_gt(db_name, output_dir)
    # Parse .bed file
    bed = parse_bed(bed_file)
    print('Finding matching samples...')
    # Find sample ID matches in the ground truth
    sample_paths = find_samples(gt, samples_dir)
    # Combine ground truth DataFrames
    gt = gt.iloc[[
            i for i, path in enumerate(gt['SAMPLE_PATH']) 
            if path in sample_paths
    ]].reset_index(drop=True)
    print(gt.shape[0])
    print('Combining samples...')
    # Combine all sample DataFrames into one big DataFrame
    df = combine_samples(samples_dir, sample_paths).reset_index(drop=True)
    print('Adding combined callers...')
    # Add x or more
    # analysis.add_x_or_more(df)
    # List of variant caller names
    callers = get_caller_names(df)
    # Clean up ground truth
    # del gt['SAMPLE_PATH']
    del gt['SAMPLE_NAME']
    # Find variants in gt not found by any caller
    print('Finding false negatives...')
    # false_negs = gt[~gt.isin(df)].dropna()
    merge_list = [
            'SAMPLE_PATH', 'CHROMOSOME', 'POSITION', 'GENE', 'DNA_CHANGE',
            'PROTEIN_CHANGE', 'SAMPLE_ID'
    ]
    false_negs = gt.merge(df, on=merge_list, how='left', indicator=True)
    false_negs = false_negs.query('_merge == "left_only"').dropna(axis=1)
    for caller in callers:
        false_negs[caller] = ['./.'] * false_negs.shape[0]
    false_negs['TOTAL_CALLERS'] = [0] * false_negs.shape[0]
    del false_negs['_merge']
    false_negs.rename(columns={
            'EXAC_FREQ_ESTIMATE_x': 'EXAC_FREQ_ESTIMATE',
            'TARGET_ALLELE_COUNT_x': 'TARGET_ALLELE_COUNT'
    }, inplace=True)
    df = df.append(false_negs, ignore_index=True)
    true_pos = df.merge(gt, on=merge_list, how='left', indicator=True)
    true_pos = true_pos.query('_merge == "both"').dropna(axis=1)
    del true_pos['_merge']

    print('Splitting panel file...')
    panel = split_panel(bed)
    covered = df.merge(panel, on=['GENE','CHROMOSOME','POSITION'], how='left',
            indicator=True).query('_merge == "both"')

    # List of variants covered by the small panel
    covered_list = df.shape[0] * [False]
    # List of variants reported in the ground truth
    reportables_list = df.shape[0] * [False]

    # Set up progress indicator
#    sys.stdout.write('Generating combined DataFrame:  0%')
#    sys.stdout.flush()
#    print('Generating combined data frame...')

    # Cycle through every variant
#    for s in range(0, df.shape[0]):
        # Progress indicator
#        progress = int(100 * (s / df.shape[0]))
#        if progress < 10:
#            sys.stdout.write('\b\b' + str(progress) + '%')
#        else:
#            sys.stdout.write('\b\b\b' + str(progress) + '%')
#        sys.stdout.flush()

        # Find variants covered by the small panel
        # If the gene is covered by the small panel
#        if df['GENE'][s] in set(bed['GENE']):
#            covered_genes = [
#                    x for x, gene in enumerate(bed['GENE'])
#                    if gene == df['GENE'][s]
#            ]
#            # List of possible locations in small panel coverage file
#            for n in covered_genes:
#                # If variant position is in a location covered by the small panel
#                if (df['CHROMOSOME'][s] == bed['CHROMOSOME'][n]
#                        and int(df['POSITION'][s]) >= bed['START'][n]
#                        and int(df['POSITION'][s]) <= bed['END'][n]):
#                # if (df['CHROMOSOME'][s] == bed['CHROMOSOME'][n]
#                #         and int(df['POSITION'][s]) == bed['POSITION'][n]):
#                    covered_list[s] = True

        # Find variants reported in the ground truth
        # This if statement weeds out unnecessary checks
#        if (df['GENE'][s] in gt['GENE'].tolist()
#                and df['DNA_CHANGE'][s] in gt['DNA_CHANGE'].tolist()
#                and df['PROTEIN_CHANGE'][s] in gt['PROTEIN_CHANGE'].tolist()):
#            # Cycle through variants in ground truth
#            for g in range(0, gt.shape[0]):
#                # If the variant info matches between the DataFrame and ground truth
#                if (df['SAMPLE_ID'][s] == gt['SAMPLE_ID'][g]
#                        and int(df['POSITION'][s]) == int(gt['POSITION'][g])
#                        and df['GENE'][s] == gt['GENE'][g]
#                        and df['DNA_CHANGE'][s] == gt['DNA_CHANGE'][g]
#                        and df['PROTEIN_CHANGE'][s] == gt['PROTEIN_CHANGE'][g]):
#                    reportables_list[s] = True
        
        # Find true positives, etc.
#        for caller in callers:
#            call = ''
#            if reportables_list[s]:
#                if df[caller][s] == './.':
#                    call = 'FN' # false negative
#                else:
#                    call = 'TP' # true positive
#            elif covered_list[s]:
#                if df[caller][s] == './.':
#                    call = 'TN' # true negative
#                else:
#                    call = 'FP' # false positive
#            else:
#                if df[caller][s] == './.':
#                    call = 'UN' # unclassified negative
#                else:
#                    call = 'UP' # unclassified positive
#            df.set_value(s, caller, call)

#    sys.stdout.write('\b\b\b100%')
#    sys.stdout.flush

    # Add covered and reportable lists to DataFrame
#    df['COVERED'] = covered_list
#    df['REPORTABLE'] = reportables_list
    print('Generating reportables...')
    df['REPORTABLE'] = [(i in true_pos.index) for i in range(0, df.shape[0])] 
    df['COVERED'] = [(i in covered.index) for i in range(0, df.shape[0])]
    print('Classifying calls...')
    for caller in callers:
        print(caller)
        print('np')
        np = ['N' if call == './.' else 'P' for call in df[caller]]
        print('tf')
        tf = [False if call == './.' else True for call in df[caller]]
        print('calls')
        calls = [str(tf[i] == df['REPORTABLE'][i])[0] + np[i] for i in range(0, df.shape[0])]
        print('adding to df')
        df[caller] = calls
    # Create parsed tab file
    df.to_csv(
            os.path.join(output_dir, 'parsed.tab'), sep='\t',
            encoding='utf-8', index=False
    )
    print('\nOutput to file ' + os.path.join(output_dir, 'parsed.tab'))
    return df

if __name__ == "__main__":
    import argparse

    # Parser / command line instructions
    # Description of script
    arg_parser = argparse.ArgumentParser(
            description='Compares sample variants to\
            ground truth file and produces false negatives, false positives,\
            and true positives.'
    )
    # First argument: mongo database name
    arg_parser.add_argument('db_name', help='mongo db name', action='store')
    # Second argument: small panel gene list, bed file
    arg_parser.add_argument('panel_file', help='.bed', action='store')
    # Third argument: directory containing sample directories
    arg_parser.add_argument(
            'directory', help='directory containing all sample directories',
            action='store'
    )

    # Parse arguments
    args = arg_parser.parse_args()
    db_name = args.db_name
    bed_file = args.panel_file
    samples_dir = args.directory

    combine(db_name, bed_file, samples_dir)
