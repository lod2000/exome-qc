import re
import os
import time
import sys

import pandas
import pymongo
import numpy
from bson import ObjectId

import utils

def replace_key(dictionary, new_key, old_key):
    dictionary[new_key] = dictionary[old_key]
    del dictionary[old_key]

# Get DataFrame of variants deemed "reportable" in final mongo database
def get_reportables(db_name, hostname):
    # Get mongo database
    # Requires mongorestore to have been run
    client = pymongo.MongoClient(hostname)
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
    for i, final_doc in enumerate(reportable_cursor):
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
                    torrent_doc['ppmp_annotation_path'] = server_doc[
                            'ppmp_annotation_path'
                    ]
            # Combine torrent and final documents
            combined.append({**torrent_doc,  **final_doc})
        # In case there isn't a matching torrent result
        except TypeError:
            print('No torrent result for ' 
                    + str(final_doc['torrent_result_id'])
            )
    # Generate combined DataFrame
    df = pandas.DataFrame(combined)
    return df

# Simplify reportable ground truth data to contain only relevant information
# Also outputs to ground_truth.csv
def get_simplified_gt(db_name, hostname, gt_file):
    df = get_reportables(db_name, hostname)
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
            ppmp_annotation_paths.append(df['ppmp_annotation_path'][i]
                    .split('Results/')[-1]
            )
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
    simple_df.to_csv(gt_file, sep='\t', encoding='utf8', index=False)
    return simple_df

# Returns parsed list of positions covered by small panel test, taken from 
# a .bed file
def parse_bed(bed):
    # Read bed file
    panel = pandas.read_csv(
            bed, names=['CHROMOSOME', 'START', 'END', 'GENE'], sep='\t'
    )
    # Split gene from location string
    genes = [
            re.search('_[A-Z0-9]*_', gene_string).group(0).split('_')[1]
            for gene_string in panel['GENE']
    ]
    panel['GENE'] = genes
    positions = []
    genes = []
    chromosomes = []
    # Split position ranges into one row per position
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

# Returns list of sample paths from ground truth that are in data/samples
def find_samples(gt, samples_dir):
    samples = next(os.walk(samples_dir))[1]
    # Remove duplicate sample paths
    potentials = [
            path for i, path in enumerate(gt['SAMPLE_PATH'])
            if path not in list(gt['SAMPLE_PATH'][:i])
    ]
    finals = []
    for potential in potentials:
        if os.path.isfile(os.path.join(samples_dir, str(potential))):
            finals.append(os.path.join(*potential.split('/')))            
    return finals

# Returns a DataFrames of samples for which matches exist in the ground truth
def combine_samples(samples_dir, sample_paths):
    # List of sample DataFrames
    df_list = []
    for sample_path in sample_paths:
        sample_id = sample_path.split(os.sep)[0]
        tab_path = os.path.join(samples_dir, sample_path)
        # Import sample CSV
        sample_df = pandas.read_csv(tab_path, sep='\t')
        # Get variant caller names
        caller_names = utils.get_og_callers(sample_df)
        # List of columns to keep
        col_list = [
                'CHROMOSOME', 'POSITION', 'TOTAL_CALLERS',
                'SNPEFF_ANNOTATED_GENE', 'SNPEFF_ANNOTATED_DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE', 'EXAC_FREQ_ESTIMATE',
                'TARGET_ALLELE_COUNT', 'TOTAL_READS'
        ] + caller_names
        # New sample DataFrame with relevant columns only
        altered_df = sample_df[col_list].rename(columns={
                'SNPEFF_ANNOTATED_GENE': 'GENE',
                'SNPEFF_ANNOTATED_DNA_CHANGE': 'DNA_CHANGE',
                'SNPEFF_ANNOTATED_PROTEIN_CHANGE': 'PROTEIN_CHANGE'
        })
        # Add SAMPLE_ID column (for when this DataFrame will be merged with 
        # others)
        altered_df['SAMPLE_ID'] = sample_df.shape[0] * [sample_id]
        altered_df['SAMPLE_PATH'] = sample_df.shape[0] * [sample_path]
        # Add sample DataFrame to list
        df_list.append(altered_df.reset_index(drop=True))
    df = pandas.concat(df_list).reset_index(drop=True)
    return df

# Find true and false positives and negatives
def classify(df, callers):
    for caller in callers:
        print(caller)
        # Change call to 'N' if negative or 'P' if positive hit
        np = ['N' if call == './.' else 'P' for call in df[caller]]
        """
        TP = positive and reportable
        FP = positive, not reportable, covered
        TN = negative, not reportable, covered,
        FN = negative, reportable
        UP or UN = unclassified positive or negative

        Example: 'N', reportable, covered
        (np[i] == 'P') is False and df['REPORTABLE'][i] is True
        (False == True) is False, so str(False)[0] is 'F'
        'F' + 'N' is 'FN'
        """
        reportables = df['REPORTABLE']
        covereds = df['COVERED']
        size = df.shape[0]
        df[caller] = [
                str((np[i] == 'P') == reportables[i])[0] + np[i] 
                if (covereds[i] or reportables[i]) else 'U' + np[i]
                for i in range(0, size)
        ]

def combine(db_name, hostname, bed_file, samples_dir):
    # File system
    main_dir = os.path.abspath(os.path.join(sys.path[0], '..'))
    output_dir = os.path.join(main_dir, 'output')
    gt_file = os.path.join(
            main_dir, 'data', 'ground_truth', 'ground_truth.csv'
    )

    # Determine whether to use csv or mongo db
    use_csv = False
    if os.path.isfile(gt_file):
        use_csv = utils.query_yes_no(
                'Found ground truth csv. Use it instead of mongo?'
        )

    # Get ground truth
    if use_csv:
        print('Reading ground truth csv...')
        gt = pandas.read_csv(gt_file, sep='\t')
        gt['SAMPLE_ID'] = gt['SAMPLE_ID'].astype(str)
    else:
        # Parse ground truth DataFrame from mongo database
        gt = get_simplified_gt(db_name, hostname, gt_file)

    # Find sample ID matches in the ground truth
    print('Finding matching samples...')
    sample_paths = find_samples(gt, samples_dir)
    # Combine ground truth DataFrames
    gt_sample_paths = gt['SAMPLE_PATH']
    gt = gt.iloc[[
            i for i, path in enumerate(gt_sample_paths) 
            if path in sample_paths
    ]].reset_index(drop=True)
    del gt['SAMPLE_NAME']
    print('Combining samples...')
    # Combine all sample DataFrames into one big DataFrame
    df = combine_samples(samples_dir, sample_paths).reset_index(drop=True)
    print('Adding combined callers...')
    # List of variant caller names
    callers = utils.get_og_callers(df)

    # Find variants in gt not found by any caller
    print('Finding false negatives...')
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

    # Find variants that are also present in ground truth
    print('Generating reportables...')
    true_pos = df.merge(gt, on=merge_list, how='left', indicator=True)
    true_pos = true_pos.query('_merge == "both"').dropna(axis=1)
    del true_pos['_merge']
    df['REPORTABLE'] = [(i in true_pos.index) for i in range(0, df.shape[0])] 

    # Find positions that are covered by the small panel
    panel = parse_bed(bed_file)
    covered = df.merge(
            panel, on=['GENE','CHROMOSOME','POSITION'], how='left',
            indicator=True).query('_merge == "both"'
    )
    df['COVERED'] = [(i in covered.index) for i in range(0, df.shape[0])]

    # Classify calls (TP, FN, etc.)
    print('Classifying calls...')
    classify(df, utils.get_og_callers(df))

    # Create parsed tab file
    df.to_csv(
            os.path.join(output_dir, 'combined.tab'), sep='\t',
            encoding='utf-8', index=False
    )
    print('Output to file ' + os.path.join(output_dir, 'combined.tab'))
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
    # Generate combined file
    combine(db_name, bed_file, samples_dir)
