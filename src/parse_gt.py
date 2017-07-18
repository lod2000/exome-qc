from bson import ObjectId
import os
import sys
import re

import pymongo
import pandas

def replace_key(dictionary, new_key, old_key):
    dictionary[new_key] = dictionary[old_key]
    del dictionary[old_key]

def get_reportables(db_name, gt_dir):
    # Get mongo database
    # Requires mongorestore to have been run
    client = pymongo.MongoClient()
    db = client[db_name]

    # Mongo collections
    final = db.mongo_final.final_results
    torrent = db.mongo_torrent.torrent_results

    reportable_cursor = final.find({'reportable': True, 'germline': False, 'ion_junk': False})
    combined = []
    print('Combining final and torrent documents...')
    for final_doc in reportable_cursor:
        torrent_doc = torrent.find_one(final_doc['torrent_result_id'])
        # Generate combined dict of torrent and final information
        try:
            replace_key(final_doc, 'final_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_result_id', '_id')
            replace_key(final_doc, 'final_result_id', '_id')
            # locus = (torrent_doc['_locus'])
            # Remove unnecessary fields
            if 'snapshot' in torrent_doc:
                del torrent_doc['snapshot']
            del final_doc['ion_junk']
            del final_doc['reportable']
            del final_doc['germline']
            # Combine dicts
            combined.append({**torrent_doc,  **final_doc})
        except TypeError:
            print('No torrent result for ' + str(final_doc['torrent_result_id']))

    df = pandas.DataFrame(combined)
    # Output to CSV
    df.to_csv(
            os.path.join(gt_dir, 'ground_truth.csv'),
            sep='\t', encoding='utf8', index=False
    )
    return df

def get_gene_id(df, index):
    headers = list(df)
    columns = [h for h in headers if re.search('gene', h)]
    for column in columns:
        if re.search('^[A-Z]', str(df[column][index])):
            return df[column][index]
    return ''

def get_simplified_gt(db_name, gt_dir):
    df = get_reportables(db_name, gt_dir)
    print('Generating simplified data frame...')
    # Simplify DataFrame
    simple_df = pandas.DataFrame()
    print('a')
    # Method 1
    chromosomes = []
    positions = []
    genes = []
    dna_changes = []
    protein_changes = []
    headers = list(df)
    gene_columns = [h for h in headers if re.search('gene', h)]
    for i in range(0, df.shape[0]):
        chromosomes.append(df['_locus'][i].split(':')[0])
        positions.append(df['_locus'][i].split(':')[1])
        if re.search('^[A-Z0-9]*\sc.[^\s]*:', str(df['reported_variant'][i])):
            genes.append(df['reported_variant'][i].split()[0])
            dna_changes.append(df['reported_variant'][i].split()[1].split(':')[0])
            try:
                protein_changes.append(df['reported_variant'][i].split()[1].split(':')[1])
            except AttributeError:
                protein_changes.append('')
        else:
            genes.append('')
            for column in gene_columns:
                if re.search('^[A-Z]', str(df[column][i])):
                    genes[i] = df[column][i]
            dna_changes.append('')
            protein_changes.append('')
    simple_df['CHROMOSOME'] = chromosomes
    simple_df['POSITION'] = positions
    simple_df['GENE'] = genes
    simple_df['DNA_CHANGE'] = dna_changes
    simple_df['PROTEIN_CHANGE'] = protein_changes
    print('b')
    # Method 2
    simple_df['CHROMOSOME'] = [locus.split(':')[0] for locus in df['_locus']]
    simple_df['POSITION'] = [locus.split(':')[1] for locus in df['_locus']]
    simple_df['GENE'] = [get_gene_id(df, i) for i in range(0, df.shape[0])]
    print('c')
    simple_df.to_csv(
            os.path.join(gt_dir, 'simple_ground_truth.csv'),
            sep='\t', encoding='utf8', index=False
    )

if __name__ == "__main__":
    get_simplified_gt(
            '2017-06-30_NgsReviewer_master',
            os.path.join(sys.path[0], '..', 'data', 'ground_truth')
    )
