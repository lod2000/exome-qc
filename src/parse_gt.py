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

def get_simplified_gt(db_name, gt_dir):
    df = get_reportables(db_name, gt_dir)
    print('Generating simplified data frame...')
    reported = []
    sample_ids = []
    chromosomes = []
    positions = []
    genes = []
    dna_changes = []
    protein_changes = []
    headers = list(df)
    for i in range(0, df.shape[0]):
        sample_ids.append('')
        chromosomes.append(df['_locus'][i].split(':')[0])
        positions.append(df['_locus'][i].split(':')[1])
        # Try to get it all from one place
        variant = str(df['reported_variant'][i])
        if re.search('^[A-Z0-9]*\sc.[^\s]*:', variant):
            reported.append(True)
            genes.append(variant.split()[0])
            dna_changes.append(variant.split()[1].split(':')[0])
            try:
                protein_changes.append(re.search('p\..*$', variant).group(0))
            except AttributeError:
                protein_changes.append('')
        else:
            reported.append(False)
            # Gene
            genes.append('')
            gene_columns = [h for h in headers if re.search('gene', h)]
            for column in gene_columns:
                if re.search('^[A-Z]', str(df[column][i])):
                    genes[i] = df[column][i]
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
            # if re.search('^[^\s]*:c.[0-9]', str(df['hgvs_coding_change'][i])):
            #     dna_changes.append(df['hgvs_coding_change'][i].split(':')[1])
            # else:
            #     dna_changes.append('')
            protein_changes.append('')
    simple_df = pandas.DataFrame({
            'SAMPLE_ID': sample_ids, 'CHROMOSOME': chromosomes,
            'POSITION': positions, 'GENE': genes, 'DNA_CHANGE': dna_changes,
            'PROTEIN_CHANGE': protein_changes, 'REPORTED': reported
    })
    simple_df.to_csv(
            os.path.join(gt_dir, 'simple_ground_truth.csv'),
            sep='\t', encoding='utf8', index=False
    )
    return simple_df

if __name__ == "__main__":
    get_simplified_gt(
            '2017-06-30_NgsReviewer_master',
            os.path.join(sys.path[0], '..', 'data', 'ground_truth')
    )
