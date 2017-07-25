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
    server = db.mongo_server.server_downloads
    # Generate cursor for reportable final documents with a locus
    reportable_cursor = final.find({
            'reportable': True, 
            'germline': False, 
            'ion_junk': False,
            '_locus': {'$ne': ':'}
    })
    combined = []
    print('Combining final and torrent documents...')
    for final_doc in reportable_cursor:
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
    # Output to CSV
    df.to_csv(
            os.path.join(gt_dir, 'ground_truth.csv'),
            sep='\t', encoding='utf8', index=False
    )
    return df

def get_simplified_gt(db_name, gt_dir):
    df = get_reportables(db_name, gt_dir)
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
            'PROTEIN_CHANGE': protein_changes, 'PATH': ppmp_annotation_paths 
    })
    # Strip out variants without DNA or protein information
    simple_df = simple_df.iloc[[
            i for i, dna in enumerate(simple_df['DNA_CHANGE']) 
            if not dna == '' and not simple_df['PROTEIN_CHANGE'][i] == ''
    ]]
    # Output CSV
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
