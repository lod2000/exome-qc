from bson import ObjectId

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

    for final_doc in reportable_cursor:
        torrent_doc = torrent.find_one(final_doc['torrent_result_id'])
        # Generate combined dict of torrent and final information
        try:
            replace_key(final_doc, 'final_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_updated_at', 'updated_at')
            replace_key(torrent_doc, 'torrent_result_id', '_id')
            replace_key(final_doc, 'final_result_id', '_id')
            locus = (torrent_doc['_locus'])
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
    # Simplify DataFrame
    df['CHROMOSOME'] = [locus.split(':')[0] for locus in df['_locus']]
    df['POSITION'] = [locus.split(':')[1] for locus in df['_locus']]
    df.drop('_locus', axis=1, inplace=True)

if __name__ == "__main__":
    get_reportables(
            '2017-06-30_NgsReviewer_master',
            os.path.join(sys.path[0], '..', 'data', 'ground_truth')
    )