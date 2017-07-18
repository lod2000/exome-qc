from bson import ObjectId

import pymongo
import pandas

# Get mongo database
# Requires mongorestore to have been run
client = pymongo.MongoClient()
db = client['2017-06-30_NgsReviewer_master']

# Mongo collections
final = db.mongo_final.final_results
torrent = db.mongo_torrent.torrent_results

reportable_cursor = final.find({'reportable':True})
combined = []

def replace_key(dictionary, new_key, old_key):
    dictionary[new_key] = dictionary[old_key]
    del dictionary[old_key]

for final_doc in reportable_cursor:
    torrent_doc = torrent.find_one(final_doc['torrent_result_id'])

    # Generate combined dict of torrent and final information
    try:
        replace_key(final_doc, 'final_updated_at', 'updated_at')
        replace_key(torrent_doc, 'torrent_updated_at', 'updated_at')
        replace_key(torrent_doc, 'torrent_result_id', '_id')
        replace_key(final_doc, 'final_result_id', '_id')
        locus = (torrent_doc['_locus'])
        if 'snapshot' in torrent_doc:
            del torrent_doc['snapshot']

        combined.append({**torrent_doc,  **final_doc})
    except TypeError:
        print('No torrent result for ' + str(final_doc['torrent_result_id']))

df = pandas.DataFrame(combined)
df['CHROMOSOME'] = [locus.split(':')[0] for locus in df['_locus']]
df['POSITION'] = [locus.split(':')[1] for locus in df['_locus']]
df.to_csv('ground_truth.csv', sep='\t', encoding='utf8', index=False)
