#!/usr/bin/python -tt

from pymongo import MongoClient
import re

from pymongo.errors import CursorNotFound


def main():


    good_sources = ["chebi_np", "chembl_np", "cmaup", "np_atlas_2019_12", "npatlas", "piellabdata", "knapsack"]

    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    aggregate = [{'$project': {'_id': 0, 'coconut_id': 1}}]
    np_ids = db.uniqueNaturalProduct.aggregate(aggregate)

    try:
        for r in np_ids:

            coconut_id = r['coconut_id']
            #print(coconut_id)
            np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})
            np = np[0]

            has_name = False
            has_literature = False
            has_organism = False
            has_good_source = False

            if 'name' in np.keys() and 'iupac_name' in np.keys() and np['name'] != np['iupac_name']:
                has_name = True

            if "textTaxa" in np.keys() and "notax" not in np["textTaxa"] and len(np["textTaxa"])>0:
                has_organism = True

            if "citationDOI" in np.keys() and len(np["citationDOI"])>0:
                has_literature = True

            if "found_in_databases" in np.keys() :
                for source in np["found_in_databases"]:
                    if source in good_sources:
                        has_good_source = True

            annotation_level = sum([has_name,has_literature, has_organism, has_good_source ])+1
            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"annotationLevel": annotation_level}})


    except CursorNotFound:
        pass

    print("all done")


if __name__ == '__main__':
    main()

