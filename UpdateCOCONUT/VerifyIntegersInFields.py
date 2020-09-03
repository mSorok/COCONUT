#!/usr/bin/python -tt

from pymongo import MongoClient

from pymongo.errors import CursorNotFound
import re

from time import sleep





def main():


    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    # for each np in the database, check correct cross-link to other ressources

    aggregate = [{'$project': {'_id': 0, 'coconut_id': 1}}]

    np_ids = db.uniqueNaturalProduct.aggregate(aggregate)
    # print(np_ids['coconut_id'])

    print("retrieved all coconut ids")
    # print(np_ids)

    client = MongoClient("localhost:27017")
    db = client['COCONUT2020-07']

    try:
        for r in np_ids:
            coconut_id = r['coconut_id']

            np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})
            np = np[0]

            #check taxid field

            if len(np['taxid'])>0:
                for tid in np['taxid']:
                    if "-" in tid:
                        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$pull": {"taxid": tid}})
                        ntid = tid.split("-")[0]
                        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id},
                                                           {"$addToSet": {"taxid": ntid}})

    except CursorNotFound:
        pass






    client.close()
    print("done")


if __name__ == '__main__':
    main()


