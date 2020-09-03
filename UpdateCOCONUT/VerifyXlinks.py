#!/usr/bin/python -tt

from pymongo import MongoClient

from pymongo.errors import CursorNotFound
import re

from time import sleep





def main():

    source_pretty_names = { "bitterdb":"BitterDB",
                            "carotenoids":"Carotenoids Database",
                            "chebi_np": "ChEBI",
                            "chebinp": "ChEBI",
                            "chembl_np": "ChEMBL",
                            "chemblnp": "ChEMBL",
                            "cmaup": "CMAUP",
                            "pubchem_tested_np": "PubChem",
                            "pubchemnp": "PubChem",
                            "drugbanknp": "DrugBank",
                            "chemspidernp": "ChemSpider",
                            "np_atlas_2019_12": "NPAtlas",
                            "npatlas": "NPAtlas",
                            "exposome-explorer": "Exposome Explorer",
                            "fooddb": "FooDB",
                            "knapsack": "KnapSack",
                            "npass": "NPASS",
                            "nubbe": "NuBBE",
                            "phenolexplorer": "Phenol Explorer",
                            "sancdb": "SANCDB",
                            "supernatural2": "SuperNatural 2",
                            "tcmdb_taiwan": "TCMDB@Taiwan",
                            "tppt": "TPPT",
                            "vietherb": "VietHerb",
                            "streptomedb": "StreptomeDB"
                            }

    pretty_list = [ "BitterDB", "Carotenoids Database", "ChEBI", "ChEMBL", "CMAUP",  "PubChem", "DrugBank", "ChemSpider", "NPAtlas" ,  "Exposome Explorer", "FooDB", "KnapSack", "NPASS", "NuBBE","Phenol Explorer", "SANCDB", "SuperNatural 2", "TCMDB@Taiwan", "TPPT", "VietHerb","StreptomeDB"]


    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    # for each np in the database, check correct cross-link to other ressources

    aggregate =[ { '$project':{  '_id':0, 'coconut_id':1 }}]

    np_ids = db.uniqueNaturalProduct.aggregate(aggregate)
    #print(np_ids['coconut_id'])

    print("retrieved all coconut ids")
    #print(np_ids)

    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']


    print("verifying Xlinks for each molecule")
    try:
        for r in np_ids:

            coconut_id = r['coconut_id']

            np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})
            np = np[0]

            if "xrefs" in np.keys() and "clean_xrefs" not in np.keys():

                for xref in np["xrefs"]:

                    #db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$pull": {"xrefs": xref}})

                    if len(xref) == 3:
                        print(coconut_id)

                        if xref[0] in source_pretty_names.keys():
                            print(xref[0])


                            clean_xref = {"source": source_pretty_names[xref[0]], "id_in_source": xref[1], "link_to_source": xref[2]}

                            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$push": {"clean_xrefs": clean_xref}})

                            print(clean_xref)


            #sleep(0.05)

    except CursorNotFound:
        pass






    client.close()
    print("done")


if __name__ == '__main__':
    main()


