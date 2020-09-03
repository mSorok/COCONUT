import pyclassyfire.client
from pymongo import MongoClient
import pandas as pd
import json

from time import sleep




def main():

    df = pd.read_csv("classification/COCONUT4MetFrag_ClassyFied.tsv", sep="\t")

    classification_dic = {}

    print("parsed data")

    print("start insertion in MongoDB")

    for index, row in df.iterrows():
        if row['InChIKey'] != "NaN" :
            inchikey = row['InChIKey']
        else:
            inchikey = ""

        if row['superclass.name'] != "NaN" :
            superclass = row['superclass.name']
        else:
            superclass = ""

        if row['class.name'] != "NaN" :
            cclass = row['class.name']
        else:
            cclass = ""

        if row['subclass.name']!= "NaN" :
            subclass = row['subclass.name']
        else:
            subclass = ""

        if row['direct_parent.name']  != "NaN" :
            direct_parent = row['direct_parent.name']
        else:
            direct_parent = ""

        classification_dic[inchikey] = [superclass, cclass, subclass, direct_parent]






    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    collection = db.uniqueNaturalProduct.aggregate(
            [{'$project': {'_id': 0, 'coconut_id': 1, 'inchikey': 1, "inchi": 1}}])


    allnp = pd.DataFrame(list(collection))

    classification_added_count = 0

    for index, row in allnp.iterrows():
        coconut_id = row['coconut_id']
        inchikey = row["inchikey"]
        inchi = row["inchi"]

        if inchikey in classification_dic.keys():
            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {
                "chemicalSuperClass": classification_dic[inchikey][0],
                "chemicalClass": classification_dic[inchikey][1],
                "chemicalSubClass": classification_dic[inchikey][2],
                "directParentClassification": classification_dic[inchikey][3]}})
            classification_added_count = classification_added_count + 1
        else:

            print("Nothing found in Classifire for " + inchikey)

            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {
                "chemicalSuperClass": "", "chemicalClass": "",
                "chemicalSubClass": "", "directParentClassification": ""}})

                # need to search with classyfire
                # try:
                #     # r = pyclassyfire.client.get_entity(inchikey, "json")
                #
                #     qid = pyclassyfire.client.structure_query(inchi)
                #     r = pyclassyfire.client.get_results(qid, 'json')
                #
                #     if len(r.keys()) > 0:
                #         cclass = ""
                #         direct_parent = ''
                #         inchikey = ""
                #         superclass = ""
                #         subclass = ""
                #
                #         if "class" in r:
                #             cclass = r["class"]["name"]
                #
                #         if "direct_parent" in r:
                #             direct_parent = cclass = r["direct_parent"]["name"]
                #
                #         if "inchikey" in r:
                #             inchikey = r["inchikey"].split("=")[1]  # "InChIKey=AAAACIDXVXFWLZ-UHFFFAOYSA-N"
                #
                #         if "superclass" in r:
                #             superclass = r["superclass"]["name"]
                #
                #         if "subclass" in r:
                #             subclass = r["subclass"]["name"]
                #
                #         db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {
                #             "chemicalSuperClass": superclass, "chemicalClass": cclass,
                #             "chemicalSubClass": subclass, "directParentClassification": direct_parent}})
                #         classification_added_count = classification_added_count + 1
                #
                #     else:
                #         db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {
                #             "chemicalSuperClass": "", "chemicalClass": "",
                #             "chemicalSubClass": "", "directParentClassification": ""}})
                # except Exception as e:
                #     print("Nothing found in Classifire for " + inchikey)
                #
                #     db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {
                #         "chemicalSuperClass": "", "chemicalClass": "",
                #         "chemicalSubClass": "", "directParentClassification": ""}})

    print("done")
    print("total classification added: "+str(classification_added_count))




if __name__ == '__main__':
    main()



