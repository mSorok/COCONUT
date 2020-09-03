import deepsmiles

from pymongo import MongoClient
import pandas as pd

def main():


    converter = deepsmiles.Converter(rings=True, branches=True)

    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    collection = db.uniqueNaturalProduct.aggregate([{'$project': {'_id': 0, 'coconut_id': 1, 'unique_smiles': 1}}])

    allnp = pd.DataFrame(list(collection))

    for index, row in allnp.iterrows():
        coconut_id = row['coconut_id']
        smiles = row["unique_smiles"]

        deep_smiles = converter.encode(smiles)
        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"deep_smiles": deep_smiles}})

        # for each row - get unique_smiles -> convert to deep smiles and save at deep_smiles

    # db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": data[2]}})



    print("done")


if __name__ == '__main__':
    main()



