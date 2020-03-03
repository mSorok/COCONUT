# Compiler for the COlleCtion of Open NatUral producTs (COCONUT)

This compiles is designed to read molecules from various file types (SMILES, SDF, MOL, csv, tsv), check them for errors and connectivity, compute a large number of molecular parameters and properties, and store everything in a Mongo database.
Unless you want to modify the code, we recomment using the compiled JAR that can be downloaded here: https://zenodo.org/record/3695455 

#### System pre-requisites:

- MongoDB installed and accessible by the default (27017) port on localhost
- Java minimum 8 version installed


#### Load COCONUT
You can download the latest version COCONUT from ZENODO (https://zenodo.org/record/3688734). In case you want to explore the whole database in MongoDB, you can load the downloaded dataset dump as following:
````bash
unzip COCONUTv2.zip
cd COCONUTv1/COCONUT/
mongorestore --db=COCONUT --noIndexRestore .
````

Note that for system compatibilities reason, it is better to restore the database without the indexes (hence the "noIndexRestore" option).
However, seen the size of the dataset, we suggest to add the indexes as following:

````bash
mongo
use COCONUT
db.sourceNaturalProduct.createIndex( {source:1})
db.sourceNaturalProduct.createIndex( {simpleInchi:1})
db.sourceNaturalProduct.createIndex( {simpleInchiKey:1})
db.uniqueNaturalProduct.createIndex( {inchi:1})
db.uniqueNaturalProduct.createIndex( {inchikey:1})
db.uniqueNaturalProduct.createIndex( {smiles:1})
db.uniqueNaturalProduct.createIndex( {clean_smiles:1})
db.uniqueNaturalProduct.createIndex( {molecular_formula:1})
db.uniqueNaturalProduct.createIndex( {name:1})
db.uniqueNaturalProduct.createIndex( {coconut_id:1})
db.fragment.createIndex({signature:1})
db.fragment.createIndex({signature:1, withsugar:-1})
````

#### Required folder structure

````bash
COCONUT
├── coconut-0.0.1-SNAPSHOT.jar # the compiled jar. It can be downloaded from ZENODO: https://zenodo.org/record/3695455
├── data # here go the files with NP molecular structures
├── fragments
│   ├── fragment_without_sugar.txt
│   ├── fragment_with_sugar.txt
├── sm # file(s) with structures of synthetic molecules 
````

#### Execution options

##### Run COCONUT compilation from scratch

```bash
java -jar coconut-0.0.1-SNAPSHOT.jar ~/Projects/NP/COCONUT/data ~/Projects/NP/COCONUT/sm/sm.tsv ~/Projects/NP/COCONUT/fragments/fragment_without_sugar.txt ~/Projects/NP/COCONUT/fragments/fragment_with_sugar.txt &
```

##### Re-run COCONUT to recompute missing molecular features 

```bash
java -jar coconut-0.0.1-SNAPSHOT.jar recomputeMissing &
```

##### Run only similarity calculation
The provided script also computes similatiry between molecules in COCONUT. This step requires a large amount of memory (this will be optimized in the future) and might need a separate calculation.
````bash
java -Xmx12288m -jar coconut-0.0.1-SNAPSHOT.jar runOnlySimilarity &
````

##### Run only addition of synthetic molecules
Synthetic molecules (SM) are required for a large number of comparisons with NPs (for example a re-calculation from scratch of NP-likeness score). You need to provide your own dataset of synthetic molecules (SM), we suggest using the ZINC15 dataset.
Note that not adding any SM will not affect COCONUT, unless there is a need of re-calculating the NP-likeness score from scratch.

```bash
java -jar coconut-0.0.1-SNAPSHOT.jar onlyAddSM ~/Projects/NP/COCONUT/sm/sm.tsv & 
```

##### Run only CNPid calculation
In case you need to produce your own CNPid (COCONUT NP identifiers) - not recommened.

```bash
java -jar coconut-0.0.1-SNAPSHOT.jar addCNPid  & 
```


