# Compiler for the COlleCtion of Open NatUral producTs (COCONUT)

This compiles is designed to read molecules from various file types (SMILES, SDF, MOL, csv, tsv), check them for errors and connectivity, compute a large number of molecular parameters and properties, and store everything in a Mongo database.
Unless you want to modify the code, we recomment using the compiled JAR that can be downloaded here: https://zenodo.org/record/3695455 

#### System pre-requisites:

- MongoDB installed and accessible by the default (27017) port on localhost
- Java minimum 8 version installed


#### Load COCONUT
You can download the latest version COCONUT from ZENODO (https://zenodo.org/record/3688734). In case you want to explore the whole database in MongoDB, you can load the downloaded dataset dump as following:
````bash
unzip COCONUT2020-07.zip
cd COCONUT2020-07/COCONUT2020-07/
mongorestore --db=COCONUT --noIndexRestore .
````

Note that for system compatibilities reason, it is better to restore the database without the indexes (hence the "noIndexRestore" option).
However, seen the size of the dataset, we suggest to add the indexes as following:

```
mongo
use COCONUT2020-07


 db.sourceNaturalProduct.createIndex( {source:1})

 db.sourceNaturalProduct.createIndex( {simpleInchi:"hashed"})

 db.sourceNaturalProduct.createIndex( {simpleInchiKey:1})
 db.sourceNaturalProduct.createIndex( {originalInchiKey:1})
 db.sourceNaturalProduct.createIndex( {originalSmiles:"hashed"})
 db.sourceNaturalProduct.createIndex( {absoluteSmiles:"hashed"})
 db.sourceNaturalProduct.createIndex( {idInSource:1})


db.uniqueNaturalProduct.createIndex( {inchi:"hashed"})
db.uniqueNaturalProduct.createIndex( {inchikey:1})
db.uniqueNaturalProduct.createIndex( {clean_smiles: "hashed"})
db.uniqueNaturalProduct.createIndex( {molecular_formula:1})
db.uniqueNaturalProduct.createIndex( {coconut_id:1})
db.uniqueNaturalProduct.createIndex( {fragmentsWithSugar:"hashed"})
db.uniqueNaturalProduct.createIndex( {fragments:"hashed"})
db.uniqueNaturalProduct.createIndex( {annotationLevel:1})
db.uniqueNaturalProduct.createIndex( {synonyms:"text", name:"text"})
db.uniqueNaturalProduct.createIndex( {npl_score:1})

db.uniqueNaturalProduct.createIndex( { pubchemBits : "hashed" } )

db.uniqueNaturalProduct.createIndex( {unique_smiles: "hashed"})

db.uniqueNaturalProduct.createIndex( { "pfCounts.bits" :1} )
db.uniqueNaturalProduct.createIndex( { "pfCounts.count" : 1 })


db.fragment.createIndex({signature:1})
db.fragment.createIndex({signature:1, withsugar:-1})



```

#### Required folder structure

```
COCONUT
├── coconut-0.0.1-SNAPSHOT.jar # the compiled jar. It can be downloaded from ZENODO: https://zenodo.org/record/3695455
├── coconut_ids_june2020.csv
├── data # here go the files with NP molecular structures
├── UpdateCOCONUT
│   ├── COCONUTupdater.py
│   ├── VerifyNames.py
├── fragments
│   ├── fragment_without_sugar.txt
│   ├── fragment_with_sugar.txt
├── sm # file(s) with structures of synthetic molecules 
```

#### Execution options

##### Run COCONUT compilation from scratch

```bash
java -Xmx16288m -jar coconut-0.0.1-SNAPSHOT.jar data sm/sm.tsv fragments/fragment_without_sugar.txt fragments/fragment_with_sugar.txt importCOCONUTids coconut_ids_june2020.csv > logs.txt &
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


##### run only import of CNPid from file
java -jar coconut-0.0.1-SNAPSHOT.jar onlyImportCoconutIds coconut_ids_june2020.csv &

## Molecule names curation
##### Import names from ChEBI, PubChem and CMAUP

```bash
python3 UpdateCOCONUT/COCONUTupdater.py
```

##### Import IUPAC names and clean

```bash
python3 UpdateCOCONUT/VerifyNames.py
```

##### Run annotation level recalculation  
```bash
java -Xmx16288m -jar coconut-0.0.1-SNAPSHOT.jar evaluateAnnotation &
```

##### Names to low case
```bash
java -Xmx16288m -jar coconut-0.0.1-SNAPSHOT.jar namesToLowerCase &
```


