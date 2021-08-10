#!/usr/bin/python -tt

from pymongo import MongoClient
import re


def nameIsWeird(name):
    inchikeyPattern = r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"

    casPattern = r"^[0-9]+-[0-9]+-[0-9]+^"

    if name == "" or len(name) <= 3 or name.startswith("MCULE") or name.startswith("STL") or "n.a." in name.lower() or name.startswith("STK") or name.startswith("ZINC")\
                or "chembl" in name.lower() or name.startswith("CTK") or name.lower().startswith("prestwick") or name.lower().startswith("eim")\
                or name.startswith("AKOS") or name.lower().startswith("mls") or name.startswith("ACMC") or name.startswith("SCHEMBL") \
                or name.startswith("NSC") or "dataset" in name.lower() or name.startswith("BMS") or name.startswith("EINECS") or "named" in name.lower() \
                or name.startswith("NCI") or name.startswith("HMS") or name.lower().startswith("cbmicro") or name.lower().startswith("oprea")  \
                or name == "ms 681c" or name == "tpu-0037-c" or name.lower().startswith("cbdive") or name == "J3.542.783J" or name == "F0848-0042" \
                or name == "CK-636" or "megapho" in name.lower() or name.startswith("sb 253") or name.startswith("CCG-") or "nadp" in name.lower() \
                or name == "jbir-28" or name == "F0848-0042" or name.startswith("DA-") or "glas#" in name or name == "jbir-67" or name == "wss2222" \
                or name.startswith("NCGC") or name.startswith("MFCD") or name.startswith("SDCCG") or "refreshn" in name.lower() or "pharmaco" in name.lower() \
                or name.startswith("BRD") or name.startswith("BIM") or name.startswith("NINDS") or name.startswith("CFF") or "ampule" in name \
                or name.lower().startswith("bio1") or name.lower().startswith("kbio") or name.startswith("F3371") or name.startswith("PDSP") \
                or "testing" in name.lower() or "specifications" in name.lower() or "melting" in name.lower() or "point" in name.lower() \
                or name.startswith("UNII") or "EP2" in name or "thermosystem" in name or "component" in name or name.startswith("CCRIS") \
                or name.startswith("BDBM") or name.startswith("ST0") or name == "C 0750" or "pharmaceutical" in name.lower() or name.startswith("EU-") \
                or name.startswith("BIDD") or  "99%" in name or "100%" in name or "98%" in name or name.startswith("DB-") or name == "No Doz" \
                or "special grade" in name or name.startswith("PubChem") or name.startswith("BRN") or name.startswith("InChI") or name.startswith("HSDB") \
                or "vial of" in name or name.startswith("bmse") or "tested" in name or name.startswith("ARONI") or name.startswith("L00") \
                or "einecs" in name.lower() or name.startswith("IDI1") or name.startswith("EC ") or name.startswith("07E4") or name.startswith("BSP") \
                or name.startswith("SBI") or "chebi" in name.lower() or "reference" in name.lower() or name.startswith("AB0") or "spectrum" in name.lower() \
                or "compound" in name.lower() or name.startswith("AS-") or name == "O926" or "synthetic" in name.lower() or "specified" in name.lower() \
                or name.startswith("ACN") or "molmap" in name.lower() or "HPLC" in name or name.startswith("CS-") or "tracecert" in name.lower() \
                or name.startswith("AI3") or "lopac" in name.lower() or name == "1l5q" or name.startswith("SR-") or name.startswith("DSS") or "microg/mL" in name.lower() \
                or "bayer" in name.lower() or name == "1l7x" or name.lower().startswith("divk") or name.lower().startswith("dtxs") or "production" in name.lower() \
                or "BAN:JAN" in name or "WLN:" in name or "mg/ml" in name.lower() or "probes" in name.lower() or "pharmagrade" in name.lower() or name == "FCC" \
                or "reagent" in name.lower() or name == "H2815" or name.startswith("GTPL") or "standard" in name.lower() or name.startswith("SBB") \
                or name == "Nix Nap" or "powder" in name.lower() or name.startswith("Tox2") or "tri-aqu" in name.lower() or "JP17" in name or name.startswith("CU-") \
                or "bioxtra" in name.lower() or name.endswith("-") or "caswell" in name.lower() or name.startswith("SMR") or name.startswith("LS-") \
                or name.startswith("\.") or "salt" in name.lower():
        return True
    else:
        return False




def reevaluateAnnotationLevel(np, db):

    hasName = False
    hasOrganism = False
    hasLiterature = False
    hasTrustedSource = False

    if db == "chebi" or db == "cmaup":
        hasTrustedSource = True

    if np["name"] is not None and np["name"] != "":
        hasName = True

    if np["textTaxa"] is not None and len(np["textTaxa"])>0 and "notax" not in np["textTaxa"] :
        hasOrganism = True

    if np["citationDOI"] is not None and len(np["citationDOI"])>0 :
        hasLiterature = True

    annotationLevel = 1 + sum([hasName,hasTrustedSource,hasOrganism, hasLiterature  ])

    return annotationLevel



def main():

    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    ##################################################################
    # updating with data from Knapsack
    print("Curating from KnapSack")

    f_knapsack = open("Knapsack_data_for_coconut.tsv", "r")

    for line in f_knapsack:
        line = line.replace("\n", "")
        data = line.split("\t")

        np = db.uniqueNaturalProduct.find({"coconut_id": data[0]})

        np = np[0]

        # add present in knapsack
        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"found_in_databases": "knapsack"}})


        # species
        species = data[2]
        if species != "":
            #species
            # remove "notax"
            db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$pull": {"textTaxa": "notax"}})

            # add the speciesc
            db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"textTaxa": species}})

        #reference
        ref = data[3]
        if ref != "":
            #text reference
            db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"citationDOI": ref}})

        # update XREFS
        if "xrefs" in np.keys():
            knapsack_id_set = False
            for l in np["xrefs"]:
                if l[0] == data[1]:
                    knapsack_id_set = True
            if not knapsack_id_set:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {
                    "xrefs":  ["knapsack", data[1], "http://www.knapsackfamily.com/knapsack_core/information.php?sname=C_ID&word="]}})



    f_knapsack.close()


    ##################################################################
    # updating with data from ChEBI
    print("Curating from ChEBI")
    f_chebi = open("coconut_chebi_total.tsv", "r")

    for line in f_chebi:
        line = line.replace("\n", "")

        data = line.split("\t")

        np = db.uniqueNaturalProduct.find({"coconut_id" : data[0]})

        np = np[0]

        #check if present in chebi
        # add present in chebi
        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"found_in_databases": "chebi_np"}})


        # UPDATE NAME
        if np['nameTrustLevel'] < 2 :
            if nameIsWeird(np["name"]) :
                if data[2] != "":
                    #update name
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": data[2]}})
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]},{"$set": {"nameTrustLevel": 2}})

                elif (np['name'] is None or np['name']=="" ) and data[3] != "":
                    #put first synonym
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": data[3].split("$$")[0]}})
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"nameTrustLevel": 2}})
            else:
                if len(np["name"]) > len(data[2]): # replace by Chebi name only if it's shorter!
                    if data[2] != "":
                        # update name
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"synonyms": np["name"]}})
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": data[2]}})
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"nameTrustLevel": 2}})

                    elif (np['name'] is None or np['name'] == "") and data[3] != "":
                        # put first synonym
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": data[3].split("$$")[0]}})
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"nameTrustLevel": 2}})
                else:
                    if data[2] != "":
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"synonyms": data[2]}})
        else:
            if data[2] != "":
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"synonyms": data[2]}})



        if data[3] != "":
            synonyms = data[3].split("$$")
            for s in synonyms:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"synonyms": s}})





        if data[4] != "":
            #species
            species = data[4].split("$$")
            # remove "notax"
            db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$pull": {"textTaxa": "notax"}})

            for s in species:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"textTaxa": s}})

        if data[5] != "":
            #taxid
            species = data[5].split("$$")
            for s in species:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"taxid": s}})

        if data[6] != "":
            #pubmed id
            pubmed = data[6].split("$$")
            for p in pubmed:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"citationDOI": p}})


        #update XREFS
        if "xrefs" in np.keys():
            chebi_id_set = False
            for l in np["xrefs"]:
                if l[0]==data[1]:
                    chebi_id_set = True
            if not chebi_id_set:
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"xrefs": [  "chebi_np", data[1],"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:" ]}})


        #update annotation level
        #annotationLevel = reevaluateAnnotationLevel(np, "chebi")
        #db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"annotationLevel":annotationLevel}})


    f_chebi.close()

    ##################################################################
    # updating with data from CMAUP
    print("Curating from CMAUP")
    # retrieving coconut-cmaup CPD


    cocodata_cpd = {}
    cmaup_data = {}
    cmaup_plant_cpd = {}
    cmaup_plant_data = {}

    # retrieving CMAUP - plant id - CPD
    cmaup_plant = open("CMAUPv1.0_download_Plant_Ingredient_Associations_allIngredients.txt", "r")
    for line in cmaup_plant:
        line = line.replace("\n", "")
        data = line.split("\t")

        if data[1] in cmaup_plant_cpd.keys():
            cmaup_plant_cpd[data[1]].append(data[0])
        else:
            cmaup_plant_cpd[data[1]] = []
            cmaup_plant_cpd[data[1]].append(data[0])
          # key = CMAUP id, value = plant IDs
    cmaup_plant.close()


    # retrieving CMAUP names
    cmaup_mol = open("CMAUPv1.0_download_Ingredients_All.txt", "r")

    for line in cmaup_mol:
        line = line.replace("\n", "")
        data = line.split("\t")
        cmaup_data[data[0]] = data[1]  # key= cmaup id, value = cmaup name
    cmaup_mol.close()




    # retrieve plant data
    cmaup_plant_file = open("CMAUPv1.0_download_Plants.txt", "r")
    for line in cmaup_plant_file:
        if not line.startswith("Plant_ID"):
            line = line.replace("\n", "")
            data = line.split("\t")

            plant_id = data[0]
            plant_name = data[1]
            plant_tax_id = data[2]
            cmaup_plant_data[plant_id] = (plant_name, plant_tax_id)
    cmaup_plant_file.close()

    cocosource = open("COCONUT_w_source_cmaup.csv", "r")
    for line in cocosource:
        if not line.startswith("uniqueNaturalProduct"):
            line = line.replace("\n", "")
            data = line.split(",")

            coconut_id = data[0]
            cmaup_id = data[1]

            np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})

            np = np[0]

            # check if present in cmaup
            # add present in cmaup
            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {"found_in_databases": "cmaup"}})

            # UPDATE NAME
            if np["name"] is None or np["name"] == "":
                if cmaup_id in cmaup_data.keys() and cmaup_data[cmaup_id] != "":
                    # update name
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": cmaup_data[cmaup_id] }})
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"nameTrustLevel": 1}})
            elif len(np["name"]) > len(cmaup_data[cmaup_id]):
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"synonyms": np["name"]}})
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"name": cmaup_data[cmaup_id]}})
                db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$set": {"nameTrustLevel": 1}})

            # UPDATE species
            # in all cases the NP is produced by a plant
            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$pull": {"textTaxa": "notax"}})
            db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"textTaxa": "plants"}})

            if cmaup_id in cmaup_plant_cpd.keys():
                for plant_id in cmaup_plant_cpd[cmaup_id]:

                    if cmaup_plant_data[plant_id][0] != "" and cmaup_plant_data[plant_id][0] != "NA" :
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]},{"$addToSet": {"textTaxa": cmaup_plant_data[plant_id][0]}})
                    if cmaup_plant_data[plant_id][1] != "" and cmaup_plant_data[plant_id][1] != "NA":
                        db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {"taxid": cmaup_plant_data[plant_id][1]}})

            # update XREFS
            if "xrefs" in np.keys():
                cmaup_id_set = False
                for l in np["xrefs"]:
                    if l[0] == data[1]:
                        cmaup_id_set = True
                if not cmaup_id_set:
                    db.uniqueNaturalProduct.update_one({'coconut_id': data[0]}, {"$addToSet": {
                        "xrefs": [ "cmaup",data[1],"http://bidd2.nus.edu.sg/CMAUP/searchresults.php?keyword_search="]}})

            # update annotation level
            #annotationLevel = reevaluateAnnotationLevel(np, "cmaup")
            #db.uniqueNaturalProduct.update_one({'coconut_id': data[0]},{"$set": {"annotationLevel": annotationLevel}})
    cocosource.close()

    ##################################################################
    # updating with data from PubChem
    print("Curating from PubChem")
    #1	CNP0341011	Acetyl-DL-carnitine	3-acetyloxy-4-(trimethylazaniumyl)butanoate	['870-77-9', '14992-62-2']	['Acetyl-DL-carnitine', 'acetylcarnitine', 'DL-O-Acetylcarnitine', 'Acetyl carnitine', 'DL-Acetylcarnitine', '(+/-)-acetylcarnitine', '3-(acetyloxy)-4-(trimethylazaniumyl)butanoate', 'Ammonium, (3-carboxy-2-hydroxypropyl)trimethyl-, hydroxide, inner salt, acetate, DL-', '870-77-9', '14992-62-2', 'O-acetylcarnitine', '1-Propanaminium, 2-(acetoxy)-3-carboxy-N,N,N-trimethyl-, hydroxide, inner salt, (+-)- (9CI)', 'L-ACETYLCARNITINE', '1-Propanaminium, 2-(acetyloxy)-3-carboxy-N,N,N-trimethyl-, inner salt', 'CAR(2:0)', 'bmse000142', 'SCHEMBL69781', 'DTXSID2048117', 'CHEBI:73024', 'LMFA07070060', 'LS-17075', 'HY-126358', 'LS-173362', '3-(acetyloxy)-4-(trimethylammonio)butanoate', 'CS-0102945', 'Q27140241', '1-Propanaminium, 2-(acetoxy)-3-carboxy-N,N,N-trimethyl-, hydroxide, inner salt, (+-)-']
    cocosource = open("Pubchem_Data_For_COCONUT.txt", "r")
    for line in cocosource:
        line = line.replace("\n", "")
        data = line.split("\t")

        pubchem_id = data[0]
        coconut_id = data[1]
        name = data[2]
        iupac_name = data[3]

        if len(data)>4:
            list_of_cas = data[4]
        if len(data)>5:
            synonyms = data[5]

        np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})
        np = np[0]

        # check if present in pubchem
        # add present in pubchem
        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {"found_in_databases": "pubchem"}})

        #check and eventually update name
        if np["name"] is None or np["name"] == "" or nameIsWeird(np["name"]):
            if name != "" and len(name) < len(np["name"]) and not nameIsWeird(name):
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {"synonyms": np["name"]}})
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"name": name}})
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"nameTrustLevel": 1}})
            elif iupac_name != "" and len(iupac_name) < len(np["name"]) :
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {"synonyms": np["name"]}})
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"name": iupac_name}})
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"nameTrustLevel": 1}})

        #if synonyms != "":
        #    s = synonyms.replace("[", "")
        #    s = s.replace("]", "")
        #   s = s.replace("\'", "")
        #    s = s.split(", ")

            #for sy in s:
                #db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {"synonyms": sy}})

        if  list_of_cas != "":
            c = list_of_cas.replace("[", "")
            c = c.replace("]", "")
            c = c.replace("\'", "")
            c = c.split(", ")

            if 'cas' in np.keys():
                if  np['cas'] is None or np['cas'] == "" :
                    for cas in c:
                        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"cas": cas}})
            else:
                for cas in c:
                    db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"cas": cas}})

        # update XREFS
        if "xrefs" in np.keys():
            pubchem_id_set = False
            for l in np["xrefs"]:
                if l[1] == pubchem_id:
                    pubchem_id_set = True
            if not pubchem_id_set:
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$addToSet": {
                    "xrefs": ["pubchem_tested_np", pubchem_id,
                              "https://pubchem.ncbi.nlm.nih.gov/compound/"]}})

        #annotationLevel = reevaluateAnnotationLevel(np, "pubchem")
        #db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"annotationLevel": annotationLevel}})

    cocosource.close()

    print("done")
    ##################################################################
    # updating with ChEMBL









if __name__ == '__main__':
    main()





