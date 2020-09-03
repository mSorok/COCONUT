#!/usr/bin/python -tt

from pymongo import MongoClient
import re

from pymongo.errors import CursorNotFound


def main():

    client = MongoClient("localhost:27018")
    db = client['COCONUT2020-07']

    cocosource = open("coconut_iupac_name.txt", "r")

    inchikeyPattern = r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"

    casPattern = r"^[0-9]+-[0-9]+-[0-9]+^"


    # retry later:  bool(re.search(inchikeyPattern, np["name"])) and  bool(re.search(casPattern, np["name"]))

    for line in cocosource:
        # 7-methyl-9-{[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy}-1H,2H,3H,4H-cyclopenta[c]chromen-4-one	CNP0220816
        line = line.replace("\n", "")
        data = line.split("\t")
        iupac = data[0]
        coconut_id = data[1]

        np = db.uniqueNaturalProduct.find({"coconut_id": coconut_id})
        np = np[0]

        oriName = np["name"]

        db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"iupac_name": iupac}})

        if np["name"] == "" or len(np["name"]) <= 3 or np["name"].startswith("MCULE") or np["name"].startswith("STL") or "n.a." in np["name"].lower() or np["name"].startswith("STK") or np["name"].startswith("ZINC")\
                or "chembl" in np["name"].lower() or np["name"].startswith("CTK") or np["name"].lower().startswith("prestwick") or np["name"].lower().startswith("eim")\
                or np["name"].startswith("AKOS") or np["name"].lower().startswith("mls") or np["name"].startswith("ACMC") or np["name"].startswith("SCHEMBL") \
                or np["name"].startswith("NSC") or "dataset" in np["name"].lower() or np["name"].startswith("BMS") or np["name"].startswith("EINECS") or "named" in np["name"].lower() \
                or np["name"].startswith("NCI") or np["name"].startswith("HMS") or np["name"].lower().startswith("cbmicro") or np["name"].lower().startswith("oprea")  \
                or np["name"] == "ms 681c" or np["name"] == "tpu-0037-c" or np["name"].lower().startswith("cbdive") or np["name"] == "J3.542.783J" or np["name"] == "F0848-0042" \
                or np["name"] == "CK-636" or "megapho" in np["name"].lower() or np["name"].startswith("sb 253") or np["name"].startswith("CCG-") or "nadp" in np["name"].lower() \
                or np["name"] == "jbir-28" or np["name"] == "F0848-0042" or np["name"].startswith("DA-") or "glas#" in np["name"] or np["name"] == "jbir-67" or np["name"] == "wss2222" \
                or np["name"].startswith("NCGC") or np["name"].startswith("MFCD") or np["name"].startswith("SDCCG") or "refreshn" in np["name"].lower() or "pharmaco" in np["name"].lower() \
                or np["name"].startswith("BRD") or np["name"].startswith("BIM") or np["name"].startswith("NINDS") or np["name"].startswith("CFF") or "ampule" in np["name"] \
                or np["name"].lower().startswith("bio1") or np["name"].lower().startswith("kbio") or np["name"].startswith("F3371") or np["name"].startswith("PDSP") \
                or "testing" in np["name"].lower() or "specifications" in np["name"].lower() or "melting" in np["name"].lower() or "point" in np["name"].lower() \
                or np["name"].startswith("UNII") or "EP2" in np["name"] or "thermosystem" in np["name"] or "component" in np["name"] or np["name"].startswith("CCRIS") \
                or np["name"].startswith("BDBM") or np["name"].startswith("ST0") or np["name"] == "C 0750" or "pharmaceutical" in np["name"].lower() or np["name"].startswith("EU-") \
                or np["name"].startswith("BIDD") or  "99%" in np["name"] or "100%" in np["name"] or "98%" in np["name"] or np["name"].startswith("DB-") or np["name"] == "No Doz" \
                or "special grade" in np["name"] or np["name"].startswith("PubChem") or np["name"].startswith("BRN") or np["name"].startswith("InChI") or np["name"].startswith("HSDB") \
                or "vial of" in np["name"] or np["name"].startswith("bmse") or "tested" in np["name"] or np["name"].startswith("ARONI") or np["name"].startswith("L00") \
                or "einecs" in np["name"].lower() or np["name"].startswith("IDI1") or np["name"].startswith("EC ") or np["name"].startswith("07E4") or np["name"].startswith("BSP") \
                or np["name"].startswith("SBI") or "chebi" in np["name"].lower() or "reference" in np["name"].lower() or np["name"].startswith("AB0") or "spectrum" in np["name"].lower() \
                or "compound" in np["name"].lower() or np["name"].startswith("AS-") or np["name"] == "O926" or "synthetic" in np["name"].lower() or "specified" in np["name"].lower() \
                or np["name"].startswith("ACN") or "molmap" in np["name"].lower() or "HPLC" in np["name"] or np["name"].startswith("CS-") or "tracecert" in np["name"].lower() \
                or np["name"].startswith("AI3") or "lopac" in np["name"].lower() or np["name"] == "1l5q" or np["name"].startswith("SR-") or np["name"].startswith("DSS") or "microg/mL" in np["name"].lower() \
                or "bayer" in np["name"].lower() or np["name"] == "1l7x" or np["name"].lower().startswith("divk") or np["name"].lower().startswith("dtxs") or "production" in np["name"].lower() \
                or "BAN:JAN" in np["name"] or "WLN:" in np["name"] or "mg/ml" in np["name"].lower() or "probes" in np["name"].lower() or "pharmagrade" in np["name"].lower() or np["name"] == "FCC" \
                or "reagent" in np["name"].lower() or np["name"] == "H2815" or np["name"].startswith("GTPL") or "standard" in np["name"].lower() or np["name"].startswith("SBB") \
                or np["name"] == "Nix Nap" or "powder" in np["name"].lower() or np["name"].startswith("Tox2") or "tri-aqu" in np["name"].lower() or "JP17" in np["name"] or np["name"].startswith("CU-") \
                or "bioxtra" in np["name"].lower() or np["name"].endswith("-") or "caswell" in np["name"].lower() or np["name"].startswith("SMR") or np["name"].startswith("LS-") \
                or "salt" in np["name"].lower():


            db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"name": iupac}})

            if coconut_id == "CNP0248861" or coconut_id == "CNP0257017" or coconut_id == "CNP0228556":
                print("Changed "+ coconut_id +" name to IUPAC")
            #print("changed "+coconut_id)


        else:
            isFormula = False
            isCAS = False
            isCID = False
            isInchiKey = False
            isBrokenInchi = False
            isWeird = False

            if re.match(r"C[0-9]+?H[0-9].+", np["name"]):
                isFormula = True
            if re.match(r"^[0-9]+-[0-9]+-[0-9]+$", np["name"]):
                isCAS = True
            if re.match(r"^cid[0-9]+$", np["name"]) or re.match(r"^GE[0-9]+$", np["name"]):
                isCID = True
            if re.match(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$", np["name"]):
                isInchiKey = True
            if re.match(r"\d\(\d+z\)|\d\(\d+z,\d+z\)", np["name"])  or np["name"].startswith("pc(") or np["name"].startswith("dg("):
                isBrokenInchi=True
            if re.match(r"[0-9]+:[0-9]+\/", np["name"]):
                isWeird = True

            a_boolean_list = [isFormula,isCAS , isCID, isInchiKey,isBrokenInchi,isWeird  ]



            if sum(a_boolean_list) != 0:

                if coconut_id == "CNP0248861":
                    print("Detected name of CNP0248861 as weird " + str(a_boolean_list) + " "+ np["name"])

                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"name": iupac}})

            if isCAS:
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"cas": oriName}})

        for synonym in np["synonyms"]:
            if synonym == "" or len(synonym) <= 3 or synonym.startswith("MCULE") or synonym.startswith("STL") or "n.a." in synonym.lower() or synonym.startswith("STK") or synonym.startswith("ZINC")\
                or "chembl" in synonym.lower() or synonym.startswith("CTK") or synonym.lower().startswith("prestwick") or synonym.lower().startswith("eim")\
                or synonym.startswith("AKOS") or synonym.lower().startswith("mls") or synonym.startswith("ACMC") or synonym.startswith("SCHEMBL") \
                or synonym.startswith("NSC") or "dataset" in synonym.lower() or synonym.startswith("BMS") or synonym.startswith("EINECS") or "named" in synonym.lower() \
                or synonym.startswith("NCI") or synonym.startswith("HMS") or synonym.lower().startswith("cbmicro") or synonym.lower().startswith("oprea")  \
                or synonym == "ms 681c" or synonym == "tpu-0037-c" or synonym.lower().startswith("cbdive") or synonym == "J3.542.783J" or synonym == "F0848-0042" \
                or synonym == "CK-636" or "megapho" in synonym.lower() or synonym.startswith("sb 253") or synonym.startswith("CCG-") or "nadp" in synonym.lower() \
                or synonym == "jbir-28" or synonym == "F0848-0042" or synonym.startswith("DA-") or "glas#" in synonym or synonym == "jbir-67" or synonym == "wss2222" \
                or synonym.startswith("NCGC") or synonym.startswith("MFCD") or synonym.startswith("SDCCG") or "refreshn" in synonym.lower() or "pharmaco" in synonym.lower() \
                or synonym.startswith("BRD") or synonym.startswith("BIM") or synonym.startswith("NINDS") or synonym.startswith("CFF") or "ampule" in synonym \
                or synonym.lower().startswith("bio1") or synonym.lower().startswith("kbio") or synonym.startswith("F3371") or synonym.startswith("PDSP") \
                or "testing" in synonym.lower() or "specifications" in synonym.lower() or "melting" in synonym.lower() or "point" in synonym.lower() \
                or synonym.startswith("UNII") or "EP2" in synonym or "thermosystem" in synonym or "component" in synonym or synonym.startswith("CCRIS") \
                or synonym.startswith("BDBM") or synonym.startswith("ST0") or synonym == "C 0750" or "pharmaceutical" in synonym.lower() or synonym.startswith("EU-") \
                or synonym.startswith("BIDD") or  "99%" in synonym or "100%" in synonym or "98%" in synonym or synonym.startswith("DB-") or synonym == "No Doz" \
                or "special grade" in synonym or synonym.startswith("PubChem") or synonym.startswith("BRN") or synonym.startswith("InChI") or synonym.startswith("HSDB") \
                or "vial of" in synonym or synonym.startswith("bmse") or "tested" in synonym or synonym.startswith("ARONI") or synonym.startswith("L00") \
                or "einecs" in synonym.lower() or synonym.startswith("IDI1") or synonym.startswith("EC ") or synonym.startswith("07E4") or synonym.startswith("BSP") \
                or synonym.startswith("SBI") or "chebi" in synonym.lower() or "reference" in synonym.lower() or synonym.startswith("AB0") or "spectrum" in synonym.lower() \
                or "compound" in synonym.lower() or synonym.startswith("AS-") or synonym == "O926" or "synthetic" in synonym.lower() or "specified" in synonym.lower() \
                or synonym.startswith("ACN") or "molmap" in synonym.lower() or "HPLC" in synonym or synonym.startswith("CS-") or "tracecert" in synonym.lower() \
                or synonym.startswith("AI3") or "lopac" in synonym.lower() or synonym == "1l5q" or synonym.startswith("SR-") or synonym.startswith("DSS") or "microg/mL" in synonym.lower() \
                or "bayer" in synonym.lower() or synonym == "1l7x" or synonym.lower().startswith("divk") or synonym.lower().startswith("dtxs") or "production" in synonym.lower() \
                or "BAN:JAN" in synonym or "WLN:" in synonym or "mg/ml" in synonym.lower() or "probes" in synonym.lower() or "pharmagrade" in synonym.lower() or synonym == "FCC" \
                or "reagent" in synonym.lower() or synonym == "H2815" or synonym.startswith("GTPL") or "standard" in synonym.lower() or synonym.startswith("SBB") \
                or synonym == "Nix Nap" or "powder" in synonym.lower() or synonym.startswith("Tox2") or "tri-aqu" in synonym.lower() or "JP17" in synonym or synonym.startswith("CU-") \
                or "bioxtra" in synonym.lower() or synonym.endswith("-") or "caswell" in synonym.lower() or synonym.startswith("SMR") or synonym.startswith("LS-") \
                or "salt" in synonym.lower():


                    #print("found some crap: "+synonym)
                    db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$pull": {"synonyms": synonym}})





    cocosource.close()

    client.close()

    print("done")


    print("names to title")

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

            if 'iupac_name' in np.keys() and np['name'] != np['iupac_name']:
                db.uniqueNaturalProduct.update_one({'coconut_id': coconut_id}, {"$set": {"name":np['name'].title() }})





            #sleep(0.05)

    except CursorNotFound:
        pass

    print("all done")




if __name__ == '__main__':
    main()

