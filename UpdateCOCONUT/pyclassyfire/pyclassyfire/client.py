"""A client for the ClassyFire API which enables efficient querying with
 chemical database files"""

import csv
import io
import json
import os
import time

import requests

url = "http://classyfire.wishartlab.com"
chunk_size = 1000
sleep_interval = 60


def structure_query(compound, label='pyclassyfire'):
    """Submit a compound information to the ClassyFire service for evaluation 
    and receive a id which can be used to used to collect results
    
    :param compound: The compound structures as line delimited inchikey or
         smiles. Optionally a tab-separated id may be prepended for each
         structure.
    :type compound: str
    :param label: A label for the query
    :type label:
    :return: A query ID number
    :rtype: int
    
    >>> structure_query('CCC', 'smiles_test')
    >>> structure_query('InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)')
    
    """
    r = requests.post(url + '/queries.json', data='{"label": "%s", '
                      '"query_input": "%s", "query_type": "STRUCTURE"}'
                                                  % (label, compound),
                      headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()['id']


def iupac_query(compound, label='pyclassyfire'):
    """Submit a IUPAC compound name to the ClassyFire service for evaluation
     and receive a id which can be used to used to collect results.
    
    :param compound: The line delimited compound names. Optionally a
         tab-separated id may be prepended for each compound.
    :type compound: str
    :param label: A label for the query
    :type label:
    :return: A query ID number
    :rtype: int
    
    >>> iupac_query('ethane', 'iupac_test')
    >>> iupac_query('C001\\tethane\\nC002\\tethanol', 'iupac_test')
    
    """
    r = requests.post(url + '/queries.json', data='{"label": "%s", '
                      '"query_input": "%s", "query_type": "IUPAC_NAME"}'
                                                  % (label, compound),
                      headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()['id']


def get_results(query_id, return_format="json"):
    """Given a query_id, fetch the classification results.
    
    :param query_id: A numeric query id returned at time of query submission
    :type query_id: str
    :param return_format: desired return format. valid types are json, csv or sdf
    :type return_format: str
    :return: query information
    :rtype: str
    
    >>> get_results('595535', 'csv')
    >>> get_results('595535', 'json')
    >>> get_results('595535', 'sdf')
    
    """
    r = requests.get('%s/queries/%s.%s' % (url, query_id, return_format),
                     headers={"Content-Type": "application/%s" % return_format})
    r.raise_for_status()
    return r.text


def get_entity(inchikey, return_format="json"):
    """Given a InChIKey for a previously queried structure, fetch the
     classification results.

    :param inchikey: An InChIKey for a previously calculated chemical structure
    :type inchikey: str
    :param return_format: desired return format. valid types are json, csv or sdf
    :type return_format: str
    :return: query information
    :rtype: str
    
    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'csv')
    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'json')
    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'sdf')
    
    """
    inchikey = inchikey.replace('InChIKey=', '')
    r = requests.get('%s/entities/%s.%s' % (url, inchikey, return_format),
                     headers={
                         "Content-Type": "application/%s" % return_format})
    r.raise_for_status()
    return r.text


def get_chemont_node(chemontid):
    """Return data for the TaxNode with ID chemontid.
    
    :param chemontid: the ChemOnt ID of the entity.
    :type chemontid: str
    :return: The classification results for the entity as json.
    :rtype: str
    
    >>> get_chemont_node('CHEMONTID:0004253')
    
    """
    chemontid = chemontid.replace("CHEMONTID:", "C")
    r = requests.get('%s/tax_nodes/%s.json' % (url, chemontid),
                     headers={"Content-Type": "application/json" })
    r.raise_for_status()
    return r.text


def tabular_query(inpath, structure_key, dialect='excel', outpath=None,
                  outfields=('taxonomy', 'description', 'substituents')):
    """Given a path to a compound set in tabular form (comma or tab delimited)
     annotate all compounds and write results to an expanded table.
    
    :param inpath: path to compound file to be annotated
    :type inpath: str
    :param structure_key: column heading which contains the compounds InChIKey
         or SMILES
    :type structure_key: str
    :param dialect: dialect for parsing table (generally 'excel' for csv,
         'excel-tab' for tsv)
    :type dialect: str
    :param outpath: Path to desired output location
    :type outpath: str
    :param outfields: Fields to append to table from ClassyFire output
    :type outfields: tuple(string)
    
    >>> tabular_query('/tabulated_data.tsv', 'structure', 'excel-tab')
    
    """
    tax_fields = ('kingdom', 'superclass', 'class', 'subclass')
    query_ids = []
    infile = open(inpath, 'rU')
    if not outpath:
        outpath = _prevent_overwrite(inpath)
    comps = []
    for line in csv.DictReader(infile, dialect=dialect):
        comps.append(line[structure_key])
        if not len(comps) % chunk_size:
            query_ids.append(structure_query('\\n'.join(comps)))
            comps = []
    if comps:
        query_ids.append(structure_query('\\n'.join(comps)))
    print('%s queries submitted to ClassyFire API' % len(query_ids))
    i = 0
    infile.seek(0)
    with open(outpath, 'w') as outfile:
        reader = csv.DictReader(infile, dialect=dialect)
        writer = csv.DictWriter(outfile, reader.fieldnames+list(outfields),
                                dialect=dialect)
        writer.writeheader()
        while i < len(query_ids):
            result = json.loads(get_results(query_ids[i]))
            if result["classification_status"] == "Done":
                for j, line in enumerate(reader):
                    if result['entities'] and str(j+1) == result['entities'][0]['identifier'].split('-')[1]:
                        hit = result['entities'].pop(0)
                        if 'taxonomy' in outfields:
                            hit['taxonomy'] = ";".join(
                                ['%s:%s' % (hit[x]['name'], hit[x]['chemont_id'])
                                 for x in tax_fields if hit[x]])
                        for field in outfields:
                            if isinstance(hit[field], list):
                                line[field] = ';'.join(hit[field])
                            else:
                                line[field] = hit[field]
                    writer.writerow(line)
                i += 1
            else:
                print("%s percent complete" % round(i/len(query_ids)*100))
                time.sleep(sleep_interval)
    infile.close()


def sdf_query(inpath, outpath=None):
    """Given a path to a compound set in a sdf file, annotate all compounds
     and write results as attributes in a sdf file.
    
    :param inpath: path to compound file to be annotated
    :type inpath: str
    :param outpath: Path to desired output location
    :type outpath: str

    >>> sdf_query('/sdf_data.sdf')
    
    """
    from rdkit.Chem import AllChem
    query_ids = []
    if not outpath:
        outpath = _prevent_overwrite(inpath)
    comps = []
    for mol in AllChem.SDMolSupplier(inpath):
        if mol:
            comps.append(AllChem.MolToSmiles(mol))
        if not len(comps) % chunk_size:
            query_ids.append(structure_query('/n'.join(comps)))
            comps = []
    if comps:
        query_ids.append(structure_query('\\n'.join(comps)))
    print('%s queries submitted to ClassyFire API' % len(query_ids))
    i = 0
    with io.open(outpath, 'w', encoding="utf-8") as outfile:
        while i < len(query_ids):
            result = json.loads(get_results(query_ids[i]))
            if result["classification_status"] == "Done":
                outfile.write(get_results(query_ids[i], return_format='sdf'))
                i += 1
            else:
                print("%s percent complete" % round(i / len(query_ids) * 100))
                time.sleep(sleep_interval)


def _prevent_overwrite(write_path, suffix='_annotated'):
    """Prevents overwrite of existing output files by appending a suffix when
    needed

    :param write_path: potential write path
    :type write_path: string
    :return:
    :rtype:
    """
    while os.path.exists(write_path):
        sp = write_path.split('.')
        if len(sp) > 1:
            sp[-2] += suffix
            write_path = '.'.join(sp)
        else:
            write_path += suffix
    return write_path
