from pyclassyfire import client
import os
import filecmp

data_dir = os.path.dirname(__file__)+'/test_data'
query = {'csv': 'CompoundID,ChemOntID,ParentName\nQ595535-1,CHEMONTID:0002500,Alkanes\n',
         'json': '{"id":595535,"label":"curl_test","classification_status":"Done","number_of_elements":1,"number_of_pages":1,"invalid_entities":[],"entities":[{"identifier":"Q595535-1","smiles":"CCC","inchikey":"InChIKey=ATUOYWHBWRKTHZ-UHFFFAOYSA-N","kingdom":{"name":"Organic compounds","description":"Compounds that contain at least one carbon atom, excluding isocyanide/cyanide and their non-hydrocarbyl derivatives, thiophosgene, carbon diselenide, carbon monosulfide, carbon disulfide, carbon subsulfide, carbon monoxide, carbon dioxide, Carbon suboxide, and dicarbon monoxide.","chemont_id":"CHEMONTID:0000000","url":"http://classyfire.wishartlab.com/tax_nodes/C0000000"},"superclass":{"name":"Hydrocarbons","description":"Organic compounds made up only of carbon and hydrogen atoms.","chemont_id":"CHEMONTID:0002837","url":"http://classyfire.wishartlab.com/tax_nodes/C0002837"},"class":{"name":"Saturated hydrocarbons","description":"Hydrocarbons that contains only saturated carbon atoms, which are linked to one another through single bonds. These includes alkanes and cycloalkanes.","chemont_id":"CHEMONTID:0004474","url":"http://classyfire.wishartlab.com/tax_nodes/C0004474"},"subclass":{"name":"Alkanes","description":"Acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","chemont_id":"CHEMONTID:0002500","url":"http://classyfire.wishartlab.com/tax_nodes/C0002500"},"intermediate_nodes":[],"direct_parent":{"name":"Alkanes","description":"Acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","chemont_id":"CHEMONTID:0002500","url":"http://classyfire.wishartlab.com/tax_nodes/C0002500"},"alternative_parents":[],"molecular_framework":"Aliphatic acyclic compounds","substituents":["Acyclic alkane","Alkane","Aliphatic acyclic compound"],"description":"This compound belongs to the class of organic compounds known as alkanes. These are acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","external_descriptors":[{"source":"CHEBI","source_id":"CHEBI:32879","annotations":["alkane"]}],"ancestors":["Alkanes","Chemical entities","Hydrocarbons","Organic compounds","Saturated hydrocarbons"],"predicted_chebi_terms":["chemical entity (CHEBI:24431)","hydrocarbon (CHEBI:24632)","alkane (CHEBI:18310)","organic molecule (CHEBI:72695)"],"predicted_lipidmaps_terms":[],"classification_version":"2.1"}]}',
         'sdf': 'Q595535-1\nMrv1568 12281606022D          \n\n  3  2  0  0  0  0            999 V2000\n    1.2375   -0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9520   -1.1270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6664   -0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n> <InChIKey>\nInChIKey=ATUOYWHBWRKTHZ-UHFFFAOYSA-N\n> <SMILES>\nCCC\n> <Kingdom>\nOrganic compounds\n> <Superclass>\nHydrocarbons\n> <Class>\nSaturated hydrocarbons\n> <Subclass>\nAlkanes\n> <Intermediate Nodes>\n\n> <Direct Parent>\nAlkanes\n> <Alternative Parents>\n\n> <Molecular Framework>\nAliphatic acyclic compounds\n> <Substituents>\nAcyclic alkane\tAlkane\tAliphatic acyclic compound\n> <Structure-based description>\nThis compound belongs to the class of organic compounds known as alkanes. These are acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.\n> <Ancestors>\nAlkanes\tChemical entities\tHydrocarbons\tOrganic compounds\tSaturated hydrocarbons\n> <External Descriptors>\nalkane (CHEBI, CHEBI:32879)\n$$$$\n'
         }
entity = {'csv': 'ChemOntID,ParentName\nCHEMONTID:0002500,Alkanes\n',
          'json': '{"smiles":"CCC","inchikey":"InChIKey=ATUOYWHBWRKTHZ-UHFFFAOYSA-N","kingdom":{"name":"Organic compounds","description":"Compounds that contain at least one carbon atom, excluding isocyanide/cyanide and their non-hydrocarbyl derivatives, thiophosgene, carbon diselenide, carbon monosulfide, carbon disulfide, carbon subsulfide, carbon monoxide, carbon dioxide, Carbon suboxide, and dicarbon monoxide.","chemont_id":"CHEMONTID:0000000","url":"http://classyfire.wishartlab.com/tax_nodes/C0000000"},"superclass":{"name":"Hydrocarbons","description":"Organic compounds made up only of carbon and hydrogen atoms.","chemont_id":"CHEMONTID:0002837","url":"http://classyfire.wishartlab.com/tax_nodes/C0002837"},"class":{"name":"Saturated hydrocarbons","description":"Hydrocarbons that contains only saturated carbon atoms, which are linked to one another through single bonds. These includes alkanes and cycloalkanes.","chemont_id":"CHEMONTID:0004474","url":"http://classyfire.wishartlab.com/tax_nodes/C0004474"},"subclass":{"name":"Alkanes","description":"Acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","chemont_id":"CHEMONTID:0002500","url":"http://classyfire.wishartlab.com/tax_nodes/C0002500"},"intermediate_nodes":[],"direct_parent":{"name":"Alkanes","description":"Acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","chemont_id":"CHEMONTID:0002500","url":"http://classyfire.wishartlab.com/tax_nodes/C0002500"},"alternative_parents":[],"molecular_framework":"Aliphatic acyclic compounds","substituents":["Acyclic alkane","Alkane","Aliphatic acyclic compound"],"description":"This compound belongs to the class of organic compounds known as alkanes. These are acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.","external_descriptors":[{"source":"CHEBI","source_id":"CHEBI:32879","annotations":["alkane"]}],"ancestors":["Alkanes","Chemical entities","Hydrocarbons","Organic compounds","Saturated hydrocarbons"],"predicted_chebi_terms":["alkane (CHEBI:18310)","chemical entity (CHEBI:24431)","organic molecule (CHEBI:72695)","hydrocarbon (CHEBI:24632)"],"predicted_lipidmaps_terms":[],"classification_version":"2.1"}',
          'sdf': 'CCC\nMrv1568 12281606022D          \n\n  3  2  0  0  0  0            999 V2000\n    1.2375   -0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9520   -1.1270    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6664   -0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n> <InChIKey>\nInChIKey=ATUOYWHBWRKTHZ-UHFFFAOYSA-N\n> <SMILES>\nCCC\n> <Kingdom>\nOrganic compounds\n> <Superclass>\nHydrocarbons\n> <Class>\nSaturated hydrocarbons\n> <Subclass>\nAlkanes\n> <Intermediate Nodes>\n\n> <Direct Parent>\nAlkanes\n> <Alternative Parents>\n\n> <Molecular Framework>\nAliphatic acyclic compounds\n> <Substituents>\nAcyclic alkane\tAlkane\tAliphatic acyclic compound\n> <Structure-based description>\nThis compound belongs to the class of organic compounds known as alkanes. These are acyclic branched or unbranched hydrocarbons having the general formula CnH2n+2 , and therefore consisting entirely of hydrogen atoms and saturated carbon atoms.\n> <Ancestors>\nAlkanes\tChemical entities\tHydrocarbons\tOrganic compounds\tSaturated hydrocarbons\n> <External Descriptors>\nalkane (CHEBI, CHEBI:32879)\n$$$$\n'
          }
annotated_line = "cpd00002	atp	ATP	C10H13N5O13P3	504	ModelSEED	[H][C@]1(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O[C@@]([H])(N2C=NC3=C(N)N=CN=C23)[C@]([H])(O)[C@]1([H])O	-3	1	0	null	0	-673.85	3.04314	26:0.84;22:2.95;14:13.03;29:1.83;30:7.72	15:3.97;14:-3.66;6:-3.48;9:-9.18	null	null	null	Organic compounds:CHEMONTID:0000000;Nucleosides, nucleotides, and analogues:CHEMONTID:0000289;Purine nucleotides:CHEMONTID:0001506;Purine ribonucleotides:CHEMONTID:0001544	This compound belongs to the class of organic compounds known as purine ribonucleoside triphosphates. These are purine ribobucleotides with a triphosphate group  linked to the ribose moiety.	Purine ribonucleoside triphosphate;Purine ribonucleoside monophosphate;Pentose phosphate;Pentose-5-phosphate;Glycosyl compound;N-glycosyl compound;6-aminopurine;Monosaccharide phosphate;Pentose monosaccharide;Imidazopyrimidine;Purine;Aminopyrimidine;Monoalkyl phosphate;Monosaccharide;N-substituted imidazole;Organic phosphoric acid derivative;Phosphoric acid ester;Imidolactam;Alkyl phosphate;Pyrimidine;Azole;Tetrahydrofuran;Imidazole;Heteroaromatic compound;Secondary alcohol;1,2-diol;Organoheterocyclic compound;Azacycle;Oxacycle;Organooxygen compound;Hydrocarbon derivative;Organic nitrogen compound;Organic oxide;Organopnictogen compound;Amine;Primary amine;Organic oxygen compound;Alcohol;Organonitrogen compound;Aromatic heteropolycyclic compound"
chemont_node = '{"name":"2-arylethylamines","chemont_id":"CHEMONTID:0004253","description":"Primary amines that have the general formula RCCNH2, where R is an organic group.","url":"http://classyfire.wishartlab.com/tax_nodes/C0004253","parent":"Primary amines","nr_of_entities":127868,"synonyms":[{"name":"2-arylethylamine","source":"CHEBI","source_id":"CHEBI:55436","mapping_scope":"Related"}]}'


def test_structure_query():
    assert isinstance(client.structure_query('CCC', 'smiles_test'), int)
    assert isinstance(client.structure_query(
        'InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)'), int)


def test_iupac_query():
    assert isinstance(client.iupac_query('ethane', 'iupac_test'), int)
    assert isinstance(client.iupac_query('C001\\tethane\\nC002\\tethanol',
                                         'iupac_test'), int)


def test_get_results():
    assert client.get_results('595535', 'csv') == query['csv']
    assert client.get_results('595535', 'json') == query['json']
    assert client.get_results('595535', 'sdf') == query['sdf']


def test_get_entity():
    inchikey = "ATUOYWHBWRKTHZ-UHFFFAOYSA-N"
    assert client.get_entity(inchikey, 'csv') == entity['csv']
    assert client.get_entity(inchikey, 'json') == entity['json']
    #assert client.get_entity(inchikey, 'sdf') == entity['sdf']


def test_get_chemont_node():
    assert client.get_chemont_node('CHEMONTID:0004253') == chemont_node


def test_tabular_query():
    try:
        client.tabular_query(data_dir+'/tabulated_data.tsv', 'structure',
                             'excel-tab')
        assert annotated_line in open(data_dir + '/tabulated_data_annotated.tsv').read()
        # assert filecmp.cmp(data_dir+'/tabulated_data_annotated.tsv',
        #                  data_dir+'/annotated.tsv')
    finally:
        os.remove(data_dir+'/tabulated_data_annotated.tsv')


def test_sdf_query():
    try:
        client.sdf_query(data_dir+'/sdf_data.sdf')
        assert open(data_dir+'/sdf_data_annotated.sdf').readlines()
    finally:
        os.remove(data_dir + '/sdf_data_annotated.sdf')
