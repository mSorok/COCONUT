from chembl_structure_pipeline import standardizer
import sys
from rdkit import Chem


def main():
    smiles_to_uncharge = sys.argv[1]
    mol = Chem.MolFromSmiles(smiles_to_uncharge)

    pmol = standardizer.get_parent_mol(mol, verbose=False)
    print(Chem.MolToSmiles(pmol[0]))



if __name__ == '__main__':
    main()

