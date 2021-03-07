from rdkit.Chem import MolFromPDBFile


def pdb2rdmol(pdb):
    return MolFromPDBFile(pdb, False, False)