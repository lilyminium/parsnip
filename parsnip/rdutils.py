from rdkit import Chem
import numpy as np

def rdview(mol):
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def get_pdbinfo(atom):
    info = atom.GetPDBResidueInfo()

    #TODO: should I copy the aton and/or info?

    if info is None:
        info = Chem.AtomPDBResidueInfo()
    return info

def get_neighbor_ints(rdmol, indices, n_neighbors=3):
    neighbor_ints = set()
    for ix in indices:
        current_neighbors = n_neighbors
        latest_indices = set(map(int, indices))

        while current_neighbors > 0:
            current_indices = set()
            for latest_ix in latest_indices:
                atom = rdmol.GetAtomWithIdx(latest_ix)
                for nb in atom.GetNeighbors():
                    nbix = nb.GetIdx()
                    if nbix not in latest_indices and nbix not in indices:
                        current_indices.add(nbix)
            neighbor_ints |= current_indices
            latest_indices = current_indices
            current_neighbors -= 1
    
    return sorted(map(int, list(neighbor_ints)))


def copy_pdbinfo(info):
    new_info = Chem.AtomPDBResidueInfo()
    for attr in ("Name", "SerialNumber", "AltLoc", "ResidueName", "ResidueNumber",
                 "ChainId", "InsertionCode", "Occupancy", "TempFactor", "IsHeteroAtom",
                 "SecondaryStructure", "SegmentNumber"):
        old = getattr(info, f"Get{attr}")()
        getattr(new_info, f"Set{attr}")(old)
    # great
    return new_info