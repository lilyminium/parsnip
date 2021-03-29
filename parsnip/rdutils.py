from rdkit import Chem

def rdview(mol):
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def get_neighbor_ints(rdmol, indices, n_neighbors=3):
    neighbor_ints = set()
    for ix in indices:
        current_neighbors = n_neighbors
        latest_indices = set(indices)

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
    
    return sorted(neighbor_ints)


