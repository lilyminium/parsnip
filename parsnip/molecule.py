
import itertools
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem
import nglview as nv
import networkx as nx
import numpy as np
import MDAnalysis as mda

from .base import ParamList, ArrayListWithIndices
from .tags import Tag, TagDict
from .atom import Atom
from .params import Bond, Angle, Dihedral

class Molecule:
    """Annoying class telling users to please don't
    mess with the attributes.

    I use schmancy List/Dict subclasses to deal with indexing and stuff.
    
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._atoms = ArrayListWithIndices()
        self._bonds = ParamList()
        self._angles = ParamList()
        self._dihedrals = ParamList()
        self._impropers = ParamList()
        self._pairs = ParamList()
        self._exclusions = ParamList()
        self._tags = TagDict()

    @property
    def atoms(self):
        return self._atoms
    
    @atoms.setter
    def atoms(self, iterable):
        self._atoms = ArrayListWithIndices(iterable)

    @property
    def bonds(self):
        return self._bonds

    @bonds.setter
    def bonds(self, iterable):
        self._bonds = ParamList(iterable)
    
    @property
    def angles(self):
        return self._angles

    @angles.setter
    def angles(self, iterable):
        self._angles = ParamList(iterable)

    @property
    def dihedrals(self):
        return self._dihedrals

    @dihedrals.setter
    def dihedrals(self, iterable):
        self._dihedrals = ParamList(iterable)

    @property
    def impropers(self):
        return self._impropers

    @impropers.setter
    def impropers(self, iterable):
        self._impropers = ParamList(iterable)

    @property
    def pairs(self):
        return self._pairs

    @pairs.setter
    def pairs(self, iterable):
        self._pairs = ParamList(iterable)

    @property
    def exclusions(self):
        return self._exclusions

    @exclusions.setter
    def exclusions(self, iterable):
        self._exclusions = ParamList(iterable)

    @property
    def tags(self):
        return self._tags
    
    @tags.setter
    def tags(self, iterable):
        self._tags = TagDict(iterable)
    
    @property
    def params(self):
        return [self.bonds, self.angles, self.dihedrals,
                self.impropers, self.pairs, self.exclusions]

    @property
    def n_atoms(self):
        return len(self.atoms)

    def remove_atoms(self, atoms):
        mask = ~np.isin(self.atoms, atoms)
        self.atoms = self.atoms[mask]

    def remove_params_within_atoms(self, atoms):
        arr = np.array(atoms)
        for params in self.params:
            params.remove_within_atoms(self.atoms)
    
    def clean(self):
        self.clean_params()
        self.clean_tags()

    def clean_params(self):
        for params in self.params:
            params.keep_within_atoms(self.atoms)

    def clean_tags(self):
        self.tags.keep_within_atoms(self.atoms)
    
    def get_tag(self, tag_name: str, tag_index: int=-1):
        return self.tags.get_tag(tag_name, tag_index)

    def add_tag(self, tag_name: str, atoms=[]):
        tag = Tag(tag_name, atoms)
        self.tags[tag_name].append(tag)

    def add_tag_from_atom_indices(self, tag_name: str, indices: list=[]):
        atoms = self.atoms[indices]
        self.add_tag(tag_name, atoms)
    
    def show_tag(self, name, index=-1):
        indices = self.get_tag(name, index).indices
        view = self.view()
        view.add_representation("ball+stick", indices,
                                opacity=0.3, color="orange",
                                radius=0.3)
        return view

    def reindex_atoms(self):
        for i, atom in enumerate(self.atoms):
            atoms.index = i
    



class RDMolecule(Molecule):

    """Sigh don't mess with this stuff either"""

    def __init__(self, rdmol, *args, **kwargs):
        super().__init__()
        self.rdmol = rdmol
    
    @property
    def rdmol(self):
        return self._rdmol

    @rdmol.setter
    def rdmol(self, mol):
        rdmol = Chem.AddHs(rdmol, explicitOnly=True, addCoords=True)
        self._rdmol = Chem.rdchem.RWMol(rdmol)
        self.atoms.clear()
        for params in self.params:
            params.clear()
        self.tags.clear()

        for i, atom in enumerate(self.rdmol.GetAtoms()):
            at = Atom(atom, name=atom.GetName())
            self.atoms.append(at)

        for i, bond in enumerate(self.rdmol.GetBonds()):
            j = bond.GetBeginAtomIdx()
            k = bond.GetEndAtomIdx()
            bd = Bond(self.atoms[[j, k]])
            self.bonds.append(bd)
        
        # infer angles, dihedrals from graph
        graph = nx.Graph()
        graph.add_edges_from(self.bonds.indices)
        # TODO: is there a less O(n^2) and iffy way to do this
        pairs = itertools.combinations(graph.nodes, 2)
        seen = []
        for a, b in pairs:
            paths = nx.all_simple_paths(graph, source=a, target=b, cutoff=4)
            paths = [x for x in paths if len(paths) > 2]
            for path in paths:
                if path in seen or path[::-1] in seen:
                    continue
                if len(path) == 3:
                    self.angles.append(Angle(self.atoms[path]))
                elif len(path) == 4:
                    self.dihedrals.append(Dihedral(self.atoms[path]))
                seen.append(path)

    def to_mda_universe(self):
        if not self.rdmol.GetNumAtoms():
            return mda.Universe.empty(0)
        return mda.Universe(self.rdmol, format="RDKIT")

    def write(self, filename):
        u = self.to_mda_universe()
        u.atoms.write(filename)

    def view(self, **kwargs):
        self.minimize()
        view = nv.show_rdkit(self.rdmol, **kwargs)
        return view

    def view_2d(self):
        return Chem.MolFromSmiles(Chem.MolToSmiles(self.rdmol))

    def minimize(self, max_iter: int=1000, nb_threshold: float=25):
        Chem.SymmetrizeSSSR(self.rdmol)
        AllChem.MMFFOptimizeMolecule(self.rdmol, maxIters=max_iter,
                                     nonBondedThresh=nb_threshold)
