
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topology as mda_top
from MDAnalysis.core import topologyattrs as mda_ta
from MDAnalysis.topology.base import change_squash

from ..lib.utils import cached, uncache
from .params import ParamContainer
from .paramlist import AtomList
from .topologyobjects import Residue
from ..get_formats import _get_writer


class TopologyBase:

    def __repr__(self):
        pms = [f"{len(self.params[x])} {x}" for x in self.params]
        rep = (f"<{type(self).__name__} with {len(self.atoms)} atoms, "
               f"{len(self.residues)} residues, {', '.join(pms)}>")
        return rep

    def __sub__(self, other):
        new = type(self)(name=self.name, nrexcl=self.nrexcl)
        res = [r for r in self.residues if r not in other.residues]
        new.residues.extend(res)
        new.params._update_nocheck(self.params - other.params)
        return new

    def __bool__(self):
        return bool(len(self.atoms) + self.params.n_params)

    def __add__(self, other):
        if not other:
            return self.copy()
        new = type(self)(name=self.name, nrexcl=self.nrexcl)
        new.residues.extend(self.residues)
        new.params._update_nocheck(self.params)
        new.extend(other)
        return new

    def __iadd__(self, other):
        self.extend(other)
        return self

    def __radd__(self, other):
        if not other:
            return self.copy()
        return other.__add__(self)

    def __getattr__(self, attr):
        if attr in self.params.params:
            return self.params.params[attr]
        raise AttributeError(f"Could not find attribute {attr}")

    def copy(self):
        new = type(self)(name=self.name, nrexcl=self.nrexcl)
        for r in self.residues:
            new.residues.append(r.copy(new))
        new.params.update(self.params)
        return new

    def deepcopy(self):
        new = type(self)(name=self.name, nrexcl=self.nrexcl)
        for r in self.residues:
            new.residues.append(r.deepcopy(new))
        newparams = self.params.deepcopy()
        new.params._update_nocheck(newparams)
        return new

    def write(self, filename, format=None):
        for v in self.params.values():
            v.sort()
        writer = _get_writer(filename, format=format)
        writer(filename).write(self)

    def _reindex(self):
        past_cg = -1
        next_cg = 0
        for i, a in enumerate(self.atoms):
            a._ix = i
            a._serial = i+1
            if a.charge_group != past_cg:
                next_cg += 1
            past_cg = a.charge_group
            a.charge_group = next_cg
        
        for i, r in enumerate(self.residues, 1):
            r.resid = i
            r.resname = r.resname
        self.params._reindex()

    def extend(self, other):
        self.residues.extend(other.residues)
        self.params.update(other.params)

    def _extend_nocheck(self, other):
        self.residues.extend(other.residues)
        self.params._update_nocheck(other.params)

    @uncache("atoms", "_serial_mapping")
    def add_residue(self, name="UNK", resid=0):
        self.residues.append(Residue(self, name=name, resid=resid))
        return self.residues[-1]

    @uncache("atoms", "_serial_mapping")
    def add_atom(self, residue=None, **kwargs):
        if residue is None:
            if len(self.residues):
                residue = self.residues[-1]
            else:
                residue = self.add_residue()
        residue.add_atom(**kwargs)
        return residue.atoms[-1]

    def reindex_atoms(self, indices):
        index_mapping = {}
        for a, i in zip(self.atoms, indices):
            index_mapping[a._ix] = i
            a._ix = i
        self.params.remap_indices(index_mapping)

    def _reindex_params(self):
        """This relies on the presumption that the indices in the params
        follows the order of the atoms in the topology"""
        # index_mapping = dict(zip(self._atom_indices, self.indices))
        self.params._reindex()

    @property
    # @cached
    def atoms(self):
        return AtomList([a for r in self.residues for a in r.atoms])

    def _create_universe(self, coordinates=None):
        # self._reindex()
        resids = self.atoms.resids
        resnames = self.atoms.resnames

        serials = self.atoms.serials

        residx, (resids, resnames) = change_squash((resids, resnames),
                                                   (resids, resnames))

        n_atoms = len(self.atoms)
        n_residues = len(resids)
        # atom_types = self.atoms.atom_types

        # attrs
        attrs = [
            mda_ta.Atomnames(self.atoms.names),
            mda_ta.Atomids(serials),
            mda_ta.Charges(self.atoms.charges),
            # mda_ta.Atomtypes(a_types),
            # mda_ta.Elements(elements),
            mda_ta.Masses(self.atoms.masses),
            mda_ta.Resids(resids),
            mda_ta.Resnums(resids.copy()),
            mda_ta.Resnames(resnames)
        ]

        topattr_params = self.params.to_mda_topologyattrs()
        # attrs.extend(topattr_params)

        elements = self.atoms.elements
        if all(elements):
            attrs.append(mda_ta.Elements(elements))

        n_atoms = len(self.atoms)
        n_residues = len(resids)
        
        top = mda_top.Topology(n_atoms=n_atoms, n_res=n_residues, n_seg=1,
                               atom_resindex=residx, residue_segindex=None,
                               attrs=attrs)
        u = mda.Universe(top)
        for tp in topattr_params:
            u._topology.add_TopologyAttr(tp)
            u._process_attr(tp)

        # ye gods PDB warnings
        for attr in ("altLocs", "icodes", "segids", "occupancies", "tempfactors"):
            u.add_TopologyAttr(attr)
        
        if coordinates is not None:
            u.load_new(coordinates)
            u.dimensions = [0, 0, 0, 90, 90, 90]
        return u

    def all_in_indices(self, indices, name=None, index_mapping={}):
        if name is None:
            name = self.name
        new = type(self)(name=name)

        atoms = self.atoms[indices]
        indices = atoms.indices

        res = [r.within_atoms(atoms, topology=new) for r in self.residues]

        new.residues.extend([r for r in res if len(r)])
        params = self.params.all_in_indices(indices)
        if index_mapping:
            self.atoms.remap_indices(index_mapping)
            params.remap_indices(index_mapping)
        new.params.update(params)
        return new

    def any_in_indices(self, indices, name=None, index_mapping={}):
        if name is None:
            name = self.name
        new = type(self)(name=name)

        atoms = self.atoms[indices]
        indices = atoms.indices

        res = [r.within_atoms(atoms, topology=new) for r in self.residues]
        new.residues.extend([r for r in res if len(r)])
        params = self.params.any_in_indices(indices)
        if index_mapping:
            self.atoms.remap_indices(index_mapping)
            params.remap_indices(index_mapping)
        new.params.update(params)
        return new

    def xor_in_indices(self, indices, name=None, index_mapping={}):

        new = type(self)()

        if not len(indices):
            return type(self)()
        indices = [a.index for a in self.atoms[indices]]
        params = self.params.xor_in_indices(indices)
        if index_mapping:
            self.atoms.remap_indices(index_mapping)
            params.remap_indices(index_mapping)
        new.params._update_nocheck(params)
        return new


class TopologyTemplate(TopologyBase):

    def __init__(self, name=None, nrexcl=3):
        self._cache = {}
        self.name = name
        self.nrexcl = nrexcl
        self.residues = []
        self.external_atoms = []
        self.params = ParamContainer(self)
        self._topology = self

    @classmethod
    def from_topology(cls, topology):
        new = cls(name=topology.name, nrexcl=topology.nrexcl)
        for r in topology.residues:
            new.residues.append(r.deepcopy())
        new.params.update(topology.params.deepcopy())
        new._reindex()
        return new


    def subset_monomer_by_index(self, monomer, remap_indices=True):
        if len(monomer.atoms) > len(self.atoms):
            raise ValueError("Missing topology atoms")
        self._reindex()
        atoms = self.atoms[monomer.atoms.indices]
        indices = [a.index for a in atoms]

        aix = self.all_in_indices(indices, name=monomer.name)
        # print("aix", aix)
        new = Topology.from_template(aix)
        new._reindex()

        return new

    def subset_monomer_by_name(self, monomer, remap_indices=True):
        name_mapping = {a.name: a for a in self.atoms}
        atoms = []
        for a in monomer.atoms:
            try:
                atoms.append(name_mapping[a.name])
            except KeyError:
                err = f"Monomer atom {a.name} not found in topology"
                raise ValueError(err) from None
        self._reindex()
        indices = [a.index for a in atoms]
        index_mapping = {x: i for i, x in enumerate(indices)}
        aix = self.all_in_indices(indices, name=monomer.name, index_mapping=index_mapping)
        # print("aix", aix)
        new = Topology.from_template(aix)
        new._reindex()
        return new 
    
    def get_indices_from_serial(self, serials):
        return tuple(self._serial_mapping[s] for s in serials)

    @property
    @cached
    def _serial_mapping(self):
        return {a.serial: i for i, a in enumerate(self.atoms)}
        


class Topology(TopologyBase, type):

    def __repr__(self):
        pms = [f"{len(self.params[x])} {x}" for x in self.params]
        rep = (f"<{type(self).__name__} '{self.name}' with {len(self.atoms)} atoms, "
               f"{len(self.residues)} residues, {', '.join(pms)}>")
        return rep

    def __getattr__(self, attr):
        if attr in self.params.params:
            return self.params.params[attr]
        raise AttributeError(f"Could not find attribute {attr}")

    def __getattribute__(self, attr):
        # well this is awful
        dct = type.__getattribute__(self, "__dict__")
        if attr in dct:
            return dct[attr]
        return object.__getattribute__(self, attr)

    def __iadd__(self, other):
        self.extend(other)
        return self

    def extend(self, other):
        self.residues.extend(other.residues)
        self.params.update(other.params)

    def get_residues_within_indices(self, indices):
        res = [r.within_indices(indices) for r in self.residues]
        return [r for r in res if len(r)]

    @classmethod
    def from_template(cls, template):
        new = cls(name=template.name, nrexcl=template.nrexcl)
        for r in template.residues:
            new.residues.append(r.deepcopy())
        new.params.update(template.params.deepcopy())
        new._reindex()
        new._template = template
        return new

    # def deepcopy(self):
    #     new = type(self)(name=self.name, nrexcl=self.nrexcl)
    #     for r in self.residues:
    #         new.residues.append(r.deepcopy())
    #     new.params._update_nocheck(self.params.deepcopy())
    #     return new

    def __new__(cls, **kwargs):
        return super().__new__(cls, cls.__name__ + 'Unit', (TopologyUnit,), {})
    
    def __init__(self, name=None, nrexcl=3):
        self._cache = {}
        self.name = name
        self.nrexcl = nrexcl
        self.residues = []
        self.external_atoms = []
        self.params = ParamContainer(self)
        self._template = self
        self._topology = self

    @property
    def template(self):
        if self._template is not self:
            return self._template

        return TopologyTemplate.from_topology(self)

    # @property
    # def _topology(self):
    #     return self
    
    @property
    def indices(self):
        return np.array([a.ix for a in self.atoms])

    def delete_atoms(self, indices):
        res_atom = {}
        for res in self.residues:
            for atom in res.atoms:
                res_atom[len(res_atom)] = (res, atom)
        n_atoms = len(res_atom)
        to_keep = [i for i in np.arange(n_atoms) if i not in indices]
        index_mapping = dict(zip(to_keep, np.arange(n_atoms)))
        for i in indices:
            res, atom = res_atom[i]
            res.atoms.remove(atom)
        self.params.remap_or_delete(index_mapping)
        self._reindex(params=False)


    @property
    # @cached
    def atoms(self):
        return AtomList([a for r in self.residues for a in r.atoms])

    
    def _reindex(self, params=True):
        past_cg = -1
        next_cg = 0
        for i, a in enumerate(self.atoms):
            a._ix = i
            a._serial = i+1
            if a.charge_group != past_cg:
                next_cg += 1
            past_cg = a.charge_group
            a.charge_group = next_cg
        
        for i, r in enumerate(self.residues, 1):
            r.resid = i
            r.resname = r.resname
        
        if params:
            self.params._reindex()
    
    

class TopologyUnit(TopologyBase):

    @classmethod
    def xor_in_indices(cls, indices):
        new = cls()
        new.params._update_nocheck(cls.params.xor_in_indices(indices))
        return new

    @classmethod
    def all_in_indices(cls, indices):
        new = cls()
        new.residues = cls.get_residues_within_indices(indices)
        new.params._update_nocheck(cls.params.all_in_indices(indices))
        return new

    @classmethod
    def any_in_indices(cls, indices):
        new = cls()
        new.residues = cls.get_residues_within_indices(indices)
        new.params._update_nocheck(cls.params.any_in_indices(indices))
        return new

    def __init__(self, *args, **kwargs):
        self.residues = []
        self.params = ParamContainer(self)
