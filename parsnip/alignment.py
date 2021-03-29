
from . import utils

from rdkit import Chem

class Alignment:

    @classmethod
    def from_tag_names(cls, unit, polymer,
                       monomer_tag_names=[], polymer_tag_names=[],
                       replace_polymer_atoms=True,
                       monomer_tag_index=-1, polymer_tag_index=-1):
        mnames = utils.asiterable(monomer_tag_names)
        pnames = utils.asiterable(polymer_tag_names)
        n_mtags = len(mnames)

        mix = utils.asiterable(monomer_tag_index)
        pix = utils.asiterable(polymer_tag_index)
        if len(pix) != n_mtags and len(pix) == 1:
            pix = pix * n_mtags
        if len(mix) != n_mtags and len(mix) == 1:
            mix = mix * n_mtags
        mtags = [unit.get_tag(x, i)
                 for x, i in zip(mnames, mix)]
        ptags = [polymer.get_tag(x, i)
                 for x, i in zip(pnames, pix)]

        return cls(unit, polymer, monomer_tags=mtags, polymer_tags=ptags,
                   replace_polymer_atoms=replace_polymer_atoms)

    def __init__(self, unit, polymer, monomer_tags=[], polymer_tags=[],
                 replace_polymer_atoms=True):
        
        self.unit = unit
        self.polymer = polymer
        self.monomer_tags = utils.asiterable(monomer_tags)
        self.polymer_tags = utils.asiterable(polymer_tags)
        n_mtags = len(self.monomer_tags)
        err = "Must provide same number of tags for monomers and polymers"
        assert len(self.polymer_tags) == n_mtags, err

        rep = utils.asiterable(replace_polymer_atoms)
        if len(rep) != n_mtags and len(rep) == 1:
            rep = rep * n_mtags
        self.replace_polymer_atoms = rep

    def process(self):
        mp_indices = []
        del_rep_atoms = []
        for mtag, ptag, rep in zip(self.monomer_tags, self.polymer_tags,
                                   self.replace_polymer_atoms):
            mp_indices.extend(list(zip(mtag.indices, ptag.indices)))

            rep = utils.asiterable(rep)
            if len(rep) != len(mtag.atoms) and len(rep) == 1:
                rep = rep * len(mtag.atoms)
            for m, p, r in zip(mtag.atoms, ptag.atoms, rep):
                atoms = [p, m]
                if not r:
                    atoms = atoms[::-1]
                del_rep_atoms.append(atoms)

        Chem.rdMolAlign.AlignMol(self.unit.rdmol, self.polymer.rdmol,
                                 atomMap=mp_indices)
        self.polymer.add_unit_only(self.unit)

        del_rep_atoms = sorted(del_rep_atoms, key=lambda x: x[0].index)
        delete_atoms = [x[0] for x in del_rep_atoms]

        self.polymer.remove_params_within_atoms(delete_atoms)

        for to_delete, to_replace in del_rep_atoms[::-1]:
            self.polymer.update_atom(to_delete, to_replace)
        
        self.polymer.remove_atoms(delete_atoms[::-1])
        self.polymer.clean()