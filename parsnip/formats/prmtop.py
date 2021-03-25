from MDAnalysis.topology.TOPParser import TOPParser, Z2SYMB
from MDAnalysis.lib.util import openany

from .. import TOPOLOGY_TYPES
from .base import BaseReader, BaseWriter

# TODO:

# %FLAG NONBONDED_PARM_INDEX
# %FLAG RESIDUE_LABEL
# %FLAG RESIDUE_POINTER
# %FLAG BOND_FORCE_CONSTANT
# %FLAG BOND_EQUIL_VALUE
# %FLAG ANGLE_FORCE_CONSTANT
# %FLAG ANGLE_EQUIL_VALUE
# %FLAG DIHEDRAL_FORCE_CONSTANT
# %FLAG DIHEDRAL_PERIODICITY
# %FLAG DIHEDRAL_PHASE
# %FLAG SCEE_SCALE_FACTOR
# %FLAG SCNB_SCALE_FACTOR
# %FLAG SOLTY
# %FLAG LENNARD_JONES_ACOEF
# %FLAG LENNARD_JONES_BCOEF
# %FLAG BONDS_INC_HYDROGEN
# %FLAG BONDS_WITHOUT_HYDROGEN
# %FLAG ANGLES_INC_HYDROGEN
# %FLAG ANGLES_WITHOUT_HYDROGEN
# %FLAG DIHEDRALS_INC_HYDROGEN
# %FLAG DIHEDRALS_WITHOUT_HYDROGEN
# %FLAG EXCLUDED_ATOMS_LIST
# %FLAG HBOND_ACOEF
# %FLAG HBOND_BCOEF
# %FLAG HBCUT
# %FLAG AMBER_ATOM_TYPE
# %FLAG TREE_CHAIN_CLASSIFICATION
# %FLAG JOIN_ARRAY
# %FLAG IROTAT
# %FLAG RADIUS_SET
# %FLAG RADII
# %FLAG SCREEN
# %FLAG IPOL

class PrmtopReader(TOPParser, BaseReader):
    formats = ["PRMTOP"]

    def __init__(self, *args, **kwargs):
        super(PrmtopReader, self).__init__(*args, **kwargs)

        self.atom_sections = {
            "ATOM_NAME": ("atomnames", 20, 0, lambda x: x),
            "CHARGE": ("charges", 5, 0, lambda x: float(x)/18.2223),
            "ATOMIC_NUMBER": ("elements", 10, 0, lambda x: int(x)),
            "MASS": ("masses", 5, 0, lambda x: float(x)),
            "ATOM_TYPE_INDEX": ("type_indices", 10, 0, lambda x: int(x)),
            "NUMBER_EXCLUDED_ATOMS": ("n_exc", 10, 0, lambda x: int(x)),


        }
        
        self.parsers = {
            "NONBONDED_PARM_INDEX": ("nb_indices", self.parse_nb_parm_index),

        }

    def parse_nb_parm_index(self):
        num = self.sys_info[1] ** 2
        per_line = 10
        numlines = (num // per_line)
        if num % per_line != 0:
            numlines += 1
        
        return self.parsesection_mapper(numlines, lambda x: int(x))

    def parse(self):
        with openany(self.filename) as self.topfile:
            header = next(self.topfile)

        if not header.startswith("%VE"):
                raise ValueError(
                    "{0} is not a valid TOP file. %VE Missing in header"
                    "".format(self.filename))
            title = next(self.topfile).split()
            # %FLAG TITLE
            if not (title[1] == "TITLE"):
                # Raise a separate warning if Chamber-style TOP is detected
                if title[1] == "CTITLE":
                    emsg = ("{0} is detected as a Chamber-style TOP file. "
                            "At this time PolyTop does not support such "
                            "topologies".format(self.filename))
                else:
                    emsg = ("{0} is not a valid TOP file. "
                            "'TITLE' missing in header".format(self.filename))
                raise ValueError(emsg)

            # %FLAG POINTERS
            while not header.startswith('%FLAG POINTERS'):
                header = next(self.topfile)
            next(self.topfile)

            topremarks = [next(self.topfile).strip() for i in range(4)]
            sys_info = [int(k) for i in topremarks for k in i.split()]

            header = next(self.topfile)
            # grab the next section title
            next_section = header.split("%FLAG")[1].strip()
            ...

    
    def parse_atomnames(self, num_per_record, numlines):
        # %FLAG ATOM_NAME
        return self.parsesection_mapper(numlines, lambda x: x)

    def parse_resnames(self, num_per_record, numlines):
        return self.parsesection_mapper(numlines, lambda x: x)

    def parse_charges(self, num_per_record, numlines):
        # %FLAG CHARGE
        vals = self.parsesection_mapper(numlines, lambda x: float(x))
        charges = np.array(vals, dtype=np.float32)
        charges /= 18.2223  # to electron charge units
        return charges

    def parse_elements(self, num_per_record, numlines):
        # %FLAG ATOMIC_NUMBER
        return self.parsesection_mapper(numlines, lambda x: int(x))

    def parse_masses(self, num_per_record, numlines):
        # %FLAG MASS
        return self.parsesection_mapper(numlines, lambda x: float(x))

    def parse_type_indices(self, num_per_record, numlines):
        """Extracts the index of atom types of the each atom involved in Lennard
        Jones (6-12) interactions.

        Parameters
        ----------
        num_per_record : int
            The number of entries for each record in section (unused input)
        numlines : int
            The number of lines to be parsed in current section
        """
        # %FLAG ATOM_TYPE_INDEX
        return self.parsesection_mapper(numlines, lambda x: int(x))

