
from rdkit import Chem

class Atom:

    def __init__(self, rdatom, name="", charge=0):
        self.index = rdatom.GetIdx()

        if not name:
            name = rdatom.GetSymbol()
        self.name = name
        self.charge = charge
        self.rdatom = rdatom
        info = rdatom.GetPDBResidueInfo()

        #TODO: should I copy the aton and/or info?

        if info is None:
            info = Chem.AtomPDBResidueInfo()
            info.SetName(name)
        
        rdatom.SetMonomerInfo(info)
        self._pdbinfo = info
    
    @property
    def indices(self):
        return self.index
    
    @property
    def resname(self):
        return self._pdbinfo.GetResidueName()

    @resname.setter
    def resname(self, name):
        self._pdbinfo.SetResidueName(name)
    
    @property
    def resid(self):
        return self._pdbinfo.GetResidueNumber()

    @resid.setter
    def resid(self, value):
        self._pdbinfo.SetResidueNumber(value)
    
    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, value):
        self._index = value
        self._pdbinfo.SetSerialNumber(value)
    
    @property
    def altloc(self):
        return self._pdbinfo.GetAltLoc()
        
    @altloc.setter
    def altloc(self, value):
        self._pdbinfo.SetAltLoc(value)

    @property
    def chainid(self):
        return self._pdbinfo.GetChainId()
        
    @chainid.setter
    def chainid(self, value):
        self._pdbinfo.SetChainId(value)

    @property
    def icode(self):
        return self._pdbinfo.GetInsertionCode()
        
    @icode.setter
    def icode(self, value):
        self._pdbinfo.SetInsertionCode(value)

    @property
    def occupancy(self):
        return self._pdbinfo.GetOccupancy()
        
    @occupancy.setter
    def occupancy(self, value):
        self._pdbinfo.SetOccupancy(value)

    @property
    def tempfactor(self):
        return self._pdbinfo.GetTempFactor()
        
    @tempfactor.setter
    def tempfactor(self, value):
        self._pdbinfo.SetTempFactor(value)