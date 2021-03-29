
from rdkit import Chem
from . import rdutils

class Atom:

    _default_resname = "UNK"
    _default_resid = 1

    def __init__(self, rdatom, name=None, charge=0):
        info = rdutils.get_pdbinfo(rdatom)
        self._pdbinfo = info
        
        self.rdatom = rdatom
        self.rdatom.SetMonomerInfo(info)
        
        self.index = rdatom.GetIdx()
        if not name:
            name = rdatom.GetSymbol()
        self.name = name
        self.charge = charge
        

    @property
    def rdatom(self):
        return self._rdatom

    @rdatom.setter
    def rdatom(self, atom):
        self._rdatom = atom

    @property
    def pdbinfo(self):
        return self._pdbinfo

    @pdbinfo.setter
    def pdbinfo(self, info):
        self._pdbinfo = info
    
    @property
    def indices(self):
        return self.index
    
    @property
    def resname(self):
        try:
            return self._pdbinfo.GetResidueName()
        except AttributeError:
            return self._default_resname

    @resname.setter
    def resname(self, name):
        self._pdbinfo.SetResidueName(name)
    
    @property
    def resid(self):
        try:
            return self._pdbinfo.GetResidueNumber()
        except AttributeError:
            return self._default_resid

    @resid.setter
    def resid(self, value):
        self._pdbinfo.SetResidueNumber(int(value))
    
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