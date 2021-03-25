
from rdkit import Chem

def isiterable(obj):
    """
    If an object is iterable and not a string
    """
    if isinstance(obj, str):
        return False
    if hasattr(obj, 'next'):
        return True
    try:
        len(obj)
    except (TypeError, AttributeError):
        return False
    return True

def asiterable(obj):
    if not isiterable(obj):
        obj = [obj]
    return obj

def rdview(mol):
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))