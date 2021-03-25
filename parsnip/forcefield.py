import inspect

from ..get_formats import _get_reader
from .topology import TopologyTemplate

def _get_format(filename, format=None):
    if inspect.isclass(format):
        return format
    if format is None and isinstance(filename, str):
        return filename.split('.')[-1].upper()

    raise ValueError(f"Cannot find format for {filename}")


class ForceField:

    registry = {}

    def __init__(self, name=None):
        if name is not None:
            self.registry[name] = self
        self.name = name
        self.atomtypes = {}
        self.moltypes = {}
        self.paramtypes = {}
        self.paramtypes_by_name = {}

    def update(self, filename, format=None):
        reader = _get_reader(filename, format=format)
        reader(filename, forcefield=self).parse()

    def _add_mol_topology(self, name, nrexcl=3):
        self.moltypes[name] = TopologyTemplate(name, nrexcl=nrexcl)
        return self.moltypes[name]

