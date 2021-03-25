
from pathlib import Path

from .. import WRITERS, READERS
from ..topology import ForceField
from ..lib import utils

class Reader(type):

    def __init__(cls, name, bases, dct):
        type.__init__(cls, name, bases, dct)
        try:
            fmt = utils.asiterable(dct['formats'])
        except KeyError:
            pass
        else:
            for f in fmt:
                # ForceField.READERS[f.upper()] = cls
                READERS[f.upper()] = cls


class BaseReader(metaclass=Reader):

    formats = []

    def __init__(self, filename, forcefield=None):
        if forcefield is None:
            name = Path(filename).name.split('.')[0]
            forcefield = ForceField(name)

        self.filename = filename
        self.forcefield = forcefield


class Writer(type):

    def __init__(cls, name, bases, dct):
        type.__init__(cls, name, bases, dct)
        try:
            fmt = utils.asiterable(dct['formats'])
        except KeyError:
            pass
        else:
            for f in fmt:
                # ForceField.WRITERS[f.upper()] = cls
                WRITERS[f.upper()] = cls


class BaseWriter(metaclass=Writer):

    formats = []

    def __init__(self, filename):
        self.filename = filename



