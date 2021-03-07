"""
parsnip
Build molecular topologies.
"""

# Add imports here
from .parsnip import *

from .monomer import Monomer

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
