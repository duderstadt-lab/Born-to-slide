# Let users know if they're missing any of our hard dependencies
hard_dependencies = ("re", "numpy", "seaborn", "scyjava", "jnius")
missing_dependencies = []

for dependency in hard_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(f"{dependency}: {e}")

if missing_dependencies:
    raise ImportError(
        "Unable to import required dependencies:\n" + "\n".join(missing_dependencies)
    )
del hard_dependencies, dependency, missing_dependencies

__all__ = ['convert']

print(f'{__name__} initialized.')

# module level doc-string
__doc__ = """
marspy - a Python implementation of our beloved MARS Fiji plugin
=====================================================================
v2.0
@author Matthias Scherr

**marspy** is a Python package providing easy use and access to MARS
functionality and data structures.
**Mars** - **M**olecule **AR**chive **S**uite is broad framework for storage and 
reproducible processing of single-molecule observations.

Mars provides a collection of ImageJ2 commands to find, fit, track and characterize 
single molecules. The algorithms provided are routinely used to analyze data arising 
from multicolor single molecule TIRF and DNA flow stretching experiments. Primary 
image analysis commands generate MoleculeArchives which contain collections of 
individual molecules and image metadata records. MoleculeArchives have a simple, 
yet powerful, API for fast and concise access to subsets of records based on 
arbitrary properties. During creation, all records are assigned universally unique 
IDs, which serve as the primary keys for retrieval and storage. This architecture 
allows for seamless virtual storage, merging, and multithreaded processing of very 
large datasets. Moreover, this framework allows for the same data structure to be 
used all the way from initial processing to the generation of final plots either 
in Fiji or in Jupyter Analysis_software using Python.

Main Features
-------------
Here are some awesome features of marspy / MARS:

  - marspy fully grants access to all MARS functionality
  - conversion of .yama files to python

"""
