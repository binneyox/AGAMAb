import os
os.aff_dll_directory("h:/u/c/agama/agama/x64/release")
from Py_agama import *            # import everything from the C++ library
from Py_agama import __version__, __doc__  # these two are not automatically imported=
del Py_agama                          # remove the C++ library from the root namespace
