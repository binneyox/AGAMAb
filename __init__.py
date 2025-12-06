from .agamab import *               # import everything from the C++ library
from .agamab import __version__, __doc__  # these two are not automatically imported=
del agamab                          # remove the C++ library from the root namespace
