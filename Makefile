# this file is for general settings such as file list, etc.
# machine-specific settings such as include paths and #defines are in Makefile.local
include Makefile.local

# this file is shared between Makefile (Linux/MacOS) and Makefile.msvc (Windows)
# and contains the list of source files and folders
include Makefile.list

LIBNAME_SHARED = Py_agama.so
LIBNAME_STATIC = Py_agama.a
OBJECTS  = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOURCES))
TESTEXE  = $(patsubst %.cpp,$(EXEDIR)/%.exe,$(TESTSRCS))
COMPILE_FLAGS_ALL += -I$(SRCDIR)

# this is the default target (build all), if make is launched without parameters
all:  lib $(TESTEXE)

# one may recompile just the shared and static versions of the library by running 'make lib'
lib:  $(LIBNAME_STATIC) $(LIBNAME_SHARED)

$(LIBNAME_STATIC):  $(OBJECTS) Makefile Makefile.local Makefile.list
	$(AR) ru $(LIBNAME_STATIC) $(OBJECTS)

$(LIBNAME_SHARED):  $(OBJECTS) Makefile Makefile.local Makefile.list
	$(LINK) -shared -o $(LIBNAME_SHARED) $(OBJECTS) $(LINK_FLAGS_ALL) $(LINK_FLAGS_LIB) $(LINK_FLAGS_LIB_AND_EXE_STATIC)

# two possible choices for linking the executable programs:
# 1. shared (default) uses the shared library agama.so, which makes the overall code size smaller,
# but requires that the library is present in the same folder as the executable files.
# 2. static (turned on by declaring an environment variable AGAMA_STATIC) uses the static library libagama.a
# together with all other third-party libraries (libgsl.a, libgslcblas.a, possibly nemo, unsio, glpk, etc.)
# but notably _excluding_ python, since it is used only in two places:
# (a) the python extension module (which is the shared library itself), and
# (b) in math::quadraticOptimizationSolve when compiled with HAVE_CVXOPT,
# but the latter function is not used by any other part of the library or test/example programs,
# except the "solveOpt" function provided by the Python interface (again).
# So it is safe to ignore python in static linking of C++/C/Fortran programs.

ifndef AGAMA_STATIC
# for each executable file, first make sure that the exe/ folder exists,
# and create a symlink named agama.so pointing to ../agama.so in that folder if needed
# (if this was an actual file and not a symlink, then delete it first and then create a symlink)
$(EXEDIR)/%.exe:  $(TESTSDIR)/%.cpp $(LIBNAME_SHARED)
	@mkdir -p $(EXEDIR)
	@[ -f $(EXEDIR)/$(LIBNAME_SHARED) -a ! -L $(EXEDIR)/$(LIBNAME_SHARED) ] && rm $(EXEDIR)/$(LIBNAME_SHARED) || true
	@[ -L $(EXEDIR)/$(LIBNAME_SHARED) ] || ln -s ../$(LIBNAME_SHARED) $(EXEDIR)/$(LIBNAME_SHARED)
	$(LINK) -o "$@" "$<" $(COMPILE_FLAGS_ALL) $(LIBNAME_SHARED) $(LINK_FLAGS_ALL) $(LINK_FLAGS_EXE_SHARED)
else
$(EXEDIR)/%.exe:  $(TESTSDIR)/%.cpp $(LIBNAME_STATIC)
	@mkdir -p $(EXEDIR)
	$(LINK) -o "$@" "$<" $(COMPILE_FLAGS_ALL) $(LIBNAME_STATIC) $(LINK_FLAGS_ALL) $(LINK_FLAGS_LIB_AND_EXE_STATIC)
endif


$(OBJDIR)/%.o:  $(SRCDIR)/%.cpp Makefile.local
	@mkdir -p $(OBJDIR)
	$(CXX) -c $(COMPILE_FLAGS_ALL) $(COMPILE_FLAGS_LIB) -o "$@" "$<"

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(EXEDIR)/*.exe $(LIBNAME_SHARED) $(EXEDIR)/$(LIBNAME_SHARED) $(LIBNAME_STATIC)

nemo:
amuse:

# auto-dependency tracker (works with GCC-compatible compilers?)
DEPENDS = $(patsubst %.cpp,$(OBJDIR)/%.d,$(SOURCES))
COMPILE_FLAGS_LIB += -MMD -MP
-include $(DEPENDS)

.PHONY: clean test lib doxy nemo amuse