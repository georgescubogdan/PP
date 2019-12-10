# Choose the head directory of Speculoos
speculoos_root = ../
# Choose the name of the user program to compile and link with Speculoos
#main           = baffledChannel.cxx
#main           = main2D_sqcavity_ns.cxx
main           = main3D_poisson.cxx

scons     = $(speculoos_root)/scons/scons.py -f $(speculoos_root)/scons/SpecCons


# Choose compilation flags:
# ========================
#
# optimize:    Decide whether to use optimization flags or not
# debug:       Decide whether to use debug flags or not
# parallel:    Decide whether to compile for parallel execution with MPI or not
# platform:    This parameter influences the file param.hxx, for the choice of name
#              mangling scheme during inclusion of BLAS/LAPACK. Choose between
#              GNU, SILICON_GRAPHICS, HEWLETT_PACKARD, DIGITAL, SUN, CRAY_T3D, BLUE_GENE
# serialCXX:   Name of the C++ compiler for serial execution
# parallelCXX: Name of the C++ compiler for parallel execution
# optimFlags:  Optimization flags to use in the case optimize=true
# debugFlags:  Debug flags to use in the case debug=true
# includePath: List of pathnames in which to search for includefiles, separated by white-spaces
# libpath:     List of pathnames in which to search for libraries, separated by white-spaces
# libs:        List of libraries to link with the executabel, separated by white-spaces

compile:
	python $(scons) \
            speculoos_root=$(speculoos_root) \
            target=$(main)                   \
            optimize=true                    \
            debug=true                     \
            parallel=true                   \
            platform=GNU                     \
            serialCXX=g++                    \
            parallelCXX=mpicxx               \
            optimFlags=-O3                   \
            debugFlags=-g                    \
            includepath=$(INCLUDE_PATH)      \
            libpath='../lib'                 \
            libs="lapack blas gfortran"

clean:
	python $(scons) -c \
            speculoos_root=$(speculoos_root) \
            target=$(main)                   \
            parallel=true
	\rm -vf `find $(speculoos_root) -name '*~'`
