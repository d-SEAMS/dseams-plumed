USE=core tools

# seams-core is not part of PLUMED, so this module stays off unless the
# caller passes its headers and libraries, for example
#   make SEAMS_CXXFLAGS="-I/path/to/seams-core/include" \
#        SEAMS_LIBS="-L/path/to/seams-core/lib -lyodaLib"
# The loadable module used by LOAD is built with meson from this directory.
CXXFLAGS:=$(CXXFLAGS) $(SEAMS_CXXFLAGS)
DYNAMIC_LIBS:=$(DYNAMIC_LIBS) $(SEAMS_LIBS)

include ../maketools/make.module
