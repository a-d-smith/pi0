#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = pi0_ERToolsPi0RecoDict
SHLIB = libpi0_ERToolsPi0Reco.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell larlite-config --includes)
INCFLAGS += $(shell seltool-config --includes)

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
LDFLAGS += $(shell seltool-config --libs)
LDFLAGS += $(shell larlite-config --libs)
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE



