# Makefile to produce a super-simple monte carlo example
# Need to control the location of the header and library
# for HrsTrkCorr, specifically $DIR_TRKCORR
# You'll also need a link to the directory "holefiles".
# e.g. ln -s ../holefiles holefiles
#
# R. Michaels, June 2021


# Choose the compiler.
GCC=g++
GLD=g++


ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lvdt
INCLUDES      = -I$(ROOTSYS)/include
CXX           = $(GCC)
CXXFLAGS      = -fno-exceptions -fpermissive -std=c++11 -fPIC $(INCLUDES)
LD            = $(GLD)
SOFLAGS       = -std=c++11  -shared 
GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 
LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) /usr/lib64/libg2c.so.0

MAKEDEPEND    = $(GCC)

# These needs to point to wherever the header HrsTrkCorr.h
# and library libtrkcorr.so  is located
DIR_TRKCORR = /home/robertmichaels/TrkCorr/HrsTrkCorr
HDR_TRKCORR = -I$(DIR_TRKCORR)
LIB_TRKCORR = $(DIR_TRKCORR)/libtrkcorr.so

INCLUDES += $(HDR_TRKCORR)
ALL_LIBS = $(LIBS) $(LIB_TRKCORR)

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g -ggdb
endif

SRC = trktest_main.C

# An executible test and the library.
# note, the library can also be used in a ROOT script like trk_script.C

all:  xtest

xtest: $(OBJS) $(HEAD) trktest_main.C
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ trktest_main.C $(OBJS) $(ALL_LIBS)

clean:
	rm -f *.o core xtest

realclean:  clean
	rm -f *.tar  *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

