# Makefile to produce library libtrkcorr.so 
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

LIB_TRKCORR = /home/robertmichaels/trkcorr/libtrkcorr.so

ALL_LIBS = $(LIBS) 

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g -ggdb
endif

SRC = HrsTrkCorr.C

HEAD = $(SRC:.C=.h) 
OBJS = $(SRC:.C=.o)

# An executible test and the library.
# note, the library can also be used in a ROOT script like trk_script.C

all:  libtrkcorr.so 

libtrkcorr.so: $(OBJS) $(HEAD)
	$(CXX) $(SOFLAGS) -O -o libtrkcorr.so $(OBJS) $(LIBS)

clean:
	rm -f *.o core libtrkcorr.so 

realclean:  clean
	rm -f *.tar  *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

