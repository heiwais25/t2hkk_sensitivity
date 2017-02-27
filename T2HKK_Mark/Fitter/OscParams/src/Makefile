CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS   = -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --libs) -lThread 
ROOTCFLAGS = $(shell root-config --cflags)
INCDIR = -I../../include 
LIBDIR = ../../lib


#CLHEPINC = $(CLHEP_BASE_DIR)/include
#CLHEPLIB = -L$(CLHEP_BASE_DIR)/lib -l$(CLHEP_LIB)

CXXFLAGS   = -Wall -mtune=core2 -O3 -g -fPIC -fopenmp $(ROOTCFLAGS) $(INCDIR)
SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS) $(LIBDIR)
SOFLAGS= -shared -Wl,-export-dynamic $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv


%.o:%.cc
	g++ -c $(CXXFLAGS) $<

libOscParams: libOscParams.a $(LIBDIR)/libOscParams.so

libOscParams.a: OscParams.o
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libOscParams.so: libOscParams.a OscParams.o
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS)

clean:
	rm -fv *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
