CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS   = -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --libs) -lThread 
ROOTCFLAGS = $(shell root-config --cflags)
INCDIR = -I../../include
LIBDIR = ../../lib


#CLHEPINC = $(CLHEP_BASE_DIR)/include
#CLHEPLIB = -L$(CLHEP_BASE_DIR)/lib -l$(CLHEP_LIB)

CXXFLAGS   = -Wall -O3 -g -fPIC $(ROOTCFLAGS) $(INCDIR)
SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS) $(LIBDIR)
SOFLAGS= -shared -Wl,-export-dynamic $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv


%.o:%.cc
	g++ -c $(CXXFLAGS) $<

libThrowParams: libThrowParams.a $(LIBDIR)/libThrowParams.so

libThrowParams.a: ThrowParams.o
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libThrowParams.so: libThrowParams.a ThrowParams.o
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS)

clean:
	rm -fv *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
