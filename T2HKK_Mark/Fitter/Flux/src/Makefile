CXX 			= g++
AR 			= ar
LD_SHARED 	= g++

ROOTLIBS   	= -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --libs) -lThread 
ROOTCFLAGS 	= $(shell root-config --cflags)
INCDIR  		= -I../../include 
LIBDIR  		= ../../lib

CXXFLAGS   	= -Wall -mtune=core2 -O3 -g -fPIC -fopenmp $(ROOTCFLAGS) $(INCDIR)
SYSLIB     	= -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS) $(LIBDIR)
SOFLAGS 		= -shared -Wl,-export-dynamic $(ROOTCFLAGS)
LDFLAGS 		=
ARFLAGS 		= rcsv


%.o:%.cc
	g++ -c $(CXXFLAGS) $<

libMakeAtmFlux: libMakeAtmFlux.a libStringParse.a $(LIBDIR)/libMakeAtmFlux.so

libMakeAtmFlux.a: MakeAtmFlux.o
	$(AR) $(ARFLAGS) $@ $^

libStringParse.a: StringParse.o
	$(AR) $(ARFLAGS) $@ $^


$(LIBDIR)/libMakeAtmFlux.so: MakeAtmFlux.o StringParse.o
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS)

	# $(LD_SHARED) -shared -Wl,-soname,$@,-export-dynamic $(ROOTCFLAGS) -o $@ $^


clean:
	rm -fv *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
