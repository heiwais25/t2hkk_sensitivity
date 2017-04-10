SRCS = $(wildcard *.cc)
TRGS = $(patsubst %.cc, ../bin/%, $(SRCS))
DICT = $(patsubst %.cc, %Dict.C, $(SRCS))

all : $(TRGS)

ROOTLIBS := $(shell root-config --glibs --libs)  -lMinuit -lMinuit2
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags) -I../include

NUPRISMLIBS := -L../lib -L$(TRACKERFILE_LIB) -L$(OAANALYSISLIBS) -L../../Prob3++ -lOscParams -lSystParams -lMakeHistograms -lMakeAtmFlux -lThreeProb_2.10 -lThrowParams

$(TRGS) : $(SRCS)
$(DICT) : $(SRCS)

../bin/% : %.cc %Dict.C 
	g++ -mtune=core2 -O3 -fopenmp -g $(ROOTCFLAGS) -I $(ROOTINCLUDES) -o $@ $(OBJECTFILES) $< $*Dict.C $(NUPRISMLIBS) $(ROOTLIBS) -Wl,-rpath,../lib 

%Dict.C :
	rootcint -f $@ -c Linkdef.h

clean:
	rm -fv ../bin/*
	rm -fv *Dict.C
	rm -fv *Dict.h
