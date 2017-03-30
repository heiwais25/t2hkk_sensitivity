CXX=g++

SRCS 	= $(wildcard *.cc)
OBJS	= $(SRCS:.cc=.o)
TARGET	= out


# NUPRISMLIBS := -L../../lib -L$(TRACKERFILE_LIB) -L$(OAANALYSISLIBS) -L../../../Prob3++ -lOscParams -lSystParams -lMakeHistograms -lThreeProb_2.10 -lThrowParams
NUPRISMLIBS := -L../../lib -L../../../Prob3++ -lOscParams -lSystParams -lMakeHistograms -lThreeProb_2.10 -lThrowParams
ROOTLIBS := $(shell root-config --glibs --libs)  -lMinuit -lMinuit2
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags) -I../../include
SOFLAGS= -shared -Wl,-export-dynamic $(ROOTCFLAGS)
ROOTLDFLAGS := $(shell root-config --ldflags)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(ROOTLDFLAGS) -o $(TARGET) $(OBJS) $(ROOTLIBS) $(NUPRISMLIBS)

$(OBJS): $(SRCS)
	$(CXX) -mtune=core2 -O3 -fopenmp -c $(ROOTCFLAGS) -I $(ROOTINCLUDES) $(SRCS)

clean:
	rm -rf $(TARGET) $(OBJS)