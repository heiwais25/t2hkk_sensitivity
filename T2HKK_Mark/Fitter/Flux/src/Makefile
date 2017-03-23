CXX 		= g++
SRCS	= $(wildcard *.cc)
OBJS	= $(SRCS:.cc=.o)
TARGET	= out

all: $(TARGET)
	$(CXX) -o $(TARGET) $(OBJS) 

$(TARGET):
	$(CXX) -c $(SRCS)

clean:
	rm -f $(TARGET) *.o
