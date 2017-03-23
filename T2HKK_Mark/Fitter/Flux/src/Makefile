CXX 		= g++
SRCS	= $(wildcard *.cc)
OBJS	= $(SRCS:.cc=.o)
TARGET	= out

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) 

$(OBJS): $(SRCS)
	$(CXX) -c $(SRCS)

clean:
	rm -f $(TARGET) *.o
