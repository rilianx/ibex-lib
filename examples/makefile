SRCS=$(wildcard *.cpp)
BINS=$(SRCS:.cpp=)

CXXFLAGS := $(shell pkg-config --cflags ibex) 
LIBS	 := $(shell pkg-config --libs  ibex)
LIBDIR	 := $(shell pkg-config --libdir  ibex) 

ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -m32 -O0 -g -pg -Wall -frounding-math  
else
CXXFLAGS := $(CXXFLAGS) -m32 -O3 -DNDEBUG -Wno-deprecated -frounding-math 
endif

all: $(BINS)

% :	%.cpp
	echo $(BINS) $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(BINS)
	
