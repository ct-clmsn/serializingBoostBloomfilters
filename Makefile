SRC_DIR=$(HOME)/git
BOOST_BLOOMFILTER=$(SRC_DIR)/boost-bloom-filters

CC=g++
CC_FLAGS=-std=c++11 -O2
LIBS=-lboost_serialization

all:
	$(CC) $(CC_FLAGS) -I$(BOOST_BLOOMFILTER) demo.cpp $(LIBS)

clean:
	rm a.out 
