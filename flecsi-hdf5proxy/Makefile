BIN := hdf5proxy

default: $(BIN)

H5C++ ?= h5c++

CXXFLAGS += -std=c++11 -g -O0

$(BIN): %:%.cc
		$(H5C++) -o $@ $(CXXFLAGS) $< $(LDFLAGS)

clean:
		rm -f *.o $(BIN)

