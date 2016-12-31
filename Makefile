MEXT := $(shell mexext)

BIN := est_divergencemex

all: $(BIN).$(MEXT)


%.$(MEXT) : %.cpp
	mex -silent -output $@ $^


#CC='gcc' CXX='g++' LD='gcc' LDFLAGS='-fopenmp'

.PHONY: clean
clean:
	rm -f $(BIN).$(MEXT)
