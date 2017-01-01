MEXT := $(shell mexext)

BIN := est_divergencemex

all: $(BIN).$(MEXT)


%.$(MEXT) : %.o
	mex -silent -output $@ $^ LD='gcc-6' CC='gcc-6' LDFLAGS='-fopenmp \$$LDFLAGS'

%.o : %.cpp
	g++-6 -fopenmp -o $@ -c $^ -I/usr/local/matlab/extern/include

.PHONY: clean
clean:
	rm -f $(BIN).$(MEXT)
