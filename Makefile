MEXT := $(shell mexext)

BIN := est_divergencemex

all: $(BIN).$(MEXT)


%.$(MEXT) : %.cpp
	mex -silent -output $@ $^


.PHONY: clean
clean:
	rm -f $(BIN).$(MEXT) *.o
