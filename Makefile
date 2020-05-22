# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
DFANAME = exec/DFA
DFADEPS = src/detrendFluctuationAnalysis.c
HENAME = exec/HE
HEDEPS = src/HurstExponent.c
SSNAME = sunSpots/exec/dataAnalysis.sh

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(DFANAME) $(DFADEPS)
	$(CC) $(CFLAGS) -o $(HENAME) $(HEDEPS)


.PHONY: sunSpots

sunSpots:
	bash $(SSNAME)

.PHONY: clean

clean: 
	$(RM) $(DFANAME)
	$(RM) $(HENAME)
