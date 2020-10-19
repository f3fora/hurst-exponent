# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
DFANAME = exec/DFA
DFADEPS = src/detrendFluctuationAnalysis.c
HENAME = exec/HE
HEDEPS = src/HurstExponent.c
DSSNAME = sunSpots/exec/dailyAnalysis.sh
MSSNAME = sunSpots/exec/monthlyAnalysis.sh
TSSNAME = sunSpots/exec/trueAnalysis.sh

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(DFANAME) $(DFADEPS)
	$(CC) $(CFLAGS) -o $(HENAME) $(HEDEPS)

.PHONY: dailySunSpots

dailySunSpots:
	bash $(DSSNAME)
	
.PHONY: monthlySunSpots

monthlySunSpots:
	bash $(MSSNAME)

.PHONY: trueSunSpots

trueSunSpots:
	bash $(TSSNAME)

.PHONY: clean

clean: 
	$(RM) $(DFANAME)
	$(RM) $(HENAME)
