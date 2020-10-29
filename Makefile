# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
DFANAME = exec/DFA
DFADEPS = src/detrendFluctuationAnalysis.c
HENAME = exec/HE
HEDEPS = src/HurstExponent.c
FONAME = exec/FO
FODEPS = src/fourierDFA.c
PRNAME = exec/PR
PRDEPS = src/profile.c
DSSNAME = sunSpots/exec/dailyAnalysis.sh
MSSNAME = sunSpots/exec/monthlyAnalysis.sh
TSSNAME = sunSpots/exec/trueAnalysis.sh
FSSNAME = sunSpots/exec/fourierAnalysis.sh

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(DFANAME) $(DFADEPS)
	$(CC) $(CFLAGS) -o $(HENAME) $(HEDEPS)
	$(CC) $(CFLAGS) -o $(FONAME) $(FODEPS)
	$(CC) $(CFLAGS) -o $(PRNAME) $(PRDEPS)

.PHONY: dailySunSpots

dailySunSpots:
	bash $(DSSNAME)
	
.PHONY: monthlySunSpots

monthlySunSpots:
	bash $(MSSNAME)

.PHONY: trueSunSpots

trueSunSpots:
	bash $(TSSNAME)

.PHONY: fourierSunSpots

fourierSunSpots:
	bash $(FSSNAME)

.PHONY: clean

clean: 
	$(RM) $(DFANAME)
	$(RM) $(HENAME)
	$(RM) $(FONAME)
	$(RM) $(PRNAME)
