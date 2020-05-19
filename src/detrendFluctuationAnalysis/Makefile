# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
NAME = ../DFA
DEPS = detrendFluctuationAnalysis.c

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(NAME) $(DEPS)

.PHONY: clean

clean: 
	$(RM) $(NAME)
