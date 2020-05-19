# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
NAME = ../HE
DEPS = HurstExponent.c

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(NAME) $(DEPS)

.PHONY: clean

clean: 
	$(RM) $(NAME)
