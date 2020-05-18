# build to hurst-exponent
CC = gcc
CFLAGS = -g -Wall -lgsl -lgslcblas -lm
NAME = hurst-exponent
DEPS = main.c

all: $(DEPS)
	$(CC) $(CFLAGS) -o $(NAME) $(DEPS)

.PHONY: clean

clean: 
	$(RM) $(NAME)
