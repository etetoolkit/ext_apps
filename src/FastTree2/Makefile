CC    = gcc
FLAGS = -Wall -O3

all:
	$(CC) $(FLAGS) -finline-functions -funroll-loops -o FastTree  FastTree.c -DOPENMP -fopenmp  -lm
clean:
	-rm -f FastTree *.o *~
