CC    = gcc
FLAGS = -Wall -O3

all:
	$(CC) $(FLAGS) -finline-functions -funroll-loops -o FastTree  FastTree.c -DOPENMP -fopenmp  -lm -DUSE_DOUBLE
clean:
	-rm -f FastTree *.o *~
