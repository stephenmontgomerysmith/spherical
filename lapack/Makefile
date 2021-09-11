CFLAGS=-W -Wall -O3 -ffast-math -march=native

all: diagonalize-sym.o

test: diagonalize-sym.c
	perl -lne 'print if !(/\/\*/ || /\*\//)' diagonalize-sym.c > test.c
	${CC} ${CFLAGS} test.c -o test -L/usr/local/lib -L/usr/local/lib/gcc-4.3.3 -llapack -lblas -lm -pthread -lgfortran
	rm test.c

clean:
	rm -f *.o *~ *core test
