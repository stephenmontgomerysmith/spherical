CFLAGS=-W -Wall -O3 -ffast-math -march=native

all: spherical

OBJECT= spherical.o derivs.o ode-adams-bash-2.o ode-adams-bash-4.o ode-rk-4.o ode-rkf-23.o ode-rkf-45.o tensor.o get-param.o evaluate-string.o psidot.o psidot-koch.o psidot-dd.o psidot-vd.o psidot-vl.o reverse-tensor.o diagonalize-sym.o rsc.o psidot-ard.o

spherical: ${OBJECT}
	${CC} ${CFLAGS} ${OBJECT} -o spherical -lm -pthread

src:
	perl expand-iterate.pl psidot.conf > psidot-unthreaded.c
	perl expand-iterate.pl psidot-koch.conf > psidot-koch-unthreaded.c
	perl expand-iterate.pl psidot-dd.conf > psidot-dd-unthreaded.c
	perl expand-iterate.pl psidot-vd.conf > psidot-vd-unthreaded.c
	perl expand-iterate.pl psidot-vl.conf > psidot-vl-unthreaded.c
	perl expand-iterate.pl psidot-ard.conf > psidot-ard-unthreaded.c
	printf "#include \"spherical.h\"\n\n" > tensor.c
	perl make-tensor.pl 2 >> tensor.c
	printf "\n" >> tensor.c
	perl make-tensor.pl 4 >> tensor.c
	printf "\n" >> tensor.c
	perl make-tensor.pl 6 >> tensor.c
	printf "#include \"spherical.h\"\n\n" > reverse-tensor.c
	perl make-reverse-tensor.pl 2 >> reverse-tensor.c
	printf "\n" >> reverse-tensor.c
	perl make-reverse-tensor.pl 4 >> reverse-tensor.c

src-threaded:
	perl expand-thread.pl psidot-unthreaded.c > psidot.c
	perl expand-thread.pl psidot-koch-unthreaded.c > psidot-koch.c
	perl expand-thread.pl psidot-dd-unthreaded.c > psidot-dd.c
	perl expand-thread.pl psidot-vd-unthreaded.c > psidot-vd.c
	perl expand-thread.pl psidot-vl-unthreaded.c > psidot-vl.c
	perl expand-thread.pl psidot-ard-unthreaded.c > psidot-ard.c

src-unthreaded:
	cp psidot-unthreaded.c psidot.c
	cp psidot-koch-unthreaded.c psidot-koch.c
	cp psidot-dd-unthreaded.c psidot-dd.c
	cp psidot-vd-unthreaded.c psidot-vd.c
	cp psidot-vl-unthreaded.c psidot-vl.c
	cp psidot-ard-unthreaded.c psidot-ard.c

src-maxima:
	env USE_MAXIMA=1 make src

clean:
	rm -f spherical *.o *~ *.aux *.dvi *.log *core s.out s-out.eps misc/*.aux misc/*.dvi misc/*.log lapack/*.o *.bak

cleansrc:
	rm -f reverse-tensor.c tensor.c psidot.c psidot-*.c *-unthreaded.c
