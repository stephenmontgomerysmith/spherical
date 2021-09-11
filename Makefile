CFLAGS=-W -Wall -O3 -ffast-math -march=native -I/usr/local/cuda/include
CC=c++
NVCC=nvcc -O3 #-deviceemu

all: spherical

OBJECT= spherical.o ode-adams-bash-2.o tensor.o get-param.o psidot.o

psidot.o: psidot.cu
	${NVCC} -c psidot.cu

ode-adams-bash-2.o: ode-adams-bash-2.cu
	${NVCC} -c ode-adams-bash-2.cu

derivs.o: derivs.cu
	${NVCC} -c derivs.cu

spherical: ${OBJECT}
	${CC} ${CFLAGS} ${OBJECT} -o spherical -lm -L/usr/local/cuda/lib -lcudart

src:
	perl expand-iterate.pl psidot.conf > psidot-unthreaded.c
	perl expand-iterate.pl psidot-koch.conf > psidot-koch-unthreaded.c
	perl expand-iterate.pl psidot-dd.conf > psidot-dd-unthreaded.c
	perl expand-iterate.pl psidot-dd-2.conf > psidot-dd-2-unthreaded.c
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
	printf "\n" >> reverse-tensor.c
	perl make-reverse-tensor.pl 6 >> reverse-tensor.c

src-threaded:
	perl expand-thread.pl psidot-unthreaded.c > psidot.c
	perl expand-thread.pl psidot-koch-unthreaded.c > psidot-koch.c
	perl expand-thread.pl psidot-dd-unthreaded.c > psidot-dd.c
	perl expand-thread.pl psidot-dd-2-unthreaded.c > psidot-dd-2.c
	perl expand-thread.pl psidot-vd-unthreaded.c > psidot-vd.c
	perl expand-thread.pl psidot-vl-unthreaded.c > psidot-vl.c
	perl expand-thread.pl psidot-ard-unthreaded.c > psidot-ard.c

src-unthreaded:
	cp psidot-unthreaded.c psidot.c
	cp psidot-koch-unthreaded.c psidot-koch.c
	cp psidot-dd-unthreaded.c psidot-dd.c
	cp psidot-dd-2-unthreaded.c psidot-dd-2.c
	cp psidot-vd-unthreaded.c psidot-vd.c
	cp psidot-vl-unthreaded.c psidot-vl.c
	cp psidot-ard-unthreaded.c psidot-ard.c

src-maxima:
	env USE_MAXIMA=1 make src

clean:
	rm -f spherical *.o *~ *.aux *.dvi *.log *core s.out s-out.eps misc/*.aux misc/*.dvi misc/*.log lapack/*.o *.bak *.linkinfo

cleansrc:
	rm -f reverse-tensor.c tensor.c psidot.c psidot-*.c *-unthreaded.c
