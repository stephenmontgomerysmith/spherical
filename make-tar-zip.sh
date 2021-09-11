#!/bin/sh

make clean

perl -i -p -e 's/\r*\n/\r\n/g' parameters*.txt

for t in *.tex; do \
  pdflatex $t
  pdflatex $t
done

rm -f *.exe

make src-unthreaded
cp spherical.c spherical-unthreaded.c
perl -i -e 'while(<>){if(/NR_THREADS/){<>;<>;<>;<>;<>;$_=<>}print if!/nr_threads/||/ignore/}' spherical-unthreaded.c
sed -e s/-pthread// -e s/spherical/spherical-unthreaded/ -e s/-march=native// Makefile | make CC=mingw32-gcc -f -
rm spherical-unthreaded.c
mv spherical-unthreaded spherical-unthreaded.exe

sleep 1
make src-threaded
sed -e s/-pthread/-lpthread/ -e s/-march=native// Makefile | make CC=mingw32-gcc -f -
mv spherical spherical.exe

make clean

PWD=`pwd`
DIR=`basename $PWD`

cd ..
tar cvfz $DIR.tar.gz $DIR
rm -f $DIR.zip
zip -r $DIR.zip $DIR
