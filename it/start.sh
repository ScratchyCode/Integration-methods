#!/bin/bash
# Coded by Pietro Squilla & Francesco Servilio
# Simula il moto di un pendolo generale (smorzato forzato)

# entro dentro la cartella con gli eseguibili
DIR=$(dirname "$0")
cd $DIR

killall gnuplot_x11 -q 
killall gnuplot_qt -q 

# compilo ed eseguo il programma
gcc integratore.c -o integratore.exe -lm -pedantic -Wall -ffast-math -O2

date

time ./integratore.exe &

wait $!
rm integratore.exe

# elimina le immagini vuote
find -size 0 -name "*.png" -delete

exit 0
