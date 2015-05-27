#!bin/sh
clean:
	rm *.qsub.o* tdqmcurrd.dat embound embound-v voltage.dat
main:
	qsub main.qsub
td:
	qsub maintd.qsub
iv:
	qsub mainiv.qsub
