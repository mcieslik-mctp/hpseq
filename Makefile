CC=gcc
CFLAGS=-g -Wall -O3 -Wno-unused-function

all:thecut

thecut:thecut.c kseq.h
		$(CC) $(CFLAGS) thecut.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out thecut *~ *.a *.dSYM session*
