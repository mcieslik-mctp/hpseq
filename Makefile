CC=gcc
CFLAGS=-g -Wall -O1 -Wno-unused-function -std=c99 -pedantic

all:hpcut hpscan_ss hpscan_sw

hpcut:hpcut.c kseq.h
	$(CC) $(CFLAGS) hpcut.c -o $@ -lz -lm

hpscan_ss:hpscan_ss.c
	$(CC) $(CFLAGS) -Wno-unused-variable hpscan_ss.c -o $@ -lz -lm

hpscan_sw:hpscan_sw.c kseq.h ksw.h
	$(CC) $(CFLAGS) ksw.c hpscan_sw.c -o $@ -lz -lm

clean:
	rm -fr gmon.out *.o ext/*.o a.out  *~ *.a *.dSYM session*


