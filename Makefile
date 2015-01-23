CC=gcc
CFLAGS=-g -Wall -O1 -Wno-unused-function -Wno-unused-variable -std=c99 -pedantic -pthread 

KLIB_PATH=klib
KLIB_INCL=-I$(KLIB_PATH)

HTSLIB_PATH=htslib
HTSLIB_LINK=-L$(HTSLIB_PATH)
HTSLIB_INCL=-I$(HTSLIB_PATH)/htslib

all:hpcut hpscan_ss hpscan_sw sam_skel


sam_skel:sam_skel.c htslib/libhts.a
	$(CC) $(CFLAGS) $(HTSLIB_INCL)  -o $@ sam_skel.c htslib/libhts.a -lz

hpcut:hpcut.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) hpcut.c -o $@ -lz -lm

hpscan_ss:hpscan_ss.c
	$(CC) $(CFLAGS) hpscan_ss.c -o $@ -lz -lm

hpscan_sw:hpscan_sw.c $(HTSLIB_PATH) $(KLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) $(KLIB_INCL) $(KLIB_PATH)/ksw.c hpscan_sw.c -o $@ -lz -lm

htslib-static:
	cd htslib && $(MAKE)

clean:
	rm -fr gmon.out *.o ext/*.o a.out  *~ *.a *.dSYM session* hpcut
	cd htslib ; make clean
