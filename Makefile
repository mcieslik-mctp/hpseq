CC=gcc
CFLAGS=-g -Wall -O1 -Wno-unused-function -Wno-unused-variable -std=c99 -pedantic -pthread 

KLIB_PATH=klib
KLIB_INCL=-I$(KLIB_PATH)

GREP_PATH=grep
GREP_INCL=-I$(GREP_PATH) -I$(GREP_PATH)/lib -I$(GREP_PATH)/src
GREP_LINK=-L$(GREP_PATH)/lib
HTSLIB_PATH=htslib
HTSLIB_LINK=-L$(HTSLIB_PATH)
HTSLIB_INCL=-I$(HTSLIB_PATH)/htslib

all:hpcut hpscan_ss hpscan_sw hpscan_cw sam_skel bed_skel grep_skel


sam_skel:sam_skel.c htslib/libhts.a
	$(CC) $(CFLAGS) $(HTSLIB_INCL) -o $@ sam_skel.c htslib/libhts.a -lz

bed_skel:bed_skel.c htslib/libhts.a
	$(CC) $(CFLAGS) $(HTSLIB_INCL) -o $@ bed_skel.c htslib/libhts.a -lz

grep_skel:grep_skel.c
	$(CC) $(CFLAGS) $(GREP_INCL) $(GREP_LINK) -o grep_skel grep_skel.c grep/src/kwset.c -lgreputils

hpcut:hpcut.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) hpcut.c -o $@ -lz -lm

hpscan_ss:hpscan_ss.c
	$(CC) $(CFLAGS) hpscan_ss.c -o $@ -lz -lm

hpscan_cw:hpscan_cw.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) $(GREP_INCL) $(GREP_LINK) -o $@ hpscan_cw.c grep/src/kwset.c -lgreputils -lz

hpscan_sw:hpscan_sw.c $(HTSLIB_PATH) $(KLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) $(KLIB_INCL) $(KLIB_PATH)/ksw.c hpscan_sw.c -o $@ -lz -lm

htslib:
	git clone
	cd htslib
	$(MAKE)

klib:
	cd klib
	$(MAKE)

grep:
	cd grep
	./configure
	$(MAKE)

clean:
	rm -fr gmon.out *.o ext/*.o a.out  *~ *.a *.dSYM session*
	rm hpcut hpscan_ss hpscan_sw sam_skel bed_skel grep_skel
