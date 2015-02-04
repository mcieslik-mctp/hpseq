CC=gcc
CFLAGS=-g -Wall -Wextra -O3 -Wno-unused-function  -std=c99 -pedantic -pthread 

KLIB_PATH=klib
KLIB_INCL=-I$(KLIB_PATH)

GREP_PATH=grep
GREP_INCL=-I$(GREP_PATH) -I$(GREP_PATH)/lib -I$(GREP_PATH)/src
GREP_LINK=-L$(GREP_PATH)/lib
HTSLIB_PATH=htslib
HTSLIB_LINK=-L$(HTSLIB_PATH)
HTSLIB_INCL=-I$(HTSLIB_PATH)/htslib

all:hpcut hpmerge hpscan_ss hpscan_sw hpscan_cw

clean:
	rm -fr gmon.out *.o ext/*.o a.out  *~ *.a *.dSYM session*
	rm hpcut hpmerge hpscan_ss hpscan_sw

hpcut:hpcut.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) hpcut.c -o $@ -lz -lm

hpmerge:hpmerge.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) hpmerge.c -o $@ -lz -lm

hpscan_ss:hpscan_ss.c
	$(CC) $(CFLAGS) hpscan_ss.c -o $@ -lz -lm

hpscan_cw:hpscan_cw.c $(HTSLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) $(GREP_INCL) $(GREP_LINK) -o $@ hpscan_cw.c grep/src/kwset.c -lgreputils -lz

hpscan_sw:hpscan_sw.c $(HTSLIB_PATH) $(KLIB_PATH)
	$(CC) $(CFLAGS) $(HTSLIB_INCL) $(KLIB_INCL) $(KLIB_PATH)/ksw.c hpscan_sw.c -o $@ -lz -lm

deps: htslib klib grep

htslib:
	git clone "https://github.com/samtools/htslib.git"
	cd htslib && $(MAKE)

klib:
	git clone "https://github.com/attractivechaos/klib.git"

grep:
	wget "http://ftp.gnu.org/gnu/grep/grep-2.21.tar.xz"
	tar xf grep-2.21.tar.xz
	rm grep-2.21.tar.xz
	mv grep-2.21 grep
	cd grep && ./configure && $(MAKE)

clean-deps:
	rm -rf htslib
	rm -rf grep
	rm -rf klib
