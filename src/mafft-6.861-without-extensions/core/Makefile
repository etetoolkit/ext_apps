PREFIX = /usr/local
LIBDIR = $(PREFIX)/libexec/mafft
BINDIR = $(PREFIX)/bin
MANDIR = $(PREFIX)/share/man/man1

#MNO_CYGWIN = -mno-cygwin

#ENABLE_MULTITHREAD = -Denablemultithread
#Uncomment this to enable multithreading (linux only)

CC = gcc
CFLAGS = -O3
#CFLAGS =  -O0  -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3

MYCFLAGS = $(MNO_CYGWIN) $(ENABLE_MULTITHREAD) $(CFLAGS)

ifeq ($(ENABLE_MULTITHREAD),-Denablemultithread)
LIBS = -lm  -lpthread
else
LIBS = -lm
endif

INSTALL = install

PROGS = dvtditr dndfast7 dndblast sextet5 mafft-distance pairlocalalign \
		pair2hat3s multi2hat3s rnatest pairash \
        splittbfast disttbfast tbfast mafft-profile f2cl mccaskillwrap contrafoldwrap countlen \
		seq2regtable regtable2seq score getlag dndpre dndpre2 setcore replaceu restoreu
SCRIPTS = mafft mafft-homologs.rb
OBJREPLACEU = mtxutl.o io.o replaceu.o defs.o mltaln9.o
OBJRESTOREU = mtxutl.o io.o restoreu.o defs.o mltaln9.o
OBJREGTABLE2SEQ = mtxutl.o io.o regtable2seq.o defs.o mltaln9.o
OBJSEQ2REGTABLE = mtxutl.o io.o seq2regtable.o defs.o
OBJCOUNTLEN = mtxutl.o io.o countlen.o defs.o
OBJF2CL = mtxutl.o io.o f2cl.o constants.o defs.o
OBJMCCASKILLWRAP = mtxutl.o io.o mccaskillwrap.o constants.o defs.o mltaln9.o
OBJCONTRAFOLDWRAP = mtxutl.o io.o contrafoldwrap.o constants.o defs.o mltaln9.o
OBJMULTI2HAT3S = mtxutl.o io.o mltaln9.o tddis.o constants.o \
				 multi2hat3s.o defs.o fft.o fftFunctions.o
OBJPAIR2HAT3S = mtxutl.o io.o mltaln9.o tddis.o constants.o \
				pair2hat3s.o defs.o fft.o fftFunctions.o
OBJRNATEST = mtxutl.o io.o mltaln9.o tddis.o constants.o Lalignmm.o \
			 rna.o rnatest.o defs.o fft.o fftFunctions.o
OBJPAIRASH = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
			 Falign.o MSalignmm.o Galign11.o MSalign11.o suboptalign11.o genGalign11.o genalign11.o Lalign11.o SAalignmm.o \
			 pairash.o defs.o fft.o fftFunctions.o
OBJPAIRLOCALALIGN = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		            Falign.o MSalignmm.o Galign11.o MSalign11.o suboptalign11.o genGalign11.o genalign11.o Lalign11.o SAalignmm.o \
					pairlocalalign.o defs.o fft.o fftFunctions.o
OBJDUMMY = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			disttbfast_dummy.o dummy.o defs.o fft.o fftFunctions.o
OBJSPLITFROMALN = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			Lalign11.o splitfromaln.o defs.o fft.o fftFunctions.o
OBJSPLITTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			Lalign11.o splittbfast.o defs.o fft.o fftFunctions.o
OBJSPLITTBFAST2 = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			Lalign11.o splittbfast2.o defs.o fft.o fftFunctions.o
OBJSPLITTBFASTP = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			Lalign11.o defs.o fft.o fftFunctions.o
OBJDISTTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			disttbfast.o defs.o fft.o fftFunctions.o addfunctions.o
OBJTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o MSalignmm.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o \
			tbfast.o defs.o fft.o fftFunctions.o addfunctions.o
OBJTBFAST2 = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o MSalignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o \
			tbfast2.o defs.o fft.o fftFunctions.o
OBJSETCORE = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			setcore.o defs.o fft.o fftFunctions.o
OBJTDITR = mtxutl.o io.o mltaln9.o tddis.o constants.o nj.o  partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		   Falign.o Falign_localhom.o Galign11.o fftFunctions.o fft.o  \
		   tditeration.o tditr.o defs.o SAalignmm.o treeOperation.o
OBJDVTDITR = mtxutl.o io.o mltaln9.o tddis.o constants.o nj.o  partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
			Falign.o Falign_localhom.o Galign11.o MSalignmm.o fftFunctions.o fft.o  \
			tditeration.o dvtditr.o defs.o SAalignmm.o treeOperation.o
OBJGETLAG = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o MSalignmm.o \
			getlag.o defs.o fft.o fftFunctions.o
OBJGAPFILL = mtxutl.o io.o constants.o gapfill.o defs.o
OBJDNDFAST5 = dndfast5.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDBLAST = dndblast.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST7 = dndfast7.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST6 = dndfast6.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST4 = dndfast4.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST6 = dndfast6.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJSEXTET5 = io.o constants.o mtxutl.o mltaln9.o tddis.o sextet5.o defs.o
OBJDISTANCE = io.o constants.o mtxutl.o mltaln9.o tddis.o mafft-distance.o defs.o
OBJTRIPLET6 = io.o constants.o mtxutl.o mltaln9.o tddis.o triplet6.o defs.o
OBJTRIPLET5 = io.o constants.o mtxutl.o mltaln9.o tddis.o triplet5.o defs.o
OBJOCTET4 = io.o constants.o mtxutl.o mltaln9.o tddis.o octet4.o defs.o
OBJDNDPRE = dndpre.o io.o constants.o mtxutl.o mltaln9.o defs.o
OBJDNDPRESCORE = dndpre2.o io.o constants.o mtxutl.o mltaln9.o defs.o
OBJGALN = io.o mtxutl.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o MSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
          SAalignmm.o Galign11.o Falign.o Falign_localhom.o fftFunctions.o fft.o mafft-profile.o defs.o
OBJSCORE = io.o mtxutl.o mltaln9.o score.o constants.o defs.o

HEADER = mltaln.h mtxutl.h
FFTHEADER = fft.h

MANPAGES = mafft.1 mafft-homologs.1



all : $(PROGS) $(SCRIPTS)
	cp $(SCRIPTS) ../scripts
	chmod 755 ../scripts/*
	cp $(PROGS) ../binaries
	chmod 755 ../binaries/*
	cp $(MANPAGES) ../binaries
	@echo done.

univscript: univscript.tmpl Makefile
	sed "s:_PROGS:$(PROGS):" univscript.tmpl  > univscript

mafft: mafft.tmpl
	sed "s:_LIBDIR:$(LIBDIR):" mafft.tmpl  > mafft

mafft-homologs.rb: mafft-homologs.tmpl
#	cp mafft-homologs.tmpl mafft-homologs.rb
	sed "s:_BINDIR:$(BINDIR):" mafft-homologs.tmpl  > mafft-homologs.rb

mltaln.h : functions.h
	touch mltaln.h

tbfast : $(OBJTBFAST)
	$(CC) -o $@ $(OBJTBFAST) $(MYCFLAGS) $(LIBS)

tbfast2 : $(OBJTBFAST2)
	$(CC) -o $@ $(OBJTBFAST2) $(MYCFLAGS) $(LIBS)

disttbfast : $(OBJDISTTBFAST)
	$(CC) -o $@ $(OBJDISTTBFAST) $(MYCFLAGS) $(LIBS)

splittbfast : $(OBJSPLITTBFAST)
	$(CC) -o $@ $(OBJSPLITTBFAST) $(MYCFLAGS) $(LIBS)

splitfromaln : $(OBJSPLITFROMALN)
	$(CC) -o $@ $(OBJSPLITFROMALN) $(MYCFLAGS) $(LIBS)

splittbfast2 : $(OBJSPLITTBFAST2)
	$(CC) -o $@ $(OBJSPLITTBFAST2) $(MYCFLAGS) $(LIBS)

dummy : $(OBJDUMMY)
	$(CC) -o $@ $(OBJDUMMY) $(MYCFLAGS) $(LIBS)

setcore : $(OBJSETCORE)
	$(CC) -o $@ $(OBJSETCORE) $(MYCFLAGS) $(LIBS)

countlen : $(OBJCOUNTLEN)
	$(CC) -o $@ $(OBJCOUNTLEN) $(MYCFLAGS) $(LIBS)

seq2regtable : $(OBJSEQ2REGTABLE)
	$(CC) -o $@ $(OBJSEQ2REGTABLE) $(MYCFLAGS) $(LIBS)

regtable2seq : $(OBJREGTABLE2SEQ)
	$(CC) -o $@ $(OBJREGTABLE2SEQ) $(MYCFLAGS) $(LIBS)

replaceu : $(OBJREPLACEU)
	$(CC) -o $@ $(OBJREPLACEU) $(MYCFLAGS) $(LIBS)

restoreu : $(OBJRESTOREU)
	$(CC) -o $@ $(OBJRESTOREU) $(MYCFLAGS) $(LIBS)

f2cl : $(OBJF2CL)
	$(CC) -o $@ $(OBJF2CL) $(MYCFLAGS) $(LIBS)

mccaskillwrap : $(OBJMCCASKILLWRAP)
	$(CC) -o $@ $(OBJMCCASKILLWRAP) $(MYCFLAGS) $(LIBS)

contrafoldwrap : $(OBJCONTRAFOLDWRAP)
	$(CC) -o $@ $(OBJCONTRAFOLDWRAP) $(MYCFLAGS) $(LIBS)

pairlocalalign : $(OBJPAIRLOCALALIGN)
	$(CC) -o $@ $(OBJPAIRLOCALALIGN) $(MYCFLAGS) $(LIBS)

pairash : $(OBJPAIRASH)
	$(CC) -o $@ $(OBJPAIRASH) $(MYCFLAGS) $(LIBS)

rnatest : $(OBJRNATEST)
	$(CC) -o $@ $(OBJRNATEST) $(MYCFLAGS) $(LIBS)

pair2hat3s : $(OBJPAIR2HAT3S)
	$(CC) -o $@ $(OBJPAIR2HAT3S) $(MYCFLAGS) $(LIBS)

multi2hat3s : $(OBJMULTI2HAT3S)
	$(CC) -o $@ $(OBJMULTI2HAT3S) $(MYCFLAGS) $(LIBS)

getlag : $(OBJGETLAG)
	$(CC) -o $@ $(OBJGETLAG) $(MYCFLAGS) $(LIBS)

tditr : $(OBJTDITR)
	$(CC) -o $@ $(OBJTDITR) $(MYCFLAGS) $(LIBS)

dvtditr : $(OBJDVTDITR)
	$(CC) -o $@ $(OBJDVTDITR) $(MYCFLAGS) $(LIBS)

mafft-profile : $(OBJGALN)
	$(CC) -o $@ $(OBJGALN) $(MYCFLAGS) $(LIBS)

gapfill : $(OBJGAPFILL)
	$(CC) -o $@ $(OBJGAPFILL) $(MYCFLAGS) $(LIBS)

dndfast4 : $(OBJDNDFAST4)
	$(CC) -o $@ $(OBJDNDFAST4) $(MYCFLAGS) $(LIBS)

dndfast5 : $(OBJDNDFAST5)
	$(CC) -o $@ $(OBJDNDFAST5) $(MYCFLAGS) $(LIBS)

dndfast6 : $(OBJDNDFAST6)
	$(CC) -o $@ $(OBJDNDFAST6) $(MYCFLAGS) $(LIBS)

dndfast7 : $(OBJDNDFAST7)
	$(CC) -o $@ $(OBJDNDFAST7) $(MYCFLAGS) $(LIBS)

dndblast : $(OBJDNDBLAST)
	$(CC) -o $@ $(OBJDNDBLAST) $(MYCFLAGS) $(LIBS)

dndfast3 : $(OBJDNDFAST3)
	$(CC) -o $@ $(OBJDNDFAST3) $(MYCFLAGS) $(LIBS)

triplet : $(OBJTRIPLET)
	$(CC) -o $@ $(OBJTRIPLET) $(MYCFLAGS) $(LIBS)

triplet3 : $(OBJTRIPLET3)
	$(CC) -o $@ $(OBJTRIPLET3) $(MYCFLAGS) $(LIBS)

sextet3 : $(OBJSEXTET3)
	$(CC) -o $@ $(OBJSEXTET3) $(MYCFLAGS) $(LIBS)

sextet4 : $(OBJSEXTET4)
	$(CC) -o $@ $(OBJSEXTET4) $(MYCFLAGS) $(LIBS)

sextet5 : $(OBJSEXTET5)
	$(CC) -o $@ $(OBJSEXTET5) $(MYCFLAGS) $(LIBS)

mafft-distance : $(OBJDISTANCE)
	$(CC) -o $@ $(OBJDISTANCE) $(MYCFLAGS) $(LIBS)

triplet5 : $(OBJTRIPLET5)
	$(CC) -o $@ $(OBJTRIPLET5) $(MYCFLAGS) $(LIBS)

triplet6 : $(OBJTRIPLET6)
	$(CC) -o $@ $(OBJTRIPLET6) $(MYCFLAGS) $(LIBS)

octet4 : $(OBJOCTET4)
	$(CC) -o $@ $(OBJOCTET4) $(MYCFLAGS) $(LIBS)

dndpre : $(OBJDNDPRE)
	$(CC) -o $@ $(OBJDNDPRE) $(MYCFLAGS) $(LIBS)

dndpre2 : $(OBJDNDPRESCORE)
	$(CC) -o $@ $(OBJDNDPRESCORE) $(MYCFLAGS) $(LIBS)

score : $(OBJSCORE)
	$(CC) -o $@ $(OBJSCORE) $(MYCFLAGS) $(LIBS)

genMtx : $(OBJGENMTX)
	$(CC) -o $@ $(OBJGENMTX) $(MYCFLAGS) $(LIBS)

gapfill.o : gapfill.c $(HEADER)
	$(CC) $(MYCFLAGS) -c gapfill.c

mltaln9.o : mltaln9.c $(HEADER)
	$(CC) $(MYCFLAGS) -c mltaln9.c

tddis.o : tddis.c $(HEADER)
	$(CC) $(MYCFLAGS) -c tddis.c

constants.o : constants.c miyata.h miyata5.h blosum.c DNA.h JTT.c $(HEADER)
	$(CC) $(MYCFLAGS) -c constants.c

defs.o : defs.c 
	$(CC) $(MYCFLAGS) -c defs.c

#A+++alignmm.o : SA+++alignmm.c $(HEADER)
#	$(CC) $(MYCFLAGS) -c SA+++alignmm.c -o A+++alignmm.o

Salignmm.o : Salignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Salignmm.c 

Halignmm.o : Halignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Halignmm.c 

Ralignmm.o : Ralignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Ralignmm.c 

Qalignmm.o : Qalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Qalignmm.c 

MSalignmm.o : MSalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c MSalignmm.c 

partSalignmm.o : partSalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c partSalignmm.c  

partQalignmm.o : partQalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c partQalignmm.c  

Lalign11.o : Lalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Lalign11.c 

genalign11.o : genalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c genalign11.c 

genGalign11.o : genGalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c genGalign11.c 

suboptalign11.o : suboptalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c suboptalign11.c 

Galign11.o : Galign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Galign11.c 

MSalign11.o : MSalign11.c $(HEADER)
	$(CC) $(MYCFLAGS) -c MSalign11.c 

SAalignmm.o : SAalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c SAalignmm.c -o SAalignmm.o

Lalignmm.o : Lalignmm.c $(HEADER)
	$(CC) $(MYCFLAGS) -c Lalignmm.c 

rna.o : rna.c $(HEADER)
	$(CC) $(MYCFLAGS) -c rna.c 

disttbfast.o : disttbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c disttbfast.c

splitfromaln.o : splitfromaln.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splitfromaln.c

splittbfast.o : splittbfast.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splittbfast.c

splittbfast2.o : splittbfast2.c $(HEADER) $(FFTHEADER) 
	$(CC) $(MYCFLAGS) -c splittbfast2.c

disttbfast_dummy.o : disttbfast_dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c disttbfast_dummy.c

dummy.o : dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c dummy.c

tbfast.o : tbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c tbfast.c

tbfast2.o : tbfast2.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c tbfast2.c

setcore.o : setcore.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c setcore.c

getlag.o : getlag.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c getlag.c

tditr.o : tditr.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c tditr.c

dvtditr.o : dvtditr.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c dvtditr.c

tditeration.o : tditeration.c $(HEADER)
	 $(CC) $(MYCFLAGS) -c tditeration.c

mafft-profile.o : mafft-profile.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c mafft-profile.c

dndfast4.o : dndfast4.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast4.c

dndfast5.o : dndfast5.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast5.c

dndfast6.o : dndfast6.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast6.c

dndfast7.o : dndfast7.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast7.c

dndblast.o : dndblast.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndblast.c

dndfast3.o : dndfast3.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c dndfast3.c

dndpre.o : dndpre.c $(HEADER)
	$(CC) $(MYCFLAGS) -c dndpre.c

countlen.o : countlen.c $(HEADER)
	$(CC) $(MYCFLAGS) -c countlen.c

seq2regtable.o : seq2regtable.c $(HEADER)
	$(CC) $(MYCFLAGS) -c seq2regtable.c

regtable2seq.o : regtable2seq.c $(HEADER)
	$(CC) $(MYCFLAGS) -c regtable2seq.c

f2cl.o : f2cl.c $(HEADER)
	$(CC) $(MYCFLAGS) -c f2cl.c

replaceu.o : replaceu.c $(HEADER)
	$(CC) $(MYCFLAGS) -c replaceu.c

restoreu.o : restoreu.c $(HEADER)
	$(CC) $(MYCFLAGS) -c restoreu.c

mccaskillwrap.o : mccaskillwrap.c $(HEADER)
	$(CC) $(MYCFLAGS) -c mccaskillwrap.c

contrafoldwrap.o : contrafoldwrap.c $(HEADER)
	$(CC) $(MYCFLAGS) -c contrafoldwrap.c

pairlocalalign.o : pairlocalalign.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pairlocalalign.c

pairash.o : pairash.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pairash.c

rnatest.o : rnatest.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c rnatest.c

multi2hat3s.o : multi2hat3s.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c multi2hat3s.c

pair2hat3s.o : pair2hat3s.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c pair2hat3s.c

dndpre2.o : dndpre2.c $(HEADER)
	$(CC) $(MYCFLAGS) -c dndpre2.c

io.o : io.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c io.c

nj.o : nj.c $(HEADER)
	$(CC) $(MYCFLAGS) -c nj.c

treeOperation.o : treeOperation.c $(HEADER)
	$(CC) $(MYCFLAGS) -c treeOperation.c

sextet5.o : sextet5.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c sextet5.c

mafft-distance.o : mafft-distance.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c mafft-distance.c

triplet5.o : triplet5.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c triplet5.c

triplet6.o : triplet6.c $(HEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c triplet6.c

fft.o : fft.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fft.c 

fftFunctions.o : fftFunctions.c $(HEADER) $(FFTHEADER)
	$(CC) $(MYCFLAGS) -c fftFunctions.c

Falign.o : Falign.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c Falign.c

Falign_localhom.o : Falign_localhom.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(MYCFLAGS) -c Falign_localhom.c

mtxutl.o : mtxutl.c 
	$(CC) $(MYCFLAGS) -c mtxutl.c

addfunctions.o : addfunctions.c $(HEADER)
	$(CC) $(MYCFLAGS) -c addfunctions.c

score.o : score.c $(HEADER)
	$(CC) $(MYCFLAGS) -c score.c

clean :
	rm -f *.o *.a *.exe *~ $(PROGS)  $(SCRIPTS)
#	rm -f ../binaries/* ../scripts/*

install : all
	mkdir -p $(LIBDIR)
	chmod 755 $(LIBDIR)
	mkdir -p $(BINDIR)
	chmod 755 $(BINDIR)
	chmod 755 $(SCRIPTS)
	$(INSTALL) $(SCRIPTS)  $(BINDIR)
	chmod 755 $(PROGS)
	$(INSTALL) $(PROGS) $(LIBDIR)
	$(INSTALL) -m 644 $(MANPAGES) $(LIBDIR)

	( cd $(BINDIR); \
rm -f linsi ginsi einsi fftns fftnsi nwns nwnsi xinsi qinsi; \
rm -f mafft-linsi mafft-ginsi mafft-einsi mafft-fftns mafft-fftnsi mafft-nwns mafft-nwnsi mafft-xinsi mafft-qinsi; \
ln -s mafft linsi; ln -s mafft ginsi; ln -s mafft fftns; \
ln -s mafft fftnsi; ln -s mafft nwns; ln -s mafft nwnsi; \
ln -s mafft einsi; \
ln -s mafft mafft-linsi; ln -s mafft mafft-ginsi; ln -s mafft mafft-fftns; \
ln -s mafft mafft-fftnsi; ln -s mafft mafft-nwns; ln -s mafft mafft-nwnsi; \
ln -s mafft mafft-einsi; ln -s mafft mafft-xinsi; ln -s mafft mafft-qinsi;\
rm -f mafft-profile mafft-profile.exe; ln -s $(LIBDIR)/mafft-profile .; \
rm -f mafft-distance mafft-distance.exe; ln -s $(LIBDIR)/mafft-distance . )

	mkdir -p $(MANDIR)
	chmod 755 $(MANDIR)
	$(INSTALL) -m 644 $(MANPAGES) $(MANDIR)
# remove incorrectly installed manpages by previous versions
#	rm -f /usr/local/man/man1/mafft.1 /usr/local/man/man1/mafft-homologs.1 
