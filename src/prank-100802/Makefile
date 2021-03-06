####### Compiler, tools and options

CC            = gcc
CXX           = g++
CFLAGS        = -pipe -O2
CXXFLAGS      = -pipe -O2
#INCPATH       = -I/usr/include -I.
INCPATH       = -I.
LINK          = g++
#LIBS          = -L/usr/lib
DEL_FILE      = rm -f

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = ancestralnode.cpp \
		ancestralsequence.cpp \
		boolmatrix.cpp \
		chaosanchors.cpp \
		characterprobability.cpp \
		dbmatrix.cpp \
		eigen.cpp \
		flmatrix.cpp \
		fullprobability.cpp \
		guidetree.cpp \
		hirschberg.cpp \
		hmmodel.cpp \
		intmatrix.cpp \
		node.cpp \
		phylomatchscore.cpp \
		postprobability.cpp \
		prank.cpp \
		progressivealignment.cpp \
		pwhirschberg.cpp \
		pwsite.cpp \
		readalignment.cpp \
		readannotation.cpp \
		readfile.cpp \
		readnewick.cpp \
		sequence.cpp \
		site.cpp \
		terminalnode.cpp \
		terminalsequence.cpp \
		translatesequences.cpp \
		treenode.cpp \
		writefile.cpp 
OBJECTS       = ancestralnode.o \
		ancestralsequence.o \
		boolmatrix.o \
		chaosanchors.o \
		characterprobability.o \
		dbmatrix.o \
		eigen.o \
		flmatrix.o \
		fullprobability.o \
		guidetree.o \
		hirschberg.o \
		hmmodel.o \
		intmatrix.o \
		node.o \
		phylomatchscore.o \
		postprobability.o \
		prank.o \
		progressivealignment.o \
		pwhirschberg.o \
		pwsite.o \
		readalignment.o \
		readannotation.o \
		readfile.o \
		readnewick.o \
		sequence.o \
		site.o \
		terminalnode.o \
		terminalsequence.o \
		translatesequences.o \
		treenode.o \
		writefile.o

TARGET        = prank

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) -o $(TARGET) $(OBJECTS) $(LIBS)

clean:clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core

####### Compile

ancestralnode.o: ancestralnode.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		writefile.h \
		hirschberg.h \
		phylomatchscore.h \
		terminalsequence.h \
		fullprobability.h \
		postprobability.h \
		characterprobability.h \
		terminalnode.h \
		readalignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ancestralnode.o ancestralnode.cpp

ancestralsequence.o: ancestralsequence.cpp ancestralsequence.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		config.h \
		hmmodel.h \
		ancestralnode.h \
		treenode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o ancestralsequence.o ancestralsequence.cpp

boolmatrix.o: boolmatrix.cpp boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o boolmatrix.o boolmatrix.cpp

chaosanchors.o: chaosanchors.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		chaosanchors.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o chaosanchors.o chaosanchors.cpp

characterprobability.o: characterprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		characterprobability.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o characterprobability.o characterprobability.cpp

dbmatrix.o: dbmatrix.cpp dbmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dbmatrix.o dbmatrix.cpp

eigen.o: eigen.cpp eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o eigen.o eigen.cpp

flmatrix.o: flmatrix.cpp flmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o flmatrix.o flmatrix.cpp

fullprobability.o: fullprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		fullprobability.h \
		phylomatchscore.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o fullprobability.o fullprobability.cpp

guidetree.o: guidetree.cpp guidetree.h \
		flmatrix.h \
		intmatrix.h \
		pwhirschberg.h \
		pwsite.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o guidetree.o guidetree.cpp

hirschberg.o: hirschberg.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		chaosanchors.h \
		hirschberg.h \
		phylomatchscore.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hirschberg.o hirschberg.cpp

hmmodel.o: hmmodel.cpp hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		eigen.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hmmodel.o hmmodel.cpp

intmatrix.o: intmatrix.cpp intmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o intmatrix.o intmatrix.cpp

node.o: node.cpp node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o node.o node.cpp

phylomatchscore.o: phylomatchscore.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		phylomatchscore.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o phylomatchscore.o phylomatchscore.cpp

postprobability.o: postprobability.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		postprobability.h \
		phylomatchscore.h \
		terminalsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o postprobability.o postprobability.cpp

prank.o: prank.cpp progressivealignment.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		ancestralsequence.h \
		translatesequences.h \
		prank.h \
		hmmodel.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o prank.o prank.cpp

progressivealignment.o: progressivealignment.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		readnewick.h \
		readfile.h \
		writefile.h \
		guidetree.h \
		progressivealignment.h \
		translatesequences.h \
		readannotation.h \
		hirschberg.h \
		phylomatchscore.h \
		terminalsequence.h \
		readalignment.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o progressivealignment.o progressivealignment.cpp

pwhirschberg.o: pwhirschberg.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		pwhirschberg.h \
		pwsite.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pwhirschberg.o pwhirschberg.cpp

pwsite.o: pwsite.cpp pwsite.h \
		flmatrix.h \
		intmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pwsite.o pwsite.cpp

readalignment.o: readalignment.cpp readalignment.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		treenode.h \
		phylomatchscore.h \
		ancestralsequence.h \
		terminalsequence.h \
		config.h \
		hmmodel.h \
		ancestralnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readalignment.o readalignment.cpp

readannotation.o: readannotation.cpp readannotation.h \
		flmatrix.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readannotation.o readannotation.cpp

readfile.o: readfile.cpp readfile.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readfile.o readfile.cpp

readnewick.o: readnewick.cpp readnewick.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		ancestralnode.h \
		ancestralsequence.h \
		node.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o readnewick.o readnewick.cpp

sequence.o: sequence.cpp sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o sequence.o sequence.cpp

site.o: site.cpp site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		hmmodel.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o site.o site.cpp

terminalnode.o: terminalnode.cpp terminalnode.h \
		treenode.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		terminalsequence.h \
		config.h \
		hmmodel.h \
		ancestralnode.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o terminalnode.o terminalnode.cpp

terminalsequence.o: terminalsequence.cpp terminalsequence.h \
		sequence.h \
		site.h \
		intmatrix.h \
		flmatrix.h \
		dbmatrix.h \
		boolmatrix.h \
		config.h \
		hmmodel.h \
		ancestralnode.h \
		treenode.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o terminalsequence.o terminalsequence.cpp

translatesequences.o: translatesequences.cpp translatesequences.h \
		config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o translatesequences.o translatesequences.cpp

treenode.o: treenode.cpp config.h \
		hmmodel.h \
		dbmatrix.h \
		flmatrix.h \
		intmatrix.h \
		ancestralnode.h \
		treenode.h \
		sequence.h \
		site.h \
		boolmatrix.h \
		ancestralsequence.h \
		writefile.h \
		hirschberg.h \
		phylomatchscore.h \
		terminalsequence.h \
		fullprobability.h \
		postprobability.h \
		characterprobability.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o treenode.o treenode.cpp

writefile.o: writefile.cpp writefile.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o writefile.o writefile.cpp
