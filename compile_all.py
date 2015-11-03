import sys
import os
from os import path
import subprocess 
from argparse import ArgumentParser

DEBUG = False

basedir = path.split(path.abspath(__file__))[0]

CONFIG = {
    "BASE": basedir,
    "LOCALDIR": path.join(basedir, 'local'),
    "BINDIR": path.join(basedir, 'bin'),
    "SRCDIR": path.join(basedir, 'src'),
    "CORES": 4,
    }

def compile_fasttree():
    cmds = """(
    rm -f %(BINDIR)s/FastTree;
    cd %(SRCDIR)s;
    cd FastTree2;
    make clean;
    make -j %(CORES)s && cp FastTree %(BINDIR)s;
    make clean; 
    ls %(BINDIR)s/FastTree;
    ) > %(BASE)s/fasttree.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_raxml():
    cmds = """(
    rm -f %(BINDIR)s/raxml*;
    cd %(SRCDIR)s;
    cd standard-RAxML;
    make clean;
    for makefile in Makefile.SSE3.gcc Makefile.SSE3.PTHREADS.gcc;
    do
       cp $makefile Makefile && rm -f *.o && make -j %(CORES)s;
    done;
    cp raxmlHPC-SSE3 raxmlHPC-PTHREADS-SSE3 %(BINDIR)s;
    ls %(BINDIR)s/raxmlHPC-SSE3 %(BINDIR)s/raxmlHPC-PTHREADS-SSE3;
    ) >%(BASE)s/raxml.log  2>&1;

    """ %CONFIG
    return system(cmds) == 0

def compile_phyml():
    cmds = """(
    rm -f %(BINDIR)s/phyml;
    cd %(SRCDIR)s;
    cd phyml-20120412-patched;
    make clean;
    ./configure; 
    make -j %(CORES)s;
    cp src/phyml %(BINDIR)s;
    ls %(BINDIR)s/phyml;
    ) >%(BASE)s/phyml.log  2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_tcoffee():
    cmds = """(
    rm -f %(BINDIR)s/t_coffee;
    cd %(SRCDIR)s/T-COFFEE_distribution_Version_10.00.r1613/t_coffee_source/ && make clean && make;
    cd %(SRCDIR)s/T-COFFEE_distribution_Version_10.00.r1613;
    ./install t_coffee m_coffee -plugins=%(BINDIR)s  -dis=%(BINDIR)s -exec=%(BINDIR)s; 
    ls %(BINDIR)s/t_coffee;
    ) >%(BASE)s/tcoffee.log  2>&1;

    """ %CONFIG
    return system(cmds) == 0

def compile_muscle():
    cmds = """(
    rm -f %(BINDIR)s/muscle;
    cd %(SRCDIR)s/muscle3.8.31/src/
    make clean;
    ./mk 
    cp muscle %(BINDIR)s/; 
    ls %(BINDIR)s/muscle;
    ) >%(BASE)s/muscle.log 2>&1;

    """ %CONFIG
    return system(cmds) == 0

def compile_argtable2():
    cmds = """(
    cd %(SRCDIR)s/argtable2-13;
    make clean;
    ./configure --prefix %(LOCALDIR)s;
    make -j %(CORES)s ;
    make install;
    ) >%(BASE)s/argtable2.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0
    

def compile_clustalo():
    for attempt in ['native', 'local_argtable']:
        if attempt == "local_argtable":
            print 'argtable2 library is missing. Attempting to compile a local version'
            CONFIG["FLAGS"] = "CFLAGS='-I%(LOCALDIR)s/include' LDFLAGS='-L%(LOCALDIR)s/lib'" %CONFIG
            compile_argtable2()
        else:
            CONFIG["FLAGS"] = ""
        cmds = """(
        rm -f %(BINDIR)s/clustalo;
        cd %(SRCDIR)s/clustal-omega-1.2.1;
        make clean;
        ./configure %(FLAGS)s;
        make -j %(CORES)s ;
        cp src/clustalo %(BINDIR)s/;
        make clean;
        ls %(BINDIR)s/clustalo;
        ) >%(BASE)s/clustalo.log 2>&1;
        """ %CONFIG
        
        if system(cmds) == 0:
            return True

    return False

def compile_mafft():
    cmds = """(
    rm -f %(BINDIR)s/mafft*;
    cd %(SRCDIR)s/mafft-6.861-without-extensions/core/;
    make clean;
    make -j %(CORES)s ;
    cp ../binaries/* %(BINDIR)s/;
    cp ../scripts/* %(BINDIR)s/;
    make clean;
    ls %(BINDIR)s/mafft
    ) >%(BASE)s/mafft.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_trimal():
    cmds = """(
    rm -f %(BINDIR)s/trimal;
    rm -f %(BINDIR)s/readal;
    rm -f %(BINDIR)s/statal;
    cd %(SRCDIR)s/trimal-1.4/source/
    make clean;
    make -j %(CORES)s ;
    cp trimal statal readal %(BINDIR)s/;
    make clean;
    ls %(BINDIR)s/readal %(BINDIR)s/trimal  %(BINDIR)s/statal;
    ) >%(BASE)s/trimal.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_dialigntx():
    cmds = """(
    rm -f %(BINDIR)s/dialign-tx;
    cd %(SRCDIR)s/DIALIGN-TX_1.0.2/source/
    make clean;
    make -j %(CORES)s ;
    cp dialign-tx %(BINDIR)s/;
    make clean;
    ls %(BINDIR)s/dialign-tx;
    ) >%(BASE)s/dialigntx.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_consel():
    cmds = """(
    rm -f %(BINDIR)s/consel;
    cd %(SRCDIR)s/consel/src/
    make clean;
    make -j %(CORES)s ;
    make install;
    cp ../bin/* %(BINDIR)s/;
    make clean;
    ls %(BINDIR)s/consel;
    ) >%(BASE)s/consel.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0


def compile_paml():
    cmds = """(
    rm %(BINDIR)s/codeml;
    cd %(SRCDIR)s/paml4.8/src/
    rm codeml;
    make codeml -j %(CORES)s ;
    cp codeml %(BINDIR)s/;
    rm codeml;    
    ls %(BINDIR)s/codeml;
    ) >%(BASE)s/paml.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0


def compile_slr():
    cmds = """(   
    rm %(BINDIR)s/Slr;
    cd %(SRCDIR)s/slr/src/
    make clean;
    rm ../bin/Slr;
    make -j %(CORES)s ;
    cp ../bin/Slr %(BINDIR)s/;
    rm ../bin/Slr;    
    ls %(BINDIR)s/Slr;
    ) >%(BASE)s/slr.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0



def compile_all(targets=None):
    if not targets:
        targets= ['fasttree', 'raxml', 'phyml', 'tcoffee', 'trimal', 'clustalo', 'muscle', 'dialigntx', 'mafft', 'consel', 'paml', 'slr']
   
    fn = globals()
    for name in targets:
        print >>sys.stderr, 'Compiling', name, "...",
        if not fn.get("compile_%s"%name)():
            if ARGS.verbose:
                print >>sys.stderr, "ERROR\nCompiling %s. Check log %s/%s.log" %(name, CONFIG["BASE"], name)
                logfile = open("%s/%s.log" %(CONFIG["BASE"], name)).read()
                print >>sys.stderr, logfile
            else:
                print >>sys.stderr, "ERROR\nCompiling %s. Check log %s/%s.log" %(name, CONFIG["BASE"], name)
        else:     
            print >>sys.stderr, "Ok"

def system(cmd):
    if ARGS.debug:
        print cmd
    return os.system(cmd)

            
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-v", dest="verbose", action="store_true")
    parser.add_argument("--debug", dest="debug", action="store_true")
    parser.add_argument(dest="targets", nargs="*")

    ARGS = parser.parse_args()
    compile_all(ARGS.targets)
