import sys
import os
from os import path
import subprocess 

basedir = path.split(path.abspath(__file__))[0]
CONFIG = {
    "BASE": basedir,
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
    return os.system(cmds) == 0

def compile_raxml():
    cmds = """(
    rm -f %(BINDIR)s/raxml*;
    cd %(SRCDIR)s;
    cd RAxML_dev;
    make clean;
    for makefile in Makefile.SSE3.gcc Makefile.SSE3.PTHREADS.gcc;
    do
       cp $makefile Makefile && make -j %(CORES)s;
    done;
    cp raxmlHPC-SSE3 raxmlHPC-PTHREADS-SSE3 %(BINDIR)s;
    ls %(BINDIR)s/raxmlHPC-SSE3 %(BINDIR)s/raxmlHPC-PTHREADS-SSE3;
    ) >%(BASE)s/raxml.log  2>&1;

    """ %CONFIG
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

def compile_tcoffee():
    cmds = """(
    rm -f %(BINDIR)s/t_coffee;
    cd %(SRCDIR)s/T-COFFEE_distribution_Version_10.00.r1613/t_coffee_source/ && make clean && make;
    cd %(SRCDIR)s/T-COFFEE_distribution_Version_10.00.r1613;
    ./install t_coffee m_coffee -plugins=%(BINDIR)s  -dis=%(BINDIR)s -exec=%(BINDIR)s; 
    ls %(BINDIR)s/t_coffee;
    ) >%(BASE)s/tcoffee.log  2>&1;

    """ %CONFIG
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

def compile_clustalo():
    cmds = """(
    rm -f %(BINDIR)s/clustalo;
    cd %(SRCDIR)s/clustal-omega-1.2.1;
    make clean;
    ./configure;
    make -j %(CORES)s ;
    cp src/clustalo %(BINDIR)s/;
    make clean;
    ls %(BINDIR)s/clustalo;
    ) >%(BASE)s/clustalo.log 2>&1;
    """ %CONFIG
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

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
    return os.system(cmds) == 0

def compile_all(targets = None):
    if not targets:
        targets= ['fasttree', 'raxml', 'phyml', 'tcoffee', 'trimal', 'clustalo', 'muscle', 'dialigntx', 'mafft', 'consel']
   
    fn = globals()
    for name in targets:
        print >>sys.stderr, 'Compiling', name, "...",
        if not fn.get("compile_%s"%name)():
            print >>sys.stderr, "ERROR\nCompiling %s. Check log %s/%s.log" %(name, CONFIG["BASE"], name)
        else:     
            print >>sys.stderr, "Ok"


if __name__ == "__main__":
    try:
        targets = sys.argv[1:]
    except: 
        compile_all()
    else: 
        compile_all(targets)
