from __future__ import absolute_import
from __future__ import print_function

import sys
import os
from os import path
import subprocess 
from argparse import ArgumentParser
import multiprocessing
from glob import glob
import shutil

ARGS = None
basedir = path.split(path.abspath(__file__))[0]
DEBUG = False

CONFIG = {
    "BASE": basedir,
    "LOCALDIR": path.join(basedir, 'local'),
    "BINDIR": path.join(basedir, 'bin'),
    "SRCDIR": path.join(basedir, 'src'),
    "CORES":  multiprocessing.cpu_count(),
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
    cd phyml;
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
    cd %(SRCDIR)s/tcoffee/t_coffee_source/ && make clean && make;
    cd %(SRCDIR)s/tcoffee;
    ./install t_coffee -plugins=%(BINDIR)s  -dis=%(BINDIR)s -exec=%(BINDIR)s; 
    ls %(BINDIR)s/t_coffee;
    ) >%(BASE)s/tcoffee.log  2>&1;

    """ %CONFIG
    return system(cmds) == 0

def compile_kalign():
    cmds = """(
    rm -rf %(BINDIR)s/kalign;
    cd %(SRCDIR)s/kalign/;
    make clean;
    ./configure && make -j %(CORES)s; 
    cp kalign %(BINDIR)s/; 
    ls %(BINDIR)s/kalign; 
    ) >%(BASE)s/kalign.log  2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_prank():
    cmds = """(
    rm -rf %(BINDIR)s/prank;
    cd %(SRCDIR)s/prank/;
    make clean;
    rm prank;
    make -j %(CORES)s; 
    cp prank %(BINDIR)s/; 
    ls %(BINDIR)s/prank; 
    ) >%(BASE)s/prank.log  2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_probcons():
    cmds = """(
    rm -rf %(BINDIR)s/probcons;
    cd %(SRCDIR)s/probcons/;
    make clean;
    make -j %(CORES)s; 
    cp probcons %(BINDIR)s/; 
    ls %(BINDIR)s/probcons; 
    ) >%(BASE)s/probcons.log  2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_pmodeltest():
    cmds = """(
    rm -rf %(BINDIR)s/pmodeltest.py;
    cp %(SRCDIR)s/pmodeltest/pmodeltest.py %(BINDIR)s;
    ls %(BINDIR)s/pmodeltest.py; 
    ) >%(BASE)s/pmodeltest.log  2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_muscle():
    cmds = """(
    rm -f %(BINDIR)s/muscle;
    cd %(SRCDIR)s/muscle/src/
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
            print('argtable2 library is missing. Attempting to compile a local version')
            CONFIG["FLAGS"] = "CFLAGS='-I%(LOCALDIR)s/include' LDFLAGS='-L%(LOCALDIR)s/lib'" %CONFIG
            compile_argtable2()
        else:
            CONFIG["FLAGS"] = ""
        cmds = """(
        rm -f %(BINDIR)s/clustalo;
        cd %(SRCDIR)s/clustal-omega;
        make clean;
        ./configure %(FLAGS)s &&
        make -j %(CORES)s ;
        cp src/clustalo %(BINDIR)s/;
        ls %(BINDIR)s/clustalo;
        ) >%(BASE)s/clustalo.log 2>&1;
        """ %CONFIG
        
        if system(cmds) == 0:
            return True

    return False

def compile_mafft():
    cmds = """(
    rm -f %(BINDIR)s/mafft*;
    cd %(SRCDIR)s/mafft/core/;
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
    cd %(SRCDIR)s/trimal/source/
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
    rm -rf %(BINDIR)s/dialigntx_conf/;
    cd %(SRCDIR)s/dialigntx/source/
    make clean;
    make -j %(CORES)s ;
    cp dialign-tx %(BINDIR)s/;
    cp -r %(SRCDIR)s/dialigntx/conf/ %(BINDIR)s/dialigntx_conf/;
    make clean;
    ls %(BINDIR)s/dialign-tx;
    ls %(BINDIR)s/dialigntx_conf/;
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
    cd %(SRCDIR)s/paml/src/
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
    make -j %(CORES)s";
    cp ../bin/Slr %(BINDIR)s/;
    ls %(BINDIR)s/Slr;
    ) >%(BASE)s/slr.log 2>&1;
    """ %CONFIG
    return system(cmds) == 0

def compile_iqtree():
    pass

def compile_phylobayes():
    pass


def compile_all(targets=None, verbose=False, cores=1):
    global CONFIG
    
    CONFIG["cores"] = cores
    if not targets:
        targets= ['tcoffee', 'clustalo', 'muscle', 'dialigntx', 'mafft', 'kalign', 'prank', 'probcons', 
                  'trimal',
                  'pmodeltest',
                  'fasttree', 'raxml', 'phyml', 'iqtree', 'phylobayes',
                  'consel', 'paml', 'slr',
        ]
        if sys.platform == "darwin":
            targets.remove("dialigntx")
        
        for fname in glob(os.path.join(CONFIG["BINDIR"], "*")):
            print("cleaning", fname)
            if os.path.isdir(fname):
                shutil.rmtree(fname)
            else:
                os.remove(fname)
   
    fn = globals()
    errors = 0
    for name in targets:
        print('Compiling', name, "...", end="")
        sys.stdout.flush()
        if not fn.get("compile_%s"%name)():
            if verbose:
                print("ERROR\nCompiling %s. Check log %s/%s.log" %(name, CONFIG["BASE"], name))
                logfile = open("%s/%s.log" %(CONFIG["BASE"], name)).read()
                print(logfile)
            else:
                print("ERROR\nCompiling %s. Check log %s/%s.log" %(name, CONFIG["BASE"], name))
            errors += 1
        else:     
            print("Ok")
    return errors
            
def system(cmd):
    if DEBUG:
        print(cmd)
    return os.system(cmd)

def _main():
    global ARGS, DEBUG
    parser = ArgumentParser()
    parser.add_argument("-v", dest="verbose", action="store_true")
    parser.add_argument("--debug", dest="debug", action="store_true")
    parser.add_argument(dest="targets", nargs="*")

    ARGS = parser.parse_args()
    if ARGS.debug:
        DEBUG = True
    
    errors = compile_all(targets=ARGS.targets, verbose=ARGS.verbose)
    #if errors:
    #    sys.exit(-1)
    
if __name__ == "__main__":
    _main()
