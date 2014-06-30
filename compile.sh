#!/bin/bash 
TARGET=$1

MYDIR="$(dirname "$(readlink -f "${0}")")"
export BASEDIR=$MYDIR/src
export PORTABLE_DIR=$MYDIR/bin
export MAKE_FLAGS="$2"

cd $BASEDIR
echo $TARGET

if [ $TARGET == "clean" ]; then
    rm $PORTABLE_DIR/*
fi 

if  [  $TARGET == "fasttree" ] || [ $TARGET == "all" ]; then
#FASTTREE
    cd FastTree2
    make clean
    make $MAKE_FLAGS && cp ./FastTree $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

if [ $TARGET == "raxml" ] || [ $TARGET == "all" ]; then
#RAXML
    cd RAxML_dev;
    make clean
    for x in Makefile.gcc Makefile.PTHREADS.gcc Makefile.SSE3.gcc Makefile.SSE3.PTHREADS.gcc Makefile.AVX.gcc Makefile.AVX.PTHREADS.gcc Makefile.AVX2.gcc Makefile.AVX2.PTHREADS.gcc;
    do 
        cp $x Makefile;
        make $MAKE_FLAGS && find -name 'raxml*' -executable -exec cp {} $PORTABLE_DIR/ \;
        make  clean;
    done;
    cd $BASEDIR
fi

if [ $TARGET == "clustalo" ] || [ $TARGET == "all" ]; then
#CLUSTAL OMEGA
    cd clustal-omega-1.2.1
    make clean
    ./configure
    make $MAKE_FLAGS && cp ./src/clustalo $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

if [ $TARGET == "tcoffee" ] || [ $TARGET == "all" ]; then    
#TCOFFEE
#cd t_coffee_9_01/src/
#make clean
#make $MAKE_FLAGS && cp ./t_coffee $PORTABLE_DIR/
#make clean
#cd $BASEDIR
    
#TCOFFEE 10
    cd T-COFFEE_distribution_Version_10.00.r1613
    ./install t_coffee m_coffee -plugins=../../bin/ -dis=../../bin/ -exec=../../bin/ 
    ./install clean
    cd $BASEDIR
fi

    
if [ $TARGET == "dialigntx" ] || [ $TARGET == "all" ]; then    
#DIALIGN-tx
    cd DIALIGN-TX_1.0.2/source/
    make clean
    make $MAKE_FLAGS && cp ./dialign-tx $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

if [ $TARGET == "muscle" ] || [ $TARGET == "all" ]; then    
# MUSCLE
    cd muscle3.8.31/src/
    make clean
    ./mk && cp muscle $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

    
if [ $TARGET == "phyml" ] || [ $TARGET == "all" ]; then    
# PHYML
    cd phyml-20120412-patched
    make clean
    ./configure && make $MAKE_FLAGS && cp src/phyml $PORTABLE_DIR;
    make clean 
    cd $BASEDIR
fi

if [ $TARGET == "mafft" ] || [ $TARGET == "all" ]; then    
# MAFFT
    cd mafft-6.861-without-extensions/core/
    make clean
    make $MAKE_FLAGS
    cp ../binaries/* $PORTABLE_DIR
    cp ../scripts/* $PORTABLE_DIR
    make clean
    cd $BASEDIR
fi

if [ $TARGET == "trimal" ] || [ $TARGET == "all" ]; then    
# TRIMAL
    cd trimal-1.4/source/
    make clean
    make $MAKE_FLAGS && cp trimal readal statal $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

if [ $TARGET == "consel" ] || [ $TARGET == "all" ]; then        
# CONSEL
    cd consel/src/
    make clean
    make $MAKE_FLAGS && make install && cp ../bin/* $PORTABLE_DIR/
    make clean
    cd $BASEDIR
fi

