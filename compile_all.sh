MYDIR="$(dirname "$(readlink -f "${0}")")"
export BASEDIR=$MYDIR/src
export PORTABLE_DIR=$MYDIR/bin
export MAKE_FLAGS="-j3"
cd $BASEDIR
rm $PORTABLE_DIR/*

#FASTTREE
cd FastTree2
make clean
make $MAKE_FLAGS && cp ./FastTree $PORTABLE_DIR/
make clean
cd $BASEDIR

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

#CLUSTAL OMEGA
cd clustal-omega-1.2.1
make clean
./configure
make $MAKE_FLAGS && cp ./src/clustalo $PORTABLE_DIR/
make clean
cd $BASEDIR

#TCOFFEE
cd t_coffee_9_01/src/
make clean
make $MAKE_FLAGS && cp ./t_coffee $PORTABLE_DIR/
make clean
cd $BASEDIR

#DIALIGN-tx
cd DIALIGN-TX_1.0.2/source/
make clean
make $MAKE_FLAGS && cp ./dialign-tx $PORTABLE_DIR/
make clean
cd $BASEDIR

# MUSCLE
cd muscle3.8.31/src/
make clean
./mk && cp muscle $PORTABLE_DIR/
make clean
cd $BASEDIR

# PHYML
cd phyml-20120412-patched
make clean
./configure && make $MAKE_FLAGS && cp src/phyml $PORTABLE_DIR;
make clean 
cd $BASEDIR

# MAFFT
cd mafft-6.861-without-extensions/core/
make clean
make $MAKE_FLAGS
cp ../binaries/* $PORTABLE_DIR
cp ../scripts/* $PORTABLE_DIR
make clean
cd $BASEDIR

# TRIMAL
cd trimal-1.4/source/
make clean
make $MAKE_FLAGS && cp trimal readal statal $PORTABLE_DIR/
make clean
cd $BASEDIR

# CONSEL
cd consel/src/
make clean
make $MAKE_FLAGS && make install && cp ../bin/* $PORTABLE_DIR/
make clean
cd $BASEDIR
