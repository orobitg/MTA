#!/bin/bash

dir=$PWD

# Environment 
SRC=./bin/plugins/src
GZ=./bin/plugins/src/gz
BIN=./bin/plugins

mkdir -p $BIN
mkdir -p $GZ

#T-Coffee
wget -P ${GZ} http://www.tcoffee.org/Packages/Stable/Latest/linux/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz
tar xvzf ${GZ}/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz -C ${SRC}
cd ${SRC}/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64/src
make -f makefile
cp t_coffee ${dir}/${BIN}

cd $dir
#ClustalW2
wget -P ${GZ} http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
tar xvzf ${GZ}/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz -C ${SRC}
cp ${SRC}/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2 ${BIN}

#NORMD
cd $dir
wget -P ${GZ} ftp://ftp-igbmc.u-strasbg.fr/pub/NORMD/norMD1_3.tar.gz
tar xvzf ${GZ}/norMD1_3.tar.gz -C ${SRC}
cd $SRC/normd_noexpat
make -f makefile
cp normd $dir/$BIN/

#MGTREE
cd $dir
cd ${SRC}/mgtree
make -f Makefile
cp bin/mgtree $dir/$BIN

#PHYLIP
cd $dir
wget -P ${GZ} http://evolution.gs.washington.edu/phylip/download/phylip-3.696.tar.gz
tar xvzf ${GZ}/phylip-3.696.tar.gz -C ${SRC}
cd $SRC/phylip-3.696/src
make -f Makefile.unx install
cd $dir
cp $SRC/phylip-3.696/exe/retree $dir/$BIN

#MAFFT
cd $dir
wget -P ${GZ} http://mafft.cbrc.jp/alignment/software/mafft-7.187-without-extensions-src.tgz
tar xvzf ${GZ}/mafft-7.187-without-extensions-src.tgz -C ${SRC}
cd $SRC/mafft-7.187-without-extensions/core
make clean
make -f Makefile
cd $dir
cp ${SRC}/mafft-7.187-without-extensions/scripts/mafft $dir/$BIN
wget -P ${BIN} http://mafft.cbrc.jp/alignment/software/newick2mafft.rb

wget -P ${BIN} http://www.clustal.org/omega/clustalo-1.2.0-Ubuntu-x86_64
mv ${BIN}/clustalo-1.2.0-Ubuntu-x86_64 ${BIN}/clustalo
chmod +x ${BIN}/clustalo

cd $dir
make -f Makefile clean
make -f Makefile.mpi clean
make -f Makefile
make -f Makefile.mpi

