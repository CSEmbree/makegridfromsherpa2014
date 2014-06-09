export BASEDIR=`pwd`
#
#-----------------------
#  setup for SLC6
#
#export platform=x86_64-slc5-gcc43-opt
export platform=x86_64-slc6-gcc46-opt
#
export external=/afs/cern.ch/sw/lcg/external
export releases=/afs/cern.ch/sw/lcg/app/releases
export contrib=/afs/cern.ch/sw/lcg/contrib
#
#export GCC_INSTALL_DIR=${contrib}/gcc/4.3.2/${platform}
#export GCC_INSTALL_DIR=${contrib}/gcc/4.6/${platform}
#source ${contrib}/gcc/4.3.2/x86_64-slc5/setup.sh
#source ${contrib}/gcc/4.3.2p3/${platform}/setup.sh
source ${contrib}/gcc/4.6/${platform}/setup.sh

#export PATH=${GCC_INSTALL_DIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${GCC_INSTALL_DIR}/lib64:${LD_LIBRARY_PATH}
#
#
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.05/${platform}/root
export LHAPDFPATH=/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/${platform}
export APPLGRID=/afs/cern.ch/user/c/cembree/public/progs
export FASTJET=/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${platform}
#
#
export PATH=${BASEDIR}/bin:${FASTJET}/bin:${APPLGRID}/bin::${ROOTSYS}/bin:$LHAPDFPATH/bin:${PATH}
export LD_LIBRARY_PATH=${BASEDIR}/lib:${FASTJET}/lib:${APPLGRID}/lib:${ROOTSYS}/lib:$LHAPDFPATH/lib:${LD_LIBRARY_PATH}
#
#export BASEDIR=$PWD
#export ARCH=`$APPLGRID/bin/rootarch.sh`
export  ARCH="-m64"
export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export CFLAGS=$ARCH
export FFLAGS=$ARCH
export LDFLAGS=$ARCH


###########################################
# link source code from other directories #
###########################################

#lhapdf sets
ln -s `lhapdf-config --pdfsets-path`  PDFsets

#R and B root files
#ln -s /afs/cern.ch/user/c/cembree/public/sherpatop21inputs .
ln -s /afs/cern.ch/work/c/cembree/public/sherpatop21inputs .

#other MyAnalysis files
export SRC=/afs/cern.ch/user/c/cembree/public/workspace/MyAnalysis

if (test -f MyCrossSection.cxx) then
 echo MyCrossSection file found
else
 ln -s  $SRC/MyCrossSection.cxx .
 ln -s  $SRC/MyCrossSection.h .
fi

if (test -f MyData.cxx) then
 echo MyData file found
else
 ln -s  $SRC/MyData.cxx .
 ln -s  $SRC/MyData.h .
fi

if (test -f MyFrame.cxx) then
 echo MyFrame file found
else
 ln -s  $SRC/MyFrame.cxx .
 ln -s  $SRC/MyFrame.h .
 ln -s  $SRC/MyFrameData.cxx .
 ln -s  $SRC/MyFrameData.h .
fi

if (test -f MyPDF.cxx) then
 echo MyPDF found
else
 ln -s  $SRC/MyPDF.cxx .
 ln -s  $SRC/MyPDF.h . 
fi

if (test -f MyBand.cxx) then
 echo MyPDF found
else
 ln -s  $SRC/MyBand.cxx .
 ln -s  $SRC/MyBand.h .
fi
