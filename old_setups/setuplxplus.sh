
#
#-----------------------
## Set an LCG platform tag if not already defined
#test -n "$LCGTAG" || LCGTAG=x86_64-slc5-gcc43-opt
#
export platform=x86_64-slc6-gcc46-opt

#source ../setup.myroot.sh

#export BASEDIR=/afs/cern.ch/user/t/tcarli/scratch1/applgrid
export BASEDIR=/afs/cern.ch/user/c/cembree/public/progs

#cd $BASEDIR
#ls -l setup_slc6.sh
#source setup_slc6.sh 
#cd -

source setup_slc6.sh

#export ARCH=`$BASEDIR/bin/rootarch.sh`
#export CXXFLAGS=$ARCH
#export F90FLAGS=$ARCH
#export CFLAGS=$ARCH
#export FFLAGS=$ARCH
#export LDFLAGS=$ARCH




#export FASTJET=/afs/cern.ch/sw/lcg/external/fastjet/3.0.3/${platform}
#export PATH=${FASTJET}/bin:${BASEDIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${FASTJET}/bin:${BASEDIR}/lib:${LD_LIBRARY_PATH}
export PATH=${BASEDIR}/bin:${PATH}
export LD_LIBRARY_PATH=${BASEDIR}/lib:${LD_LIBRARY_PATH}

#export PATH=${FASTJET}/bin:${BASEDIR}/bin:${PATH}
#export LD_LIBRARY_PATH=${FASTJET}/bin:${BASEDIR}/lib:${LD_LIBRARY_PATH}

# to redo the makefiles etc


