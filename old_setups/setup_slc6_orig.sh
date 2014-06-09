
#
#-----------------------
#source  /home/tcarli/setup.sh

#export BASEDIR=/afs/cern.ch/user/t/tcarli/scratch0/applgrid
#export BASEDIR=/home/admin/Desktop/applgrid-1.3.3
export BASEDIR=/home/tcarli/applgrid
#export BASEDIR=/home/tcarli/applgridnew
#export ARCH=`$BASEDIR/bin/rootarch.sh`
export CXXFLAGS=$ARCH
export F90FLAGS=$ARCH
export CFLAGS=$ARCH
export FFLAGS=$ARCH
export LDFLAGS=$ARCH

echo basedir is $BASEDIR

export MYDIR=$PWD

source /home/tcarli/setup.sh
cd  $BASEDIR
source setup.sh
cd  $MYDIR


#export PATH=/usr/local/bin:${PATH}

#export LHAPATH=`lhapdf-config --pdfsets-path` 
#export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}



#ln -s `lhapdf-config --pdfsets-path` PDFsets

#
# link source code from other directories
#
export SRC=/home/tcarli/MyAnalysis

if (test -d $SRC) then
 echo src exists $SRC
else
 echo check out $SRC 
 svn co svn+ssh://svn.cern.ch/reps/atlasusr/tcarli/MyAnalysis $SRC
fi

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

  





