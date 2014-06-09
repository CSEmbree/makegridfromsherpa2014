STD = -std=c++0x
#STD = -std=c++0x -stdlib=libstdc++ #<--for Cam's Mac on OSX Mavricks 10.9

CXX = g++ $(STD)
#CXX = g++
F77 = gfortran

FFLAGS   += -O3 -fPIC 
CXXFLAGS += -O3 -fPIC -fpermissive 
LDFLAGS +=  -fno-automatic -fno-f2c -O3 -g # -I./src/Inc -Iobj

# root
ROOTINCS = $(shell root-config --cflags) 
ROOTLIBS = $(shell root-config --glibs) 
ROOTARCH = $(findstring -m64, $(ROOTINCS) )

#LHAPDF
LHAPDFINCS = -I$(shell lhapdf-config --prefix)/include
LHAPDFDIR  = $(shell lhapdf-config --libdir)
LHAPDFLIBS = -L$(LHAPDFDIR) -lLHAPDF

# applgrid
APPLCXXFLAGS = $(shell applgrid-config --cxxflags)
APPLCLIBS    = $(shell applgrid-config --ldcflags)
APPLFLIBS    = $(shell applgrid-config --ldflags)

# hoppet
HOPPETLIBS =  $(shell hoppet-config --libs)

# fastjet
FASTJET     = $(shell fastjet-config --prefix)
#FASTJETLIBS = -L$(FASTJET)/lib -lfastjet
FASTJETLIBS = $(shell fastjet-config --libs)
FASTJETINCS = $(shell fastjet-config --cxxflags)
CXXFLAGS   += $(FASTJETINCS)

# get the fotran runtime library for linking fortran 

FRTLLIB = $(shell gfortran $(CXXFLAGS) -print-file-name=libgfortran.a)
FRTLIB  = -L$(subst /libgfortran.a, ,$(FRTLLIB) ) -lgfortran

HOPPETINCS=$(shell hoppet-config --cxxflags)
HOPPETLIBS=$(shell hoppet-config --libs)

# now set up the compile and link flags and libs
#CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS) $(FASTJETINCS) $(HOPPETINCS)
CXXFLAGS += $(ROOTARCH) $(ROOTINCS) $(APPLCXXFLAGS) $(LHAPDFINCS) $(FASTJETINCS) $(HOPPETINCS)
#LDFLAGS  += $(ROOTARCH)
#FFLAGS   += $(ROOTARCH)

CXXFLAGS += -g
#FFLAGS += -g 

CLIBS += $(APPLCLIBS) $(HOPPETLIBS)  $(FASTJETLIBS) $(ROOTLIBS) $(FRTLIB) $(LHAPDFLIBS)
#CLIBS += $(APPLCLIBS) $(HOPPETLIBS) $(LHAPDFLIBS) $(FASTJETLIBS) $(ROOTLIBS) 
#CLIBS += $(ROOTLIBS) $(LHAPDFLIBS) $(HOPPETLIBS) $(APPLCLIBS) $(FASTJETLIBS) 
FLIBS += $(LHAPDFLIBS) $(HOPPETLIBS) $(APPLFLIBS) $(ROOTLIBS) $(FRTLIB)

#OFILE= MyFrame.o t3.o fjClustering.o MyData.o MyFrameData.o MyEvent.o MyGrid.o MySubProcess.o sherpa_pdf.o MyCrossSection.o

#<-- added MyBand.o for references
OFILE= AtlasStyle.o MyFrame.o t3.o fjClustering.o MyHisto.o MyData.o MyFrameData.o MyBand.o MyEvent.o MyObservable.o MyGrid.o MyPDF.o MyCrossSection.o 

install : all 
#all : makegridfromsherpa addgrid
all : makegridfromsherpa  plot

makegridfromsherpa: makegridfromsherpa.o  $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

testgridfromsherpa: testgridfromsherpa.o  $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

plot: plot.o  $(OFILE)
	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)

#addgrid: addgrid.o  $(OFILE)
#	$(CXX) $(LDFLAGS) -o $@ $(OFILE) $<  $(CLIBS)


.SUFFIXES : .cxx .o .f .c .C

.f.o :
	$(F77) $(FFLAGS)   -c $<

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 


clean:
	rm -rf ./.libs ./.obj *.lo *.o *.la  *~
