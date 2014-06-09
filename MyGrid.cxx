//
//   for sherpa
//
#include <iostream>
using namespace std;
#include <string>
//#include "appl_grid/appl_igrid.h"

#include "MyGrid.h"
/******************************************************************
 ** Method Implementations
 ******************************************************************/
MyGrid::MyGrid(string name, string version)
{

    debug=false;
    bookrefsubprocess=true;

    phasespacerun=false;
    //alluncorrevents=  0; // uncorrelated events passing cuts (subevents are counted as one event)
    //allevents      =  0; // all events passing cuts
    uncorrevents.clear();
    events.clear();
    newevent=false;
    eventid=-99;

    // Grid architecture
    nXbins  = 30;
    xLow    = 1.0e-8;
    xUp     = 1.0;
    xorder  = 5;
    //xorder  = 0;

    nQ2bins = 2; 
    q2Low   = 30624.; // this should be automated for the moment mtop=175*175
    q2Up    = 30628.; //
    qorder  = 1;      // there seems to be a bug, if only one bins is used

    iorder  = 0;

    // lowest order in alphas
    lowest_order = 2;
    // how many loops
    nloops = 1;

    processnumber=1;

    pdf_function="";
    subprocesssteername="";
    ntupname="";
    steername=name;
    versionname=version; //optional argument

    fixedinputgrid=false;
    bool incljets=false;

    cout<<" MyGrid steering file name= "<<steername<<endl;

    if (!file_exists(name)) {
     cout<<" MyGrid steering file name not found= "<<name.c_str()<<endl;
     exit (1);
    }

    gridname.clear();
    this->ReadSteering(name);

    pdf_function=subprocesssteername; //subprocesssteername was read into from reading steering file

    this->Initialize();

    const int nSub=nsub;

    weight= new double[nSub]; // needs to be put after ReadSteering

    istylemarkerdata.clear();
    icolormarkerdata.clear();
    mydata.clear();
    //isBooked.clear();

    //save important vector of histograms into a new vector for later retreival when adding another to this one
    std::vector<vhref> referenceHistogramVectors;
    if (bookrefsubprocess) {
     referenceHistogramVectors.clear();
     referenceHistogramVectors.push_back(hreference);        //save pointer to hreference
     referenceHistogramVectors.push_back(hreferencefine);    //save pointer to hreferencefine
    }
    return;
}

void MyGrid::Initialize() {

    if (debug) cout<<" MyGrid::Initialize start"<<endl;
    if (debug) cout<<" MyGrid::Initialize Number of grids to produce "<<gridname.size()<<endl;

    cout << " MyGrid::Initialize subprocesssteername: " <<subprocesssteername <<  endl;

    // Loop over grid
    for (int  igrid = 0; igrid < gridname.size(); igrid++) {

        cout<<" MyGrid::Initialize Grid Name "<<this->GetGridFileName(igrid)<<endl;
        //this->SetStyleMarker (20+igrid);
        //this->SetColorMarker (1+igrid);
        this->SetStyleMarker (20);
        this->SetColorMarker (1);

        MyData *mydatatmp= new MyData;
        if (!mydatatmp) cout<<" MyGrid::Initialize MyData could not be created "<<this->GetGridFileName(igrid)<<endl;

        TString dataname=TString(this->GetDataName(igrid));
        dataname+=".txt";
        if (debug) cout<<" MyGrid::Initialize name= "<<TString(this->GetDataName(igrid)).Data()<<endl;
        else cout<<" MyGrid::Initialize MyData created name= "<<dataname.Data()<<endl;

        //if (debug)
        cout<<" MyGrid::Initialize name= "<<dataname.Data()<<endl; 
        cout<<" MyGrid::Initialize gridnamedir= "<<gridnamedir<<endl;
         
	//mydatatmp->ReadData(dataname.Data(),gridnamedir.c_str());
	//TString mydir=gridnamedir.c_str();
        string mydataname=string(dataname.Data());
	mydatatmp->ReadData(mydataname,gridnamedir);

        mydatatmp->SetMarker(20,1);
        mydata.push_back(mydatatmp);

        uncorrevents.push_back(0);
        events.push_back(0);

	if (debug) if (fixedinputgrid) cout<<" MyGrid::Initialize fixed input file is on ! "<<endl;

        this->book_grid(igrid);
        isBooked.push_back(true);

        if (bookrefsubprocess) {

            TString hnameref="href";
            //if (fixedinputgrid) hnameref+=this->GetGridVersionName();
            TH1D *htmp=this->GetHisto(igrid,hnameref.Data());
            MyHisto *myhisttmp= new MyHisto(htmp);
            hreference.push_back(myhisttmp);

            TString hnamereffine="hreffine"; //+this->GetGridVersionName();
            //if (fixedinputgrid) hnamereffine+=this->GetGridVersionName();
            htmp=this->GetHisto(igrid,hnamereffine.Data());
            if (!htmp) cout<<" MyGrid::Initialize Histo not found "<<hnamereffine.Data()<<endl;

            int nbins=htmp->GetNbinsX();
            double xmin=htmp->GetBinLowEdge(1);
            double xmax=htmp->GetBinWidth(nbins)/2.+htmp->GetBinCenter(nbins);
            nbins*=10;
            //cout<<" MyGrid::Initialize nbins= "<<nbins<<" xmin= "<<xmin<<" xmax= "<<xmax<<endl;
            TString htit=htmp->GetTitle()+TString(" fine bins");
            hnamereffine+=this->GetGridVersionName();
            TH1D *hfine= new TH1D(hnamereffine,htit,nbins,xmin,xmax);
            MyHisto *myhisttmpfine= new MyHisto(hfine);
            hreferencefine.push_back(myhisttmpfine);

            hrefLOsubprocesseshistos.clear();
            hrefsubprocesseshistos.clear();

            vhref vtmphreflo;
            vhref vtmphref;

            for ( int isub=0; isub<nsub; isub++) {
                TString hname="hrefLOsub";
                hname+=isub;
                //cout<<" hname= "<<hname<<endl;
                TH1D *htestlo=this->GetHisto(igrid,hname.Data());
                TString htit=TString(hname)+TString(" LO ")+TString(GetGridVersionName().c_str());//TString(this->GetSubProcessName(isub)); //<--**
                htestlo->SetTitle(htit);
                MyHisto *myhisttmplo= new MyHisto(htestlo);
                vtmphreflo.push_back(myhisttmplo);

                TString hname2="hrefsub";
                hname2+=isub;
                TH1D *hnlo=this->GetHisto(igrid,hname2.Data());
                TString htit2=TString(hname2)+TString(" ")+TString(GetGridVersionName().c_str());//TString(this->GetSubProcessName(isub)); //<--**
                hnlo->SetTitle(htit2);
                MyHisto *myhistsubnlo= new MyHisto(hnlo);
                vtmphref.push_back(myhistsubnlo);

            }
            hrefLOsubprocesseshistos.push_back(vtmphreflo);
            hrefsubprocesseshistos.push_back(vtmphref);
        }
    }


    if (debug) 
     cout<<" MyGrid::Initialize dynamic_cast using subprocesssteername: "<<subprocesssteername<<endl;
    
    /*
    if(!mypdf)
        mypdf = new basic_pdf();
    else
        mypdf=dynamic_cast<basic_pdf*>( appl::appl_pdf::getpdf("basic") );
    */

    //mypdf = dynamic_cast<basic_pdf*>( appl::appl_pdf::getpdf("basic") ); //TEST-basic
    //mypdf = dynamic_cast<generic_pdf*>( appl::appl_pdf::getpdf(subprocesssteername) ); //TEST-generic
    mypdf = dynamic_cast<lumi_pdf*>( appl::appl_pdf::getpdf(subprocesssteername) ); //TEST-lumi
    
    if(!mypdf)
        cout<<" MyGrid::Initialize Warning: mypdf not found"<<endl;

    //mypdf->PrintSubprocess(); //<--** was used to check that subprocesses had been received
    

    //nsub=mypdf->GetSubProcessNumber(); //TEST-basic,generic
    nsub=mypdf->size(); //TEST-lumi
    std::cout<<" MyGrid::Initialize nsub= "<<nsub<<std::endl;
    
    //nsub=mygrid.back()->subProcesses(0);//<--**appl_grid method. What was does parameter mean? The 0, 1, ...
    if(nsub<1)
        cout<<" MyGrid::Initialize WARNING: nsub: "<<nsub<<endl;

    //cout<<" finished MyGrid ini "<<endl;


    return;

}


void MyGrid::PrintRefHistos(int igrid) {
    if (!bookrefsubprocess) return;
    if(!hreference[igrid]) cout<<" MyGrid::PrintRefHistos reference histogram igrid= "<<igrid<<" not found "<<endl;
    hreference[igrid]->Print();
}

void MyGrid::PrintSubprocessRefHistos(int igrid) {
    if (!bookrefsubprocess)  return;
    cout<<" MyGrid::PrintSubprocessRefHistos reference histos "
        <<hrefsubprocesseshistos[igrid].size()<<endl;
    for (int i=0; i<nsub; i++) {
        if (!hrefsubprocesseshistos[igrid][i]) cout<<" reference histogram i= "<<i<<" not found "<<endl;
        else {
	  cout<<" MyGrid:PrintSubProcessRefHistos: Histo= "<<hrefsubprocesseshistos[igrid][i]->GetHisto()->GetTitle()<<endl;
            hrefsubprocesseshistos[igrid][i]->Print();
        }
    }

    return;
}

TH1D *MyGrid::GetReference(int igrid){
 if ((int)mygrid.size()<=igrid) {
  cout<<" MyGrid::GetReference: WARNING: mygrid vector too small igrid="<<igrid<<endl;
  exit (1);
 }
 if (!mygrid[igrid]) cout<<" MyGrid::GetReference mygrid["<<igrid<<"] not filled ! "<<endl;
  TH1D *href=(TH1D*)mygrid[igrid]->getReference();
  if (!href) {
   cout<<" MyGrid::GetReference: WARNING: reference histo not found igrid="<<igrid<<endl;
  }

  return href;
}


TH1D *MyGrid::TH1NormStatError(TH1D *hsum, TH1D *hsum2, double norm) {
    TH1D* h1new=(TH1D*) hsum->Clone(hsum->GetName());

    Double_t x, y, y2, ey;
    for (Int_t i=0; i<=hsum->GetNbinsX(); i++) {
        // sherpa formula
        //y= hsum ->GetBinContent(i)/norm;
        //y2=hsum2->GetBinContent(i)/norm;
        //ey=sqrt((y2-y*y)/(norm-1));

        y= hsum ->GetBinContent(i);
        y2=hsum2->GetBinContent(i);
        ey=sqrt((y2/norm-(y*y)/(norm*norm))/(norm-1)); //<--**think about


        y/=norm;
        ey/=norm;

        h1new->SetBinContent(i,y);
        h1new->SetBinError(i,ey);

        //if (debug) cout <<"bincenter= "<<hsum->GetBinCenter(i)<<" Binw = " << x << " y= " << y << endl;
    }
    return h1new;
}


void MyGrid::NormaliseRefHistos(int igrid, int nev) {
 if (!bookrefsubprocess)  return;

 if ((int)hreference.size()<igrid) {
  cout<<" MyGrid::NormaliseRefHistos: WARNING: hreference vector too short "<<endl;
  return;
 }


 //debug=true;
 double yfac=mydata[igrid]->GetUnitfbFactor();
 double xfac=mydata[igrid]->GetUnitGeVFactor();
 if (debug) cout<<" MyGrid::NormaliseRefHistos: normalise xfac= "<<xfac<<" yfac= "<<yfac<<endl;

 bool binw=mydata[igrid]->DivideByBinWidth();
 if (debug) cout<<" MyGrid::NormaliseRefHistos: divided by bin width  "<<endl;

 // bool normtot= false;
 bool normtot=mydata[igrid]->isNormTot();
 if (normtot&&debug) cout<<" MyGrid::NormaliseRefHistos: normalise to total cross section  "<<endl;

 if (nev==0) {
  cout<<" MyGrid::NormaliseRefHistos: WARNING: nev is zero  "<<endl;
  return;
 } else if (debug)
  cout<<" MyGrid::NormaliseRefHistos:  nev= "<<nev<<endl;

 if (!hreference[igrid]) {
  cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hreference' for igrid: "
      <<igrid<<" not found!"<<endl;
   exit(0);
 } else {
   this->Normalise(hreference[igrid]->GetHisto(),yfac/nev,xfac,binw,normtot);
 }

 if (!hreferencefine[igrid]) {
  cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hreferencefine' for igird: "
      <<igrid<<" not found!"<<endl;
  exit(0);
 } else {
   this->Normalise(hreferencefine[igrid]->GetHisto(),yfac/nev,xfac,binw,normtot);
 }

  for (int iproc=0; iproc<nsub; iproc++) {
   if (!hrefsubprocesseshistos[igrid][iproc]) {
    cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hrefsubprocesseshistos' for igrid: "
        <<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
    exit(0);
   } else {
     this->Normalise(hrefsubprocesseshistos[igrid][iproc]->GetHisto(),yfac/nev,xfac,binw,normtot);
   }

   if (!hrefLOsubprocesseshistos[igrid][iproc]) {
    cout<<" MyGrid::NormaliseInternalRefHistos: WARNING: Reference histogram 'hrefLOsubprocesseshistos' for igrid: "
        <<igrid<<", iproc: "<<iproc<<" not found!"<<endl;
     exit(0);
   } else {
    this->Normalise(hrefLOsubprocesseshistos[igrid][iproc]->GetHisto(),yfac/nev,xfac,binw,normtot);
   } 
  }
 }


//
void MyGrid::ReadSteering(string fname) {

    steername=fname;
    if (debug) cout<<" MyGrid::ReadSteering steering= "<<steername<<endl;

    ifstream infile(steername.c_str(), ios::in);
    if(!infile) { // Check open
        cerr << " MyGrid::ReadSteering: Can't open " << steername <<"\n";
        infile.close();
        exit (1);
    } else {
        cout <<" MyGrid::ReadSteering read data file: " << steername << endl;
    }

// cout<<" filed opened "<<fname<<endl;

    int iline=0;
    int nsyst=1;

    char line[256];
    char lineFirstWord[256]; //same length as 'line' string buffer
    char text[100];
    char name[100];
    float mys;

    while (infile.good()) {
        if (debug) cout << " good: " << infile.good() << " eof: " << infile.eof() << endl;
        //if (!infile.good()) break;
        infile.getline(line,sizeof(line),'\n');
        if (debug) cout<< "line= "<< line << "\n";

        sscanf(line," %s",lineFirstWord);
        std::string cpp_lineFirstWord(lineFirstWord);

        if(line[0] != '%') {
            /*if(line[0] != 'g' && line[0] != 'n' && line[0] != 's' && line[0] != 'p' &&
               line[0] != 'q'
               ){*/
            //  if (mygridoptions->isKnownOption(cpp_lineFirstWord)==false) {
                // do something
                //} else if (strstr(line,"nsub")!=0) {
                //char text[100];
                //sscanf(line," %s %d ",text, &nsub);
                //if (debug) printf("**MyGrid Read: Nsub= %d   \n",nsub);
           if (strstr(line,"nprocessnumber")!=0) {
                sscanf(line," %s %d ",text, &processnumber);
                if (debug) printf("**MyGrid::Read processnumber= %d   \n",processnumber);
           } else if (strstr(line,"subprocesssteername")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                subprocesssteername=string(name);
                if (debug) printf("**MyGrid::Read subprocesssteername= %s \n",subprocesssteername.c_str());
                //cout<<"     subprocesssteername= "<<subprocesssteername<<endl;
                //   } else if (strstr(line,"pdffunction")!=0) {
                //sscanf(line," %s %[^\n] ",text, name);
                //pdf_function=string(name);
                //if (debug) printf("**MyGrid Read: pdf_function= %s \n",pdf_function.c_str());
                //cout<<"   pdffunction= "<<pdf_function<<endl;

           } else if (strstr(line,"ntupdiroutput")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                ntupdiroutput=name;
           } else if (strstr(line,"fixedinputgridname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                fixedinputgrid=true;
                if (fixedinputgrid) cout<<" MyGrid::Read fixedinputgrid is on "<<endl;
                fixedinputgridname=name;
           } else if (strstr(line,"ntupdirinput")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                ntupdirinput=name;
           } else if (strstr(line,"ntupname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                cout<<" MyGrid::Read "<<text<<" "<<name<<endl;
                ntupname=name;
           } else if (strstr(line,"nXbins")!=0) {
                sscanf(line," %s %f ",text, &mys);
                cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                nXbins=mys;
           } else if (strstr(line,"xLow")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                xLow=mys;
           } else if (strstr(line,"xUp")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                xUp=mys;
           } else if (strstr(line,"xorder")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<mys<<endl;
                xorder=mys;
           } else if (strstr(line,"nQ2bins")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                nQ2bins=mys;
           } else if (strstr(line,"q2Low")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                q2Low=mys;
           } else if (strstr(line,"q2Up")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                q2Up=mys;
           } else if (strstr(line,"qorder")!=0) {
                sscanf(line," %s %f ",text, &mys);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<mys<<endl;
                qorder=mys;
           } else if (strstr(line,"griddir")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<"MyGrid::Read "<<text<<" "<<name<<endl;
                gridnamedir=name;
           } else if (strstr(line,"gridname")!=0) {
                sscanf(line," %s %[^\n] ",text, name);
                if (debug) cout<<" MyGrid::Read "<<text<<" "<<name<<endl;
                string myname=name;
                gridname.push_back(myname);
           };
        };
        memset(text, '\0', 100); //reset option input buffers for next run
        memset(name, '\0', 100);
    };
    return;
};

void MyGrid::Print() {

    cout<<" "<<endl;
    cout<<"\n MyGrid >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
    cout<<" steering: "<<steername<<endl;
    cout<<" ntupdirinput= "<<ntupdirinput<<" ntupname=  "<<ntupname<<endl;
    cout<<" ntupdiroutput= "<<ntupdiroutput<<endl;
    cout<<" gridnamedir=  "<<gridnamedir<<endl;
    cout<<" subprocesssteername: "<<subprocesssteername<<endl;
    //cout<<" pdf_function=  "<<pdf_function<<endl;

    cout<<"\n Grid architecture: "<<endl;
    cout<<"\t  nXbins= "<<nXbins<<"  xLow= "<<xLow<<" xUp= "<<xUp<<" xorder= "<<xorder<<endl;
    cout<<"\t nQ2bins= "<<nQ2bins<<" q2Low= "<<q2Low<<" q2Up= "<<q2Up<<" qorder= "<<qorder<<endl;

    cout<<"\t lowest_order= "<<lowest_order<<" nloops= "<<nloops<<endl;

    //cout<<"\n Selected Events: "<<endl;
    //cout<<"all uncorrellated events= "<<alluncorrevents
    //     <<" all events= "<<allevents<<endl;

    cout<<"\n Number of grids N= "<<gridname.size()<<endl;
    for (int  i = 0; i <   gridname.size(); i++)
    {
        cout<<" "<<gridname[i]
            <<" uncorrelated events= "<<uncorrevents[i]
            <<" total events= "<<events[i]
            <<" style: "<<this->GetStyleMarker(i)
            <<" color: "<<this->GetColorMarker(i)
            <<endl;
        //mydata[i]->Print();
    }


    cout<<" MyGrid <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<endl;

}

bool MyGrid::file_exists(const string& s) {

    if ( FILE* testfile=fopen(s.c_str(),"r") ) {
     fclose(testfile);
     //cout<<" file ok "<<s<<endl;
     return true;
    }
    else {
      //cout<<" file not ok "<<s<<endl;
     return false;
    }
}

void MyGrid::book_grid(int igrid)  // inital grid booking
{

    time_t _t;
    time(&_t);

    cout<<" ***********************************************"<<endl;
    cout<<" MyGrid::book_grid("<<igrid<<")"<<endl;
    cout<<" booking the grids " << ctime(&_t) <<endl;

    //double apramval=5.;
    //appl::igrid::transformvar(apramval);
    //cout<<" MyGrid::book_grid "<<glabel<<endl;

    bool create_new = false;
    // if the file does not exist, create a new grid...
    //
    /*
    string filename=this->GetGridFullFileName(igrid);//+this->GetGridVersionName();
    //string filename=this->GetGridName(igrid);//+this->GetGridVersionName();

    string _filename = filename;
    string newFileName;
    if (debug) cout<<" MyGrid::book_grid before fixed filename= "<<filename.c_str()<<endl;
    */
   /*
    if (debug&&fixedinputgrid) cout<<" MyGrid::book_grid fixedinputgrid is ON "<<endl;

    if (debug) cout<<" MyGrid::book_grid fixedinputgridname= "<< fixedinputgridname.c_str()<<endl;
    bool emptyfixedinput=TString(fixedinputgridname).IsWhitespace();
    if (debug&&fixedinputgrid) cout<<" MyGrid::book_grid fixedinputgridname is empty "<<endl;

    if (fixedinputgrid) 
      if (emptyfixedinput) newFileName+= ".root";
      else                 newFileName=TString(fixedinputgridname)+".root";
    else
   */
/*
    newFileName=this->GetGridVersionName()+".root";
    _filename.replace( _filename.find(".root"), 5, newFileName.c_str());
    filename=_filename;
*/
    TString vers=TString(this->GetGridVersionName().c_str());
    cout<<" MyGrid book vers= "<<vers.Data()<<endl;

    string filename=this->GetGridFullFileName(igrid);
    if (vers.Contains("-B") ||vers.Contains("-R") ) {
     if (debug) cout<<" filename= "<<filename.c_str()<<endl;
     TString newname=this->GetGridFullFileName(igrid);
     newname.ReplaceAll(".root","");
     newname+=vers+".root";
     if (debug) cout<<" newname= "<<newname.Data()<<endl;
     string command;
     command+="cp ";
     command+=filename;
     command+=" ";
     command+=newname;
     cout<<" MyGrid executing "<<command.c_str()<<endl;
     int i=system (command.c_str());
     filename=newname;
    }

    if (debug) cout<<" MyGrid::book_grid filename= "<<filename.c_str()<<endl;

    /*
    if (debug) if (fixedinputgrid) cout<<" MyGrid::book_grid fixedinputgrid true "<<endl;
    if (fixedinputgrid) {
     filename=fixedinputgridname;
     cout<<" MyGrid::book_grid working with fixed grid name ..create_new overwritten to false "<<endl;  
     create_new = false;
    }
    */
    if (debug) cout<<" MyGrid::book_grid igrid= "<<igrid<<" checking if this file exists: "<<filename<<endl;

    if (!file_exists(filename)) {
     cout<<" MyGrid::book_grid file named '"<<filename<<"' does not exist. set create_new=true"<<endl;
     create_new = true;
     phasespacerun=true;
    }
   /*
    // or if it does exists but root file is a zombie...
    if (!create_new) {
     cout<<" MyGrid::book_grid create test file "<<endl;
     TFile testFile( (TString(filename)));
     if ( testFile.IsZombie() ) {
      cout<<" MyGrid::book_grid test file is a zombie! Create_new=true"<<endl;
      create_new = true;
     }
     testFile.Close();
    }
   */
    //if (debug) 
     cout<<" MyGrid::book_grid final create_new is: "<<create_new<<" for file: "<<filename<<endl;
     if (debug) if (phasespacerun) cout<<" MyGrid::book_grid This is a phasespace run ! "<<endl;

    //Print();

    //mydata[igrid]->Print();
    int nObsBins=mydata[igrid]->GetNBins();
    double *obsBins=mydata[igrid]->GetBins();

    //  for (int i=0; i<nObsBins+1; i++)
    // cout<<i<<" obsBins= "<<obsBins[i]<<endl;

    if (create_new) {
        cout<<" >>>>---------------------------------------------- "<<endl;
        cout<<" MyGrid::book Creating NEW grid... " << filename << "\tigrid " << igrid << endl;
        cout<<" nObsBins= "<<nObsBins<<endl;
        cout<<" nQ2bins= "<<nQ2bins<<" q2Low= "<<q2Low<<" q2Up= "<<q2Up<<" qorder= "<<qorder<<endl;
        cout<<" nXbins= "<<nXbins<<" xLow= "<<xLow<<" xUp= "<<xUp<<" xorder= "<<xorder<<endl;
        cout<<" pdf_function= "<<pdf_function<<endl;
        cout<<" lowest_order= "<<lowest_order<<" nloops= "<< nloops<<endl;
        cout<<" <<<<---------------------------------------------- "<<endl;
 
        /*
        appl::grid *tmpgrid = new appl::grid( nObsBins, obsBins,      // obs bins
                                              nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
                                              nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
                                              pdf_function, lowest_order, nloops ); //TEST-generic
        */
                                              
        /*
        appl::grid *tmpgrid = new appl::grid( nObsBins, obsBins,      // obs bins
                                              nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
                                              nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
                                              "basic", lowest_order, nloops );  //TEST-basic
        */
        
        
        
        ///******RUN WITH STEERING FILE: atlas2012_top-config.txt that will make pdf_function=top.config
        if (debug) printf(" MyGrid::fill pdf_function= %s \n",pdf_function.c_str());

        cout<<" MyGrid::book_grid pdf_function= "<<pdf_function<<endl;

        appl::grid *tmpgrid = new appl::grid( nObsBins, obsBins,      // obs bins
                                              nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
                                              nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
                                              pdf_function, lowest_order, nloops );  //TEST-lumi
        
        //exit(-1); //TEST
        if (debug) cout<<" MyGrid::book_grid tmpgrid booked "<<endl;
        
        mygrid.push_back(tmpgrid);

        //cout<<"Mygrid::book_grid printing appl_grid..."<<endl;
        //tmpgrid->print();

        //double sqrts=7000.; //put here MyData
        //mygrid[igrid]->setCMScale(sqrts ); // done later

        //int nsub=mypdf->GetSubProcessNumber(); //<--**
        ///int nsub=tmpgrid->subProcesses(0); //<--**appl_grid method. What was does parameter mean? The 0, 1, ...
        //int nsub=mygrid.back()->subProcesses(0);//<--**appl_grid method. What was does parameter mean? The 0, 1, ...

        //cout << " MyGrid::book_grid nsub = " << nsub << endl;
        cout << " MyGrid::book_grid reference histo name = " << mygrid[igrid]->getReference()->GetName() << endl;

    } else {
        cout <<" MyGrid::book_grid Using existing grid file "<<filename<<endl;
        appl::grid *tmpgrid = new appl::grid(filename); //optimise grid x,Q2 bins
        mygrid.push_back(tmpgrid);
        mygrid[igrid]->getReference()->Reset();
        mygrid[igrid]->optimise(nQ2bins, nXbins);
        cout<<*(mygrid[igrid])<<endl;
    }


    //appl_pdf::printmap();
    //appl_pdf* mypdf = appl_pdf::getpdf("topsubprocesses.dat"); //<--**


    //std::cout << "grid name " << mypdf << std::endl; //<--**
    //std::cout << "grid name " << mypdf->name() << std::endl; //<--**

    int subprocessesi=0;

    nsub=mygrid[igrid]->subProcesses(subprocessesi); //<--
    cout<<" MyGrid::book_grid subprocesses= "<<mygrid[igrid]->subProcesses(subprocessesi)<<endl;//<--

    cout<<" ***********************************************"<<endl;

}



void  MyGrid::fill(MyEvent *myevent, double obs )
{

 // counter of number of events

    //this->SetEventId(myevent->GetEventId());
    bool countevent=myevent->GetCountEvent(); 

    bool newevent=myevent->NewEvent();
    if (countevent) myevent->SetEventId(myevent->GetEventId()); //turn off newevent


    if (debug) {
     cout<<" MyGrid::fill After event ID set, newevent= "<<(newevent?"TRUE":"FALSE")
                                        <<" countevent= "<<(countevent?"TRUE":"FALSE")
         <<" events= "<<events[0]<<endl;
    }

    double xsec=myevent->GetXSection();
    double x1=myevent->GetX1();
    double x2=myevent->GetX2();
    double q2=myevent->GetQ2();
    double mewgt=myevent->GetWeight();

    if (debug)
     cout<<" MyGrid::fill x1= "<<x1<<" x2= "<<x2<<" q2= "<<q2
         <<" mewgt= "<<mewgt<<" xsec= "<<xsec<<endl;

    if (x1<xLow) cout<<" MyGrid::fill x1 outside range x1= "<<x1<<" xLow= "<<xLow<<endl;
    if (x2<xLow) cout<<" MyGrid::fill x2 outside range x2= "<<x1<<" xLow= "<<xLow<<endl;
    if (x1>xUp)  cout<<" MyGrid::fill x1 outside range x1= "<<x1<<" xUp= "<<xUp<<endl;
    if (x2>xUp)  cout<<" MyGrid::fill x2 outside range x2= "<<x1<<" xUp= "<<xUp<<endl;
    if (q2<q2Low)cout<<" MyGrid::fill q2 outside range q2= "<<q2<<" q2Low= "<<q2Low<<endl;
    if (q2>q2Up) cout<<" MyGrid::fill q2 outside range q2= "<<q2<<" q2Up= "<<q2Up<<endl;

    int Ngrids=this->GetNGrid();
    if (debug) cout<<" MyGrid::fill Ngrids= "<<Ngrids<<endl;
    int id1=myevent->GetID(0);
    int id2=myevent->GetID(1);
    if (debug) cout<<" MyGrid::fill incoming partons id1= "<<id1<<" id2= "<<id2<<endl;

//
// Need to prepare weight vector in subprocesses
//

    int nproc=this->GetProcessNumber();
    if (debug) cout<<" MyGrid::fill nproc= "<<nproc<<endl;
    

    if (nproc==1 || nproc==2) { // process 1 and 2 is W
        int icharge=myevent->GetLeptonCharge();
        if (icharge== 1) nproc=1;
        if (icharge==-1) nproc=2;
        if (debug) cout<<" MyGrid::fill nproc= "<<nproc<<" icharge= "<<icharge<<endl;
    }

    //if (!mypdf) cout<<" MyGrid::fill mypdf not initialized "<<endl; //<--**
    //mypdf->SetCurrentProcess(nproc); //<--** Currently does nothing? nothing seen last code checked: 19-Jun-2013

    //
    // to be replaced by iorder ?
    // decided on basis of lowest_order insteed of nproc
    int myorder=-1;
    if (nproc==1||nproc==2) { // Wminus and Wplus
        myorder=myevent->GetOrder();
        if (myorder==1) myorder=0;
        if (myorder==2) myorder=1;
        if (debug) cout<<" MyGrid::fill myorder= "<<myorder<<endl;
    }

    if (nproc==3 || nproc==4) { // 3: inclusive jet 4: top
        myorder=myevent->GetOrder();
        if (myorder==2) myorder=0;
        if (myorder==3) myorder=1;
        if (debug) cout<<" MyGrid::fill myorder= "<<myorder<<endl;
    }
        

    if (myorder==-1) {
        cout<<" MyGrid::fill do not know what to do myorder= -1"<<endl;
        return;
    }

    //cout<<" MyGrid::fill Dynamically making generic_pdf for function 'decideSubProcess'"<<endl; //<--** test output
    //generic_pdf* mypdf  = dynamic_cast<generic_pdf*>( appl::appl_pdf::getpdf(pdf_function) ); //<--**

    //cout<<"MyGrid::fill decideSubProcess using: id1: "<<id1<<", id2: "<<id2<<endl;
    int iproc=mypdf->decideSubProcess(id1,id2); //<--**
    if (debug) cout<<" MyGrid::fill decideSubProc: "<<iproc<<endl;
    if (iproc==-1)  {
        cout<<" MyGrid::fill do not know what to do "<<endl;
        cout<<" MyGrid::fill incoming partons where id1= "<<id1<<" id2= "<<id2<<endl;
        myevent->Print();
        //exit(0); //TEST
        return;
    }

    if (debug) cout<<" MyGrid::fill valid iproc: "<<iproc<<endl;

    this->ResetWeight();
    
    this->SetWeight(iproc,mewgt);
    if (debug)  
     this->PrintWeight();

    double *weight=this->GetWeight();
    
    for (int igrid = 0; igrid < Ngrids; igrid++) {

        //
        // get cuts
        //
        if (debug) cout<<" MyGrid::fill filling igrid= "<<igrid<<endl;
        if (newevent&&countevent) {
          uncorrevents[igrid]++;
          if (debug) 
           cout<<" MyGrid::fill new event encountered: updating uncorrevents= "<<uncorrevents[igrid]<<endl;
        }
        if (countevent) events[igrid]++;
	//
	if (mewgt==0.) continue;
	//

        if (debug) cout<<" MyGrid::fill events= "<<events[igrid]<<" iproc= "<<iproc<<endl;


	if (debug) {
         cout<<" MyGrid::fill uncorrelated events= "<<uncorrevents[igrid]<<endl;
         cout<<" MyGrid::fill events= "<<events[igrid]<<" iproc= "<<iproc<<endl;
         this->PrintWeight();
         cout<<" MyGrid::fill Filling with: x1= "<<x1<<", x2= "<<x2<<", q2= "<<q2
             <<", mewgt= "<<mewgt<<", xsec= "<<xsec<<", myorder= "<<myorder<<endl;
        }

        //double obsmax=this->GetObsMax(igrid);
        //cout<<" obsmax= "<<obsmax<<endl;
        //if (obs>obsmax) {
        //cout<<" cut by obsmax bins= "<<obsmax<<" obs= "<<obs<<endl;
        //return; //speeds up code
        //}

      /*
        incljets=mydata[igrid]->isObsInclJets();
        if (debug) {
            if (incljets) cout<<" MyGrid::fill inclusive jet observable "<<endl;
            else          cout<<" MyGrid::fill inclusive jet not set "<<endl;
        }

        double ptmin =mydata[igrid]->GetJetPtCut();
        double rapmin=mydata[igrid]->GetJetMinY();
        double rapmax=mydata[igrid]->GetJetMaxY();
        if (debug) cout<<" MyGrid::fill ptmin= "<<ptmin<<" rapmin= "<<rapmin<<" rapmax= "<<rapmax<<endl;

        std::vector<fastjet::PseudoJet> myseljets=myevent->GetSelectedJets(ptmin,rapmin,rapmax);
        int nseljets=1;
        // if Observable is call inclusive jets loop over jets, otherwise fill grid only once
        double obs=0;
        if (incljets) {
            nseljets=myevent->GetNumberSelJets();
            if (debug) {
                cout<<" MyGrid::fill incljets number of selected jets= "<<nseljets<<endl;
                myevent->PrintSelJets2();
            }
        } else {
            if (debug) {
                cout<<" MyGrid::fill calling >PrintSelJets2 "<<endl;
                cout<<" MyGrid::fill number of selected jets= "<<myevent->GetNumberSelJets()<<endl;
                myevent->PrintSelJets2();
            }
        }

        //checking values are set correctly
        if (debug) myevent->Print(); //TEST

        for (int ijets=0; ijets<nseljets; ijets++) {
            if (incljets) obs=myseljets[ijets].pt();
            else if (nproc!=4)   obs=myevent->GetLeadingJet();
            else { // need to think about better implementation
                if (nproc==4) {
                    if (myevent->GetNumberSelJets()<2)
                        cout<<" MyGrid::fill something is wrong nsel= "
                            <<myevent->GetNumberSelJets()<<endl;
                    else {
		      if (debug) cout<<" MyGrid::fill calling GetInvariantMass12 "<<endl;
                        if (debug) cout<<" MyGrid::fill number of selected jets= "<<myevent->GetNumberSelJets()<<endl;
                        obs=myevent->GetInvariantMassTop12();
                    }
                    if (debug) cout<<" MyGrid::fill obs= "<<obs<<endl;
                }
            }
            if (debug) cout<<" MyGrid::fill "<<" ijets= "<<ijets<<" obs= "<<obs<<endl;
      */
            if (debug) cout<<" MyGrid::fill obs= "<<obs<<endl;

            if (mygrid[igrid]->isOptimised()) {
                if (debug) cout<<" MyGrid::fill filling grid and reference!\n"<<endl;
                if (!mygrid[igrid]) cout<<"MyGrid::fill grid not found igrid= "<<igrid<<endl;
                mygrid[igrid]->fill(x1, x2, q2, obs, weight, myorder );
            } else {
                if (debug) cout<<" MyGrid::fill filling phasespace!\n"<<endl;
                if (!mygrid[igrid]) cout<<" MyGrid::fill grid not found igrid= "<<igrid<<endl;
                mygrid[igrid]->fill_phasespace(x1, x2, q2, obs, weight, myorder );
            }

	    if (debug)
                cout<<" MyGrid::fill Filled with: x1= "<<x1<<" x2= "<<x2<<" obs= "<<obs<<" q2= "<<q2
                <<" mewgt= "<<mewgt<<" xsec= "<<xsec<<" myorder= "<<myorder<<endl;

            mygrid[igrid]->getReference()->Fill( obs,xsec);

 
            if (debug) 
             cout<<" MyGrid::fill filled refHisto with xsec: "<<xsec
                 <<", filled grid with (weight,iproc): ("<<mewgt<<","<<iproc<<")"<<endl;

            if (debug) {
	     cout<<" MyFill print reference events= "<<events[igrid]<<endl;
             mygrid[igrid]->getReference()->Print("all");
            }

            if (bookrefsubprocess) {
             if (myorder==0) {
	       hrefLOsubprocesseshistos[igrid][iproc]->Fill(obs,xsec,newevent);
             }

             hrefsubprocesseshistos[igrid][iproc]->Fill(obs,xsec,newevent);
             hreference[igrid]->Fill(obs,xsec,newevent);
             hreferencefine[igrid]->Fill(obs,xsec,newevent);

             //cout<<" MyGrid::fill subProcs for igrid: "<<igrid<<", obs: "<<obs<<", weight[iproc:"<<iproc<<"]: "<<weight[iproc]<<endl;


            }
    } //loop over grid
    return;
}

void MyGrid::Normalise(TH1D* h1, double yscale, double xscale=1., bool binw=false, bool normtot=false)
{

    double x, y, ey;
    double sigtot=0.;
   
    if (!h1) cout<<" MyGrid::Normalise histogram not found !"<<endl;
    double nev=h1->GetEntries();

    //this->CorrectStatUncertainty(h1);

    /*
    if (debug) {
     cout<<" MyGrid::Normalise histogram  ! before "<<endl;
     h1->Print("all");
    }
    */ 

    for (Int_t i=0; i<=h1->GetNbinsX()+1; i++) {
        y=h1->GetBinContent(i)*yscale;
        x=h1->GetBinWidth(i);
        //x=1.;
        sigtot+=y*x;
    }

//  cout<<" MyGrid::Normalise sigtot= "<<sigtot<<endl;

    for (Int_t i=0; i<=h1->GetNbinsX()+1; i++) {
        x =h1->GetBinWidth(i);
        y =h1->GetBinContent(i)*yscale*x;
        ey=h1->GetBinError(i)  *yscale*x;
        x =h1->GetBinWidth(i)  *xscale;
        if (x!=0) h1->SetBinContent(i,y/x);
        else      h1->SetBinContent(i,0.);

        if (x!=0) h1->SetBinError(i,ey/x);
        else      h1->SetBinError(i,0.);
    }
    if (normtot) {
        if (sigtot!=0.) 
            h1->Scale(1. / sigtot);
    }

    h1->SetEntries(nev);

    if (debug) {
     cout<<" MyGrid::Normalise histogram  ! after "<<endl;
     h1->Print("all");
    }

    return;
}

//
// just normalise to bin width
//
void MyGrid::DivideByBinWidth(TH1D* h)
{

    if (!h) cout<<" MyGrid::DivideByBinWidth histogram not found !"<<endl;
    if (debug) cout<<" MyGrid::DivideByBinWidth called "<<h->GetName()<<endl;
    
    double nev=h->GetEntries();
    for ( int ibin=0 ; ibin<=h->GetNbinsX() ; ibin++ )
    {
     double width = h->GetBinWidth(ibin);
     if (debug) cout<<" MyGrid::DivideByBinWidth width= "<<width
                    <<" bin content= "<<h->GetBinContent(ibin)
                    <<endl;
     if (width==0) cout<<" MyGrid::DivideByBinWidth something is wrong width= "<<width<<endl;
     h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
     h->SetBinError  ( ibin, h->GetBinError(ibin)/width );
    }

    h->SetEntries(nev);
    return;
}

//extern "C" void write_grid_(double& xstotal)   // writes out grid after some events
void MyGrid::write_grid()   // writes out grid after some events
{

 int Ngrids=this->GetNGrid();

 //if (debug) 
 cout<<" MyGrid::write_grid Write out the grid ... Ngrids " << Ngrids << endl;

 for (int igrid = 0; igrid < Ngrids; igrid++)
 {

  //mygrid[igrid]->getReference()->Print("all");
  if (debug) cout<<" MyGrid::write_grid saving grid N=" << igrid+1 << "\tof total " << Ngrids << endl;
  //mygrid[igrid]->setNormalised( true );

  if (debug) cout<<" MyGrid:write_grid uncorrevents= "<<uncorrevents[igrid]<<endl;

  int nev=uncorrevents[igrid]; // for weight2
  //int nev=events[igrid]; // for weight

  if (nev==0) {
   cout<<" MyGrid::write_grid run igrid= "<<igrid<<" no event found ! "<<endl;
   exit (1);
  }

  bool binw=mydata[igrid]->DivideByBinWidth();
  TH1D *href=mygrid[igrid]->getReference();
  this->Normalise(href,1./nev,1.,binw,false);

  cout<<" MyGrid::write_grid run igrid= "<<igrid<<" nev= "<<nev<<endl;
  mygrid[igrid]->run() = nev;

  mygrid[igrid]->trim();

  if (debug) {
   int trim_size = mygrid[igrid]->size();
   cout<<" MyGrid::write_grid trimmed "<<endl;
   mygrid[igrid]->untrim();
   int untrim_size = mygrid[igrid]->size();
   cout<<" MyGrid::write_grid Saved Space ratio: "<<(trim_size/untrim_size*1.)<<endl;
  }

  string filename=this->GetGridFullFileName(igrid);//+this->GetGridVersionName();
  string _filename = filename;

  string newFileName;
  // if (phasespacerun) {
  // cout<<" MyGrid::write_grid Writing to phasespacerun ON "  << endl;

  if (fixedinputgrid) { 
   bool emptyfixedinput=TString(fixedinputgridname).IsWhitespace();
   if (emptyfixedinput) newFileName+= ".root";
   else                 newFileName=TString(fixedinputgridname)+".root";
  } else
   newFileName=this->GetGridVersionName()+".root";

  _filename.replace( _filename.find(".root"), 5, newFileName.c_str());
  filename=_filename;

  //if (debug) 
  cout<<" MyGrid::write_grid Writing to filename= " << filename << "\tigrid: " << igrid << endl;
  if (mygrid[igrid])   mygrid[igrid]->Write(filename);
  else cout<<" MyGrid::write_grid mygrid"<<igrid<<"] not found ! "<<endl;

  if (bookrefsubprocess) {
   /*
   _filename = filename;
   newFileName = "-ghistos.root";
   _filename.replace( _filename.find(".root"), 5, newFileName.c_str());
  
   TFile *f = new TFile( _filename.c_str(),"recreate");
   //if (debug) f->Print();
   */

   TFile rootfile(filename.c_str(),"update");

   if (debug) {
    cout<<" MyGrid:Write print applgrid reference "<<endl;
    mygrid[igrid]->getReference()->Print("all");
   }
    
   this->NormaliseRefHistos(igrid,nev);
   //hreference[igrid]->GetHisto()->Print("all");
   hreference[igrid]->GetHisto()->Write();
   hreferencefine[igrid]->GetHisto()->Write();

   for (int iproc=0; iproc<nsub; iproc++) {
    hrefLOsubprocesseshistos[igrid][iproc]->GetHisto()->Write();
    hrefsubprocesseshistos[igrid][iproc]->GetHisto()->Write();
   }
  }
  //f->Close();
 }
 

 time_t _t;
 time(&_t);

  //cout<<" ***********************************************"<<endl;
  //cout<<" saved grids " << ctime(&_t) << endl;
  //cout<<" ***************************************************"<<endl;
 return;
}


//
// ----------------------------------------------
//    analysis
// ----------------------------------------------
//

/*
void MyGrid::AddGrid(MyGrid *myOtherGrid) {
    cout<<"MyGrid::AddGrid: called..."<<endl;


    //get all of the other grids reference histograms
    std::vector<std::vector<TH1D*> > *otherRefLOsubprocessesHistos;
    otherRefLOsubprocessesHistos = myOtherGrid->GetRefLOsubprocessesHistos();

    std::vector<std::vector<TH1D*> > *otherRefsubprocessesHistos;
    otherRefsubprocessesHistos = myOtherGrid->GetRefsubprocessesHistos();

    std::vector<TH1D*> *otherReferenceHistos;
    otherReferenceHistos = myOtherGrid->GetReferenceHistos();

    std::vector<TH1D*> *otherReferenceFineHistos;
    otherReferenceFineHistos = myOtherGrid->GetReferenceFineHistos();


    //make sure all histos can be added together
    if (hreference.size()!=otherReferenceHistos->size())
     cout<<" MyGrid::AddGrid: WARNING: hreference are the not the same but will attemp to be added!"<<endl;

    if (hreferencefine.size()!=otherReferenceFineHistos->size())
     cout<<" MyGrid::AddGrid: WARNING: hreferencefine are the not the same but will attemp to be added!!"<<endl;

    if (hrefLOsubprocesseshistos.size()!=otherRefLOsubprocessesHistos->size())
     cout<<" MyGrid::AddGrid: WARNING: hrefLOsubprocesseshistos are the not the same but will attemp to be added!"<<endl;

    if (hrefsubprocesseshistos.size()!=otherRefsubprocessesHistos->size())
     cout<<" MyGrid::AddGrid: WARNING: hrefsubprocesseshistos are the not the same but will attemp to be added!"<<endl;


    //add this grid's and the other grid's reference histos
    for(int href_x=0; href_x<hreference.size(); href_x++) {
        hreference.at(href_x)->Add(otherReferenceHistos->at(href_x));
    }

    for(int hreffine_x=0; hreffine_x<hreferencefine.size(); hreffine_x++) {
        hreferencefine.at(hreffine_x)->Add(otherReferenceFineHistos->at(hreffine_x));
    }

    for(int hrefLO_x=0; hrefLO_x<hrefLOsubprocesseshistos.size(); hrefLO_x++) {
        for(int hrefLO_y=0; hrefLO_y<hrefLOsubprocesseshistos.at(hrefLO_x).size(); hrefLO_y++) {
            hrefLOsubprocesseshistos.at(hrefLO_x).at(hrefLO_y)->Add(otherRefLOsubprocessesHistos->at(hrefLO_x).at(hrefLO_y));
        }
    }

    for(int hrefsub_x=0; hrefsub_x<hrefsubprocesseshistos.size(); hrefsub_x++) {
        for(int hrefsub_y=0; hrefsub_y<hrefsubprocesseshistos.at(hrefsub_x).size(); hrefsub_y++) {
            hrefsubprocesseshistos.at(hrefsub_x).at(hrefsub_y)->Add(otherRefsubprocessesHistos->at(hrefsub_x).at(hrefsub_y));
        }
    }



    //add event counters
    //cout<<"TEST: THISalluncorrevents: "<<alluncorrevents<<", OTHERalluncorrevents: "<<myOtherGrid->GetUncorrellatedEventNumber()<<std::endl;
    //cout<<"TEST: THISallevents: "<<allevents<<", OTHERmyOtherGrid->GetTotalEventNumber(): "<<myOtherGrid->GetTotalEventNumber()<<std::endl;
    // may be do this per grid
    //alluncorrevents+=myOtherGrid->GetUncorrellatedEventNumber();
    //allevents+=myOtherGrid->GetTotalEventNumber();
    //cout<<"TEST: THISalluncorrevents: "<<alluncorrevents<<std::endl;
    //cout<<"TEST: THISallevents: "<<allevents<<std::endl;


    //check if adding grids is safe
    int numGrids = myOtherGrid->GetNGrid();
    if ( numGrids != this->GetNGrid() ) {
     cout<<" MyGrid::AddGrid: WARNING: Number of grids are the not the same but will attemp to be added!"<<endl;
     exit(0);
    }

    //add grids
    for(int igrid=0; igrid<numGrids; igrid++) {
        *(GetGrid(igrid)) += *(myOtherGrid->GetGrid(igrid));
    }


    return;
}
*/
/*
void MyGrid::CorrectStatUncertainty(TH1D *h1){
  
  //cout<<" MyGrid::CorrectStatUncertainty before: "<<endl;
  //h1->Print("all");
  int nev=h1->GetEntries();
  for ( int ibin=0 ; ibin<=h1->GetNbinsX() ; ibin++ ){
    double  x = h1->GetBinContent(ibin);
    double ex = h1->GetBinError(ibin); //sqrt(Sum w^2)
    double err=ex*ex-x*x;
    err=sqrt(err)/sqrt(nev-1);
    h1->SetBinError(ibin,err);
  }
  //cout<<" MyGrid::CorrectStatUncertainty after: "<<endl;
  //h1->Print("all");

  return;
};
*/
