#ifndef MyGrid_H
#define MyGrid_H

#include <stdlib.h> // exit()
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector>
#include <map>

using std::string;
#include <sys/time.h> 

#include "root.h"
#include "MyHisto.h"
#include "MyEvent.h"
#include "MyData.h"
#include <appl_grid/appl_grid.h>
//#include <appl_grid/appl_igrid.h>
//#include "appl_grid/generic_pdf.h" //TEST-generic
//#include "appl_grid/basic_pdf.h" //TEST-basic
#include "appl_grid/lumi_pdf.h" //TEST-lumi

//#include "appl_grid/appl_pdf.h"

//static const int __nf__   = 5;
//static const int __nf2__  = 11;


class MyGrid {
 private:

  bool debug; 
  bool bookrefsubprocess;

  bool phasespacerun;
  //long unsigned int alluncorrevents;
  //long unsigned int allevents;
  std::vector<long unsigned int> uncorrevents; // per grid
  std::vector<long unsigned int> events;

  bool newevent;
  int eventid;
  string glabel;

  string ntupdiroutput;
  string ntupdirinput;
  string ntupname;
  //  string filename;
  string subprocesssteername;

  string gridnamedir;
  string steername;
  string versionname; //<--added as optional argument
  
  bool fixedinputgrid;
  string fixedinputgridname;  

  int processnumber; // process number foresee setting from reaction
  int nsub; // number of subprocesses
 
  double *weight;
  // grid architecture -> could be in steering file
  //
  int    nXbins;
  double xLow;
  double xUp;
  int    xorder;
  int    nQ2bins;
  double q2Low; 
  double q2Up;
  int    qorder;
  int    iorder;

  bool incljets;
  
  // lowest order in alphas	
  int lowest_order;
  // how many loops
  int nloops;

  string pdf_function;

  std::vector<int> istylemarkerdata; // marker of theory prediction
  std::vector<int> icolormarkerdata; // marker of theory prediction
  std::vector<appl::grid*> mygrid;     // grid vector
  std::vector<bool> isBooked;          // flag grid already booked 

  std::vector<string> gridname;

  typedef std::vector<MyHisto*> vhref;

  std::vector<vhref> hrefLOsubprocesseshistos;      //saved 
  std::vector<vhref> hrefsubprocesseshistos;        //saved 
  vhref hreference;                                 //saved 
  vhref hreferencefine;                             //saved 
 
  std::vector<MyData*> mydata;        // information about data from steering file

 lumi_pdf *mypdf; //TEST-lumi
 //basic_pdf *mypdf; //TEST-basic
 //generic_pdf *mypdf; //TEST-generic
  
 

 public:
 //MyGrid(char name[100]);
  MyGrid(string name, string version="");
  
  //void AddGrid(MyGrid *myOtherGrid); //add the reference histograms of another grid to this one.

  //int GetSubProcess() { return mypdf->GetCurrentSubProcess();};
  //string GetSubProcessName(int isub) { return mypdf->GetSubProcessName(isub);}; //not implimented by basic_pdf
  //void SetSubProcess(generic_pdf *subpro) { mypdf=subpro; return;}; //TEST-generic
  //void SetSubProcess(basic_pdf *subpro) { mypdf=subpro; return;}; //TEST-basic
  void SetSubProcess(lumi_pdf *subpro) { mypdf=subpro; return;}; //TEST-lumi
  int GetNSubProcess(int igrid)
  { 
    int subproci=0; // default subprocesses argument
    return mygrid[igrid]->subProcesses(subproci); 
  };

  //new GetNSubProcess for overloaded choice in applgrid's subProcesses method
  int GetNSubProcess(int igrid, int subproci)
  { 
    return mygrid[igrid]->subProcesses(subproci); 
  };

  int GetNSubProcessPairs(int iproc){ 
    if(mypdf)
        return (*mypdf)[iproc].size();
    else {
        std::cout<<"MyGrid::GetNSubProcessPairs: ERROR: mypdf not set-up yet!"<<std::endl;
        exit(0);
    } };
    
  int GetDecideSubProcess(int iflav1, int iflav2) {
    if(mypdf)
        return mypdf->decideSubProcess( iflav1, iflav2 );
    else {
        std::cout<<"MyGrid::GetDecideSubProcess: ERROR: mypdf not set-up yet!"<<std::endl;
        exit(0);
    } };
  

  void Initialize();
  //void ReadSteering(char fname[100]);
  void ReadSteering(string);
  bool file_exists(const string& s);
  void book_grid(int igrid);  // inital grid booking  
  void fill(MyEvent *myevent, double obs);
  //bool eventcuts(MyEvent *myevent,int igrid);
  void write_grid();   // writes out grid after some events
  void DivideByBinWidth(TH1D* h);
  void Normalise(TH1D* h1, double yscale, double xscale, bool binw, bool normtot); 
  void NormaliseRefHistos(int igrid, int nev);
  TH1D *TH1NormStatError(TH1D *hsum, TH1D *hsum2, double norm);
  //void CorrectStatUncertainty(TH1D *h1);
  //void ScaleInternalRefHistos(int igrid);
  //void getObservable(MyEvent *myevent);
  int GetNGrid(){return gridname.size();};
  string *GetSubProcessSteername() { std::cout<<" MyGrid::GetSubProcessSteername: TEST: subprocesssteername: '"<<subprocesssteername<<"'"<<std::endl; &subprocesssteername; };
  string GetNtupDirOutput(){ return ntupdiroutput;};
  string GetInputNtupDir() { return ntupdirinput;};

  void SetInputNtupDir(TString dirinput){ ntupdirinput=dirinput; return;};
  string GetNtupName(){ return ntupname;};
  
  bool IsPhaseSpaceRun() {return phasespacerun;};
  
  void SetQ2Low (double _q2Low)     { q2Low=_q2Low; }; 
  void SetQ2Up  (double _q2Up)      { q2Up=_q2Up; }; 
  void SetQ2bins(double _nQ2bins)   { nQ2bins=_nQ2bins; };

  //  string GetGridName(int i){
  // TString name=TString(gridname[i]); 
  // name.ReplaceAll(".txt","");
  // return string(name.Data());
  //};
/*
  string SetFixedGridFileName(string filename){
   if (debug) cout<<" MyGrid::SetFixedGridFileName set filename= "<<filename.c_str()<<endl;
   fixedinputgrid=true;
   fixedinputgridname=filename;
  }
  void SetFixedGridFileName(){
   fixedinputgrid=true;
   if (debug) cout<<" myGrid SetFixedGridFileName fixedinputgrid is ON "<<endl;
   //fixedinputgridname=filename;
  }
*/
  string GetGridFileName(int i){
   TString name=TString(this->GetDataName(i));
   //name.ReplaceAll(".txt","");
   name+=".root";
   //return string(this->GetGridName(i).Data())+".root";
   return string(name);
  }
  string GetGridVersionName(){ //<--added accessor method to versionname
    return versionname;
  }
  
  string SetGridVersionName(string version){ //<--added mutator method to versionname
    return versionname=version;
  }

  string GetGridFullFileName(int i){
    return this->GetNtupDirOutput()+"/"+this->GetGridFileName(i);
  }

  MyData *GetMyData(int i){ return mydata[i];};

  int GetProcessNumber() {return processnumber;};
  void SetProcessNumber(int nproc) {processnumber=nproc; return;};
  
  //int GetNumberSubProcess() {return nsub;};
  //void SetNumberSubProcess(int nproc) {nsub=nproc; return;};

  string GetGridDir(){ return gridnamedir;};
  //string GetDataName(int i){ return gridnamedir+"/"+gridname[i];};
  string GetDataName(int i){ 
   TString name=gridname[i];
   name.ReplaceAll(".txt","");
   return string(name.Data());
  };

  int GetStyleMarker(int i) {
   if (istylemarkerdata.size()<1) {cout<<" MyGrid::GetStyleMarker No style marker set "<<endl; return 0;}
   return istylemarkerdata[i];
  };
  void SetStyleMarker (int i) {istylemarkerdata.push_back(i); return;};
  int GetColorMarker(int i) {
   if (icolormarkerdata.size()<1) {cout<<" MyGrid::GetColorMarker No color marker set "<<endl; return 0;}
   return icolormarkerdata[i];
  };
  void SetColorMarker (int i) {icolormarkerdata.push_back(i); return;};

  double *GetWeight(){return weight;};
  double GetWeight(int i){return weight[i];};
  void SetWeight(int i,double wgt){weight[i]=wgt; return;};
  void PrintWeight(){ 
   cout<<" MyGrid::PrintWeight: nsub= "<<nsub<<" (zero suppressed print out)"<<endl;
   for (int i=0; i<nsub; i++) 
    if (weight[i]!=0.)
     cout<< " weight["<<i<<"]="<<weight[i]<<endl; 
   return;
  }  

  void ResetWeight(){ for (int i=0; i<nsub; i++) weight[i]=0.; return;};

  appl::grid *GetGrid(int igrid){
   if (!mygrid[igrid]) cout<<" MyGrid::GetGrid mygrid["<<igrid<<"] not filled ! "<<endl;
   return mygrid[igrid];
  }


  
 
  //START - Accessor methods for all reference histograms
  std::vector<std::vector<MyHisto*> > *GetRefLOsubprocessesHistos(){
   if (hrefLOsubprocesseshistos.size()<=0) {
    cout<<" MyGrid::GetRefLOsubprocessesHistos: WARNING: no histograms found!"<<endl;
   }
   return &hrefLOsubprocesseshistos; 
  }
  
  std::vector<std::vector<MyHisto*> > *GetRefsubprocessesHistos(){
   if (hrefsubprocesseshistos.size()<=0) {
    cout<<" MyGrid::hrefsubprocesseshistos: WARNING: no histograms found!"<<endl;
   }
   return &hrefsubprocesseshistos; 
  }

  TH1D *GetReference(int igrid);  
   
  TH1D* GetReferenceHisto(int igrid) {
   if (!bookrefsubprocess) {
    cout<<"MyGrid::GetReferenceHistos: WARNING: no hreference booked "<<endl;
    return 0;
   }
   if ((int)hreference.size()<=igrid) {
    cout<<" MyGrid::GetReferenceHistos: WARNING: number of histograms found smaller igrid!"<<endl;
   }
   return hreference[igrid]->GetHisto(); 
  }

  TH1D* GetReferenceFineHisto(int igrid) {
   if ((int)hreferencefine.size()<=igrid) {
    cout<<" MyGrid::GetReferenceFineHistos: WARNING: number of histograms found smaller igrid!"<<endl;
   }
   return hreferencefine[igrid]->GetHisto(); 
  }
 
  /*
  std::vector<TH1D*> *GetReferenceHistos() {
   if (hreference.size()<=0) {
    cout<<" MyGrid::GetReferenceHistos: WARNING: no histograms found!"<<endl;
   }
   return &hreference; 
  }

  std::vector<TH1D*> *GetReferenceFineHistos() {
   if (hreferencefine.size()<=0) {
    cout<<" MyGrid::GetReferenceFineHistos: WARNING: no histograms found!"<<endl;
   }
   return &hreferencefine; 
  }
  //END - Accessor methods for all reference histograms
  */
  
 double GetObsMax(int igrid){
  if (!mydata[igrid])  cout<<" MyGrid::GetObsMax: grid not found= mydata["<<igrid<<"]"<<endl;
  //int nobs=mydata[igrid]->GetNBins(); 
  //cout<<" GetObsMax nobs= "<<nobs<<endl;
  //double *xbins=mydata[igrid]->GetBins();
  double obsmax=mydata[igrid]->GetMaxX(); 
  //if (debug) cout<<" obsmax= "<<obsmax<<endl;
  return obsmax;
 }

  TH1D *GetHisto(int igrid,string name="htest"){
   if (!mydata[igrid]) cout<<" MyGrid::GetHisto mydata not filled ! "<<endl;
   int nobs=mydata[igrid]->GetNBins(); 
   //if (debug) cout<<" MyGrid::GetHisto: nobs= "<<nobs<<", for igrid: "<<igrid<<endl;
   double *xbins=mydata[igrid]->GetBins();

   //cout<<" dataname= "<<this->GetDataName(igrid)<<" version= "<<this->GetGridVersionName()<<endl; 
   TString hname=TString(name)+TString("_")+TString(this->GetDataName(igrid))+TString(this->GetGridVersionName()); 
   if (debug) cout<<" MyGrid::GetHisto: hname set to: "<<hname.Data()<<endl;

   TH1D *htest = new TH1D(hname.Data(),hname.Data(),nobs,xbins);
   if (!htest) cout<<" htest["<<igrid<<"] not found "<<endl;
   //htest->Print("all");
   return htest;
  }; 


  //double GetUncorrellatedEventNumber(){return alluncorrevents;};
  double GetUncorrellatedEventNumber(int igrid){return uncorrevents[igrid];};
  //double GetTotalEventNumber(){return allevents;};
  double GetTotalEventNumber(int igrid){return events[igrid];};

  void Print();
  void PrintSubprocessRefHistos(int igrid);
  void PrintRefHistos(int igrid);
  //void NormRefHistos(int igrid, double norm);

};


#endif
