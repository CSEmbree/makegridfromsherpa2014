#include <iostream>
using namespace std;
#include <string>
#include "MyHisto.h"
/******************************************************************
 ** Method Implementations
 ******************************************************************/
//MyHisto::MyHisto() 
 MyHisto::MyHisto(TH1D *h1) 
{

  if (debug) cout<<" MyHisto::constructor: book "<<h1->GetName()<<endl;
  this->Initialise();

  if (!h1) cout<<" MyHisto::constructor: histogram not found "<<endl;

  TString name="MyHisto";
  name+=h1->GetName();
  TString htit=h1->GetTitle();

  TString name2=name;

  hsum=(TH1D*)h1->Clone(name);
  hsum->SetTitle(htit);

  name2+="sumw2";
  TString htit2=htit;
  htit2=htit;
  htit2+="sum of weight squares";
  hsum2=(TH1D*)h1->Clone(name2);
  hsum2->SetTitle(htit2);

  name+="tmp";
  htmp=(TH1D*)h1->Clone(name);

  name2+="tmp";
  htit2+="tmp ";
  htmp2=(TH1D*)h1->Clone(name2);
  htmp2->SetTitle(htit2);

  return;
}

 MyHisto::MyHisto(string name, string htit, int nbin, double xmin, double xmax)
{
  if (debug) cout<<" MyHisto constructor: book "<<name.c_str()<<endl;
  this->Initialise();

  string name2=name;
  string htit2=htit;

  hsum= new TH1D(name.c_str(),htit.c_str(),nbin,xmin,xmax);
  name2+="sumw2";
  htit2+="sum of weight squares";
  hsum2= new TH1D(name2.c_str(),htit2.c_str(),nbin,xmin,xmax);

  name+="tmp";
  htmp= new TH1D(name.c_str(),htit.c_str(),nbin,xmin,xmax);
  htmp->Sumw2();

  name2+="tmp";
  htit2+="tmp ";
  htmp2= new TH1D(name2.c_str(),htit2.c_str(),nbin,xmin,xmax);

  return;
}

void  MyHisto::Initialise(){

  debug=false; 
  if (debug) cout<<" MyHisto constructor: Initialise "<<endl;

  nev=0;  nevall=0;
  newevent=false;

  wev=0.; 
  wsum=0.; wsum2=0.; 
 
}
void MyHisto::Fill(double obs, double w, bool newevent)
{
  //cout<<" MyHisto::Fill pointer hsum "<<hsum<<endl;
  // set
  // http://sherpa.hepforge.org/doc/SHERPA-MC-2.0.0.html#Structure-of-ROOT-NTuple-Output
  // to understand the strategy how to deal with correllated events
  //
  nevall++;

  if (debug) {
   cout<<" MyHisto::Fill    ..filling w = "<<w<<endl;
   cout<<" MyHisto::Fill    ..filling we= "<<wev<<endl;
   cout<<" MyHisto::Fill oldevent nev= "<<nev<<endl;
   cout<<" MyHisto::Fill newevent w= "<<w  <<" wev= "<<wev<<endl;
  }

  if (newevent) {
   nev++;
   if (debug) cout<<" MyHisto::Fill newevent nev= "<<nev<<endl;

   //if (nev!=2) { //? derived from Stefan code
   if (nev>1) { 
    wsum+=wev;    
    wsum2+=wev*wev; 

    if (debug) {
     cout<<" before adding "<<endl;
     hsum ->Print("all");
     hsum2->Print("all");
     htmp->Print("all");
    }
 //
    hsum ->Add(htmp);
    htmp->Multiply(htmp,htmp);
    hsum2->Add(htmp);
    if (debug) {
     cout<<" before adding "<<endl;
     hsum ->Print("all");
     hsum2 ->Print("all");
    }
   }

   htmp ->Reset();
   wev=0.;

  } 
  wev+=w; 
  htmp ->Fill(obs,w);

  if (debug) 
   cout<<" MyHisto::Fill          wsum= "<<wsum<<" wsum2= "<<wsum2<<endl;

  if (debug) {
   cout<<" print hsum "<<endl;
   hsum->Print("all");
   cout<<" print htmp "<<endl;
   htmp->Print("all");
  }
};

void MyHisto::AddLastEvent(){
  wsum+=wev;    
  wsum2+=wev*wev; 

  hsum ->Add(htmp);
  htmp->Multiply(htmp,htmp);
  hsum2->Add(htmp);

  wev=0.; 
  htmp->Reset();

  if (debug) {
   cout<<" MyHisto::AddLastEvent nev= "<<nev<<endl;
   hsum->Print("all");
   hsum2->Print("all");
  }

  for ( int ibin=1; ibin<=hsum->GetNbinsX(); ibin++ ){
   double  x  = hsum->GetBinContent(ibin)/nev;
   //double ex = hsum->GetBinError(ibin); 
   double  x2 = hsum2->GetBinContent(ibin)/nev;
   double err=x2*x2-x*x;
   if (debug) 
    cout<<" MyHisto::AddLastEvent ibin= "<<ibin<<" x= "<<x<<" x2= "<<x2<<endl;

   if (err<0) cout<<" MyHisto::LastBin something is wrong err<0 err= "<<err<<endl;
   err=sqrt(err/(nev-1));
   hsum->SetBinError(ibin,err);
  }
  return;
}

double MyHisto::GetSumError(){

 debug=true;
 if (debug)    cout<<" MyHisto::GetSumError nev= "<<nev<<" nevall= "<<nevall<<endl; 
 if (nev<2) { 
  cout<<" MyHisto::GetSumError nev too small nev= "<<nev<<endl; 
  return -1.;
 }
 //if (debug) { 
  cout<<" MyHisto::GetSumError wsum= "<<wsum<<" wsum2= "<<wsum2<<" nev= "<<nev<<endl; 
  cout<<" MyHisto::GetSumError <wsum>= "<<wsum/nev<<" <wsum2>= "<<wsum2/nev<<" nev= "<<nev<<endl; 
  //}
 
 double xev=nev;
 double sigerr=wsum2/(xev)-pow(wsum/xev,2);
 if (sigerr<0)   cout<<" MyHisto::GetSumError something is wrong sigerr= "<<sigerr<<endl; 
 sigerr=sqrt(sigerr/(xev-1.));

 if (debug) 
  cout<<" MyHisto::GetSumError sigerr= "<<sigerr<<endl; 

 return sigerr;
};
void MyHisto::Print(){
 
 cout<<" MyHisto::Print hsum: "<<hsum->GetName()<<endl;
 hsum->Print("all"); 

};
