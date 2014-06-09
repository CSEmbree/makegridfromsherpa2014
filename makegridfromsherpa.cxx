/*
 * Title:   makegridfromsherpa
 * Purpose: Given SHERPA ntuples, we test appl_grid LO and NLO convolute
 * Authors: T. Carli,  C. Embree
 * with help of M. Sutton and Stefan Hoeche
 *
 * removed R,B logic since Sherpa2.1 does everything in one step
 */
//general
//
#include <iostream>
#include <string>
#include <vector>
#include <climits>

#include "root.h"
#include "t3.h"
//#include "utils.h"
//#include "normtmp.h"
#include "appl_grid/appl_grid.h"
#include "fjClustering.h"

#include "MyEvent.h"
#include "MyData.h"
#include "MyFrameData.h"
#include "MyObservable.h"
#include "MyGrid.h"
#include "MyHisto.h"

//root
//#include <TCanvas.h>
//#include <TH1D.h>
//#include <TFile.h>
//#include <TPad.h>

// extra
#include "LHAPDF/LHAPDF.h"

#define PI 3.141592653589793238462

/*
 * EXAMPLE Execution:
 *                   ./makegridfromsherpa
 *                   ./makegridfromsherpa <filename> <numevents>
 *                   ./makegridfromsherpa steering/atlas2012_top-config.txt 1000000
 */

//Constructs for more easily running over NTuples.
//
string ntup_names="";
long int nevmax=LONG_MAX; //allow for halt of execution when a user defined maximum number of events is reached
double wgt_fac=-1; //factor needed for conversion from SHERPA to ApplGrid when setting the weight 
static const int GLUON = 6; //used for more readable code

bool debug=false;
bool debug2=false;
bool NLOdebug=false;   //for activating NLO debug statements
bool nloworks=true;    //for testing each NLO piece 
//
//bool nloworks=false;
bool nloworks2=true;

extern "C" void evolvepdf_(const double& , const double& , double* );
extern "C" double alphaspdf_(const double& Q);

/*
//skeleton to replace hard coded and repetative NLO fillings
double PrepAndFillEvent(MyEvent *myevent, MyGrid *mygrid, double w, double x1, double x2, double x1p, double x2p, double id1, double id2, double xsec=0) {

    myevent->SetX1( x1 );
    myevent->SetX2( x2 );

    myevent->SetID( 0, id1 );   //setting parton A's ID
    myevent->SetID( 1, id2 );   //setting parton B's ID

    myevent->SetXSection( 0 ); //hack for not filling ref histos


    //determine subprocess that may be needed to divide by
    int currentSubProc = mygrid->GetDecideSubProcess( id1, id2 );
    int npairs         = mygrid->GetNSubProcessPairs( currentSubProc );

    //general weight calculation
    double wgt = ( w * wgt_fac ) / npairs;

    return wgt;
}
*/

/* to be done...

//calculates mins & maxs of ren & fac scales to get the correct Q2
void GetRenAndFacMaxAndMins(string NtupName, double *facMin, double *facMax, double *renMin, double *renMax) {

    TChain *fChaintmp= new TChain("t3");

    fChaintmp->Add(TString(NtupName));
    if(debug) fChaintmp->Print();

    t3 ttmp(fChaintmp);

    // determine Q2 and x boundaries from fac and ren scales
    Long64_t nentriestmp = fChaintmp->GetEntries();
    for ( Long64_t jentry=0 ; jentry<nentriestmp && jentry<nevmax ; jentry++ ) {
        ttmp.GetEntry(jentry);
        if (ttmp.fac_scale < *facMin) *facMin = ttmp.fac_scale;
        if (ttmp.fac_scale > *facMax) *facMax = ttmp.fac_scale;
        if (ttmp.ren_scale < *renMin) *renMin = ttmp.ren_scale;
        if (ttmp.ren_scale > *renMax) *renMax = ttmp.ren_scale;
    }

    if(debug) std::cout<<" makegridfromsherpa::GetRenAndFacMaxAndMins: ntuple:"<<NtupName
                           <<", facMin= "<<*facMin
                           <<", facMax= "<<*facMax
                           <<", renMin= "<<*renMin
                           <<", renMax= "<<*renMax<<std::endl;

    delete fChaintmp; //cleanup
}
*/

//makes conversion of id1, id2 and fills f[]
void getpdf(const double& x, const double& Q2, double* f) {

 evolvepdf_(x, Q2, f);
 //
 // ...for tests 
 // for(int id=1; id<12; id++) f[id]=x; // for testing
 // cout<<" getPDF x= "<<x<<" Q2= "<<Q2<<endl;
 // for(int id=0; id<13; id++) cout<<" getpdf: f["<<id<<"]= "<<f[id]<<endl;
}

//allows for getting value of environment variables
string GetEnv( const string & var ) {

 const char* res= getenv( var.c_str() );

 std::string s = res!=NULL? res:"";
 return s;
}

//method from Mark to divide two histograms nicely for ratio comparisons between two histos
TH1D* divide( const TH1D* h1, const TH1D* h2 ) {
 //
 //cout<<" makegridfromsherpa::divide: dividing two histos..."<<endl;
 //
 bool DBG=false;
 if ( h1==NULL || h2==NULL ) return NULL;

 TH1D* h = (TH1D*)h1->Clone();

 if (DBG) cout << " makegridfromsherpa::divide:histograms h1: " << h1->GetTitle() 
               << ", h2: " << h2->GetTitle() << endl;


 for ( int i=1 ; i<=h1->GetNbinsX() ; i++ ) {
  double b  = h2->GetBinContent(i);
  double be = h2->GetBinError(i);
  double t  = h1->GetBinContent(i);
  double te = h1->GetBinError(i);


  double r  = ( b!=0 ? t/b : 0. );
  double re = 0;

  if (DBG) cout << " makegridfromsherpa::divide b= "<<b<<" r= "<<r<<endl;
  h->SetBinContent( i, r );
  h->SetBinError( i, re ) ;
 }

 double hmin = h->GetBinContent(1);
 double hmax = h->GetBinContent(1);

 for ( int i=2 ; i<=h->GetNbinsX() ; i++ ) {
  double d = h->GetBinContent(i);
  if ( hmin>d ) hmin=d;
  if ( hmax<d ) hmax=d;
 }

 if (DBG) cout << " makegridfromsherpa::divide: \tmin ratio = " << hmin << "\tmax ratio = " << hmax << endl;
 cout<<" makegridfromsherpa::divide: h->GetMaximum(): "<<h->GetMaximum()<<", h->GetMinimum(): "<<h->GetMinimum()<<endl;

 /*
 if ( h->GetMaximum()<1.01 ) //h->SetMaximum(1.01);
 h->SetMaximum(0.99);
 if ( h->GetMinimum()>0.99 ) //h->SetMinimum(0.99);
 h->SetMinimum(1.01);
 */
 return h;
}

int main(int argc, char** argv) {

 // use a default atlas inclusive grid
 //attempt to use EnvVar to find steeringfile, otherwise use default folder path steering/ in current dir
 string steeringName = "atlas2012_top-config.txt"; //*-config for lumi, without for generic
 string steeringPath = "steering";
 string steeringFile = steeringPath + "/" + steeringName;

 string steeringDefaultPath = GetEnv("STEERINGPATH");

 if ( steeringDefaultPath.size() > 0 ) {
  steeringFile=steeringDefaultPath+"/"+steeringName;
  if (debug) cout<<" makegridfromsherpa::main: STEERINGPATH environment variable found, using path: "<<steeringFile<<endl;
 } else {
  if (debug) cout<<" makegridfromsherpa::main: STEERINGPATH environment varaible not set, using default: "<<steeringFile<<endl;
 }

 //attempt to use EnvVar to find PDFsets, otherwise use default folder path PDFsets/ in current dir
 // only needed for testing
 string pdfSetName = "CT10.LHgrid"; //"MSTW2008nlo90cl.LHgrid";
 // string pdfSetName = "MSTW2008nlo68cl.LHgrid";
 string pdfSetPath = "PDFsets";
 string pdfSetFile = pdfSetPath + "/" + pdfSetName;

 string pdfSetDefaultPath = GetEnv("LHAPATH");

 if (pdfSetDefaultPath.size() > 0) {
  pdfSetFile = pdfSetDefaultPath + "/" + pdfSetName;
  if (debug) cout<<" makegridfromsherpa::main: LHAPATH environment varaible found, using path: "<<pdfSetFile<<endl;
 } else {
  if (debug) cout<<" makegridfromsherpa::main: LHAPATH environment variable not set, using default: "<<pdfSetFile<<endl;
 }

 //allow for passing of another steering file name (default=inputname) and a different number of events(default uses all_events)
 if ( argc>1 ) {
  steeringFile = string(argv[1]);
  //if (debug) 
  cout << " makegridfromsherpa::main: Reading steering file " << steeringFile << endl;
 }

 if ( argc>2 ) {
  ntup_names = argv[2];
  //if (debug) 
  cout << " makegridfromsherpa::main: version= " << ntup_names.c_str() << endl;
 }

 if ( argc>3 ) {
  nevmax = atoi(argv[3]);
  //if (debug) 
  cout << " makegridfromsherpa::main: Reading Number of events " << nevmax << endl;
 }

//Create uniquely named grids for each Type: 0-B, 1-R, 2-RthenB
 TString newname=TString(ntup_names);
 if (newname.Contains("-*")) newname.ReplaceAll("-*","");
 string version=string(newname.ReplaceAll("_*",""));
 
 if (debug) cout<<" makegridfromsherpa::main: Creating grid using steeringFile: "<<steeringFile<<", version: "<<version<<endl;


 MyGrid *mygrid = new MyGrid(steeringFile, version);

 const int NGrid=mygrid->GetNGrid(); //NGrid happens to be the same for all Types

 cout<<"\n makegridfromsherpa  MyGrid booked NGrid= "<<NGrid<<" \n"<<endl;

 MyObservable *myobs[NGrid];
 for ( int igrid=0 ; igrid<NGrid ; igrid++ ) {
  MyData *mydata=mygrid->GetMyData(igrid);
  myobs[igrid]= new MyObservable(mydata);
  if (!myobs[igrid]) cout<<igrid<<" makegridfromsherpa MyObservable not created !"<<endl;
 }


 //create histograms to store tests:
 TH1D *htest1[NGrid];     // scaled by all events events and weight
 TH1D *htest2[NGrid];     // scaled by uncorr events and weight2
 TH1D *htest2lo[NGrid];   // LO only

 TH1D *htest3[NGrid];     // filled in NLO, reference for htest2 to make sure NLO is filled correctly
 TH1D *href[NGrid];       // holds histogram references from mygrid[x] to ensure they are equal to htest2[]
 TH1D *hreference[NGrid]; // holds histogram hreference from MyGrid
 MyHisto *hsigtot[NGrid];

 std::vector<int> nSubCount = std::vector<int>(25,0); //replace 25 by nprocs from lumi

 //maintain a count of the number of events run over for each histogram type we are testing
 int allevents=0;
 int alleventslo=0;
 int alluncorrevents=0;
 int alluncorreventslo=0;

 int idold=-1;
 double m_fsq=0.; 

 double sigtot=0.;  double sigtotsq=0.; 
 double sigtot2=0.; //double sigtot2sq=0.; 


 //
 // Loop over different histograms to test different *-like* approaches
 //
 //if (debug) 
 cout<<" makegridfromsherpa::main: prepare event loop "<<" \n"<<endl;
 TChain *fChain= new TChain("t3");

 //read in each NTuple for each test histogram
 string NtupName=mygrid->GetInputNtupDir()+"/"+mygrid->GetNtupName();
 NtupName += ntup_names;
 NtupName += ".root";

 //if (debug) 
 cout<<" makegridfromsherpa::main: Opening "<<NtupName<<" \n"<<endl;
        
 fChain->Add( TString( NtupName ) );
 if (debug) fChain->Print();

 t3 t(fChain);

 // set-up jet algorithm via fastjet wrapper class
 if (debug) cout<<" makegridfromsherpa::main: Set up jet algorithm:"<<endl;
 fjClustering* jetclus = new fjClustering( fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best );

 //create a new event for each event. This could be done better no?
 MyEvent* myevent = new MyEvent();

 Long64_t nentries = fChain->GetEntries();
 if (debug) cout<<" makegridfromsherpa::main: Number of events= "<<nentries<<endl;
 //
 // set-up test histogram
 //

 for ( int igrid=0 ; igrid<NGrid ; igrid++ )
 { //loop over grids

  string htest1name = "htest1" + version;
  htest1[igrid] = mygrid->GetHisto( igrid, htest1name );
  htest1[igrid]->SetTitle( htest1name.c_str() );
  htest1[igrid]->SetName ( htest1name.c_str() );
  htest1[igrid]->Sumw2();;
  if (debug) cout<<" makegridfromsherpa::main: Got histogram: "<<htest1name<<endl;

  string htest2name = "htest2" + version;;
  htest2[igrid] = mygrid->GetHisto( igrid, htest2name );
  htest2[igrid]->SetTitle( htest2name.c_str() );
  htest2[igrid]->SetName ( htest2name.c_str() );
  htest2[igrid]->Sumw2();;
  if (debug) cout<<" makegridfromsherpa::main: Got histogram: "<<htest2name<<endl; 

  string htest2namelo = "htest2lo" + version;;
  htest2lo[igrid] = mygrid->GetHisto( igrid, htest2namelo );
  htest2lo[igrid]->SetTitle( htest2namelo.c_str() );
  htest2lo[igrid]->SetName ( htest2namelo.c_str() );
  if (debug) cout<<" makegridfromsherpa::main: Got histogram: "<<htest2namelo<<endl;
            
  string htest3name = "htest3" + version;;
  htest3[igrid] = mygrid->GetHisto( igrid, htest3name );
  htest3[igrid]->SetTitle( htest3name.c_str() );
  htest3[igrid]->SetName ( htest3name.c_str() );
  htest3[igrid]->Sumw2();
  if (debug) cout<<" makegridfromsherpa::main: Got histogram: "<<htest3name<<endl;
 
  TString name="hsigtot"+ version;
  TH1D *htmp=(TH1D*)htest3[igrid]->Clone(name);
  //TH1D *htmp= new TH1D(name," tmp",1,0.,1.e6);
  htmp->SetTitle("total cross section with correct statistical uncertainties");
  hsigtot[igrid]= new MyHisto(htmp);

 }

 if (debug) cout<<" makegridfromsherpa::main: Loop over events for: "<<version<<endl;

 //
 // Event loop over all events
 //
 Long64_t nbytes = 0, nb = 0;
 LHAPDF::initPDFSet( pdfSetFile.c_str(), 0 ); //only needed for testing

 for ( Long64_t jentry=0 ; jentry<nentries && jentry<nevmax; jentry++ )
 {// loop over events
  t.GetEntry(jentry);
  nbytes += nb;
  //t.Show();
  Int_t np=(Int_t)t.nparticle;
  //inform user every NUMBER of events completed to illustrate progress.
  if (nevmax>nentries) nevmax=nentries;

  if (debug) cout<<" nevmax "<<nevmax<<" nentries= "<<nentries<<endl;

  if (jentry%100000==0 && jentry>0) {
   cout<<" Type: "<<version<<" events processed "<<jentry<<"  of  ";
   cout<<" TOTAL EVENTS: "<<nentries<<" done in "<<(float(jentry)/float(nevmax))*100.<<" %"<<endl;
  }

  if (debug) {
   // Documentation see
   // http://sherpa.hepforge.org/doc/SHERPA-MC-2.0.0.html
   //

   cout<<"\nmakegridfromsherpa::main: SHERPA EVENT INFO >>>>>>>>" <<endl;
   cout<<"  Event id = "<<t.id<<endl; //Event ID to identify correlated real sub-events.

   cout<<"  Incoming parton id1= "<<t.id1<<" id2= "<<t.id2<<endl; // PDG code of incoming parton 1/2
   cout<<"  Incoming parton number of particles= "<<np<<endl;
   for (int i=0; i<np; i++) {
    cout<<"  Outcoming parton kf["<<i<<"]= "<<t.kf[i]; //
    cout<<" px= "<<t.px[i]<<" pyin= "<< t.py[i]<<" pzin= "<<t.pz[i]<<" Ein= "<<t.E[i]<<endl;
   }
   cout<<"  x1 = "<<t.x1<<", x2 = "<<t.x2<<endl;     // Bjorken-x of incoming parton 1/2
   cout<<"  x1p= "<<t.x1p<<", x2p= "<<t.x2p<<endl;   // x’ for I-piece of incoming parton 1/2
   // Factorisation and Factorisaton scale.
   cout<<"  fac_scale= "<<t.fac_scale<<", ren_scale= "<<t.ren_scale
       <<" alphas= "<<t.alphas<<endl; //" alphasPower= "<<(Char_t)t.alphasPower<<endl;
   //
   // nuwgt Number of additional ME weights for loops and integrated subtraction terms.
   // usr_wgt[nuwgt] Additional ME weights for loops and integrated subtraction terms.
   //8.8.4.1 Computing (differential) cross sections of real correction events with statistical errors
   cout<<"  nuwgt= "<<t.nuwgt<<endl;
   for (int i=0; i<t.nuwgt; i++) cout<<" usr_wgts["<<i<<"]= "<<t.usr_wgts[i];
   cout<<endl;
   // units are GeV and pb
   // weight    Event weight, if sub-event is treated independently.
   // weight2   Event weight, if correlated sub-events are treated as single event.
   // see
   cout<<"  weight= "<<t.weight<<", weight2= "<<t.weight2<<endl;
   // me_wgt    ME weight (w/o PDF), corresponds to ’weight’.
   // me_wgt2   ME weight (w/o PDF), corresponds to ’weight2’.
   cout<<"  me_wgt= "<<t.me_wgt<<", me_wgt2= "<<t.me_wgt2<<endl;
   cout<<"  b_part= "<<t.part<<endl;
   cout<<"  b_nparticle= "<<t.nparticle<<endl; //Number of outgoing partons.
   // cout<<" np= "<<np<<endl;
   //Int_t np=t.b_nparticle;
   cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<endl;
  }

  int iorder=int(t.alphasPower);
  int itype =int(t.part[0]);

  //double fillweight=1.; 
  //if (t.id>50) fillweight=-1.; 
  double fillweight=t.weight2;  
  double fillme_wgt=t.me_wgt2;
  //double fillweight=t.weight;  
  //double fillme_wgt=t.me_wgt;

  jetclus->ClearJets();
  myevent->ClearEvent();;
  myevent->SetCMS(7000.);

  //prepare coefficients
  int Wsize=13;
  double xf1[Wsize];
  double xf2[Wsize];

  for (int i=0; i<Wsize; i++) {
   xf1[i]=0.0; //reset
   xf2[i]=0.0; //reset
  }
    
  // use evolvepdf_ to get weights in f1, f2
  getpdf( t.x1, t.fac_scale, xf1 );
  getpdf( t.x2, t.fac_scale, xf2 );

  //particle num convention conversion from SHERPA
  int id1 = t.id1;
  int id2 = t.id2;
  if (t.id1==21) id1=0;
  id1 = id1+6;
  if (t.id2==21) id2=0;
  id2 = id2+6;

  if (debug) cout<<" makegridfromsherpa::main: iorder: "<<iorder<<endl;

  //******FORCING REJECTION FOR TESTING
  //if (iorder<=2) continue; //Accepting only NLO
  //if (iorder!=2) continue; //Accepting only LO
  // if (!(id1==6 && id2==6)) continue; //Accepting only certain particles
  // if (id1==6 || id2==6) continue;
  //******FORCING REJECTION FOR TESTING

  //weights are now set explicitly bellow depending on event order being LO or NLO
  //myevent->SetWeight(wgt);
  myevent->SetCountEvent(true); 
  myevent->SetXSection( fillweight );
  myevent->SetOrder( iorder );
  myevent->SetType( itype );

  if (debug) cout<<" event-id= "<<t.id <<endl;

  myevent->SetEventId( t.id );

  allevents++; 
  sigtot  +=t.weight;
  sigtotsq+=t.weight*t.weight;

  sigtot2  +=fillweight;
  if (t.id!=idold) {
   idold=t.id;
   if (debug) cout<<" new event id= "<<t.id<<endl;
   sigtot2=fillweight;
  }
  //cout<<t.id<<" weight= "<<fillweight<<" tweight= "<<sigtot2<<endl;

  if (iorder==2) alleventslo++; 
  //cout<<" New Event ?  "<<myevent->NewEvent()<<endl;
  //hsigtot->Fill(obs,fillweight,myevent->NewEvent());
  bool newevent=myevent->NewEvent();
  if (myevent->NewEvent()) {
   alluncorrevents++; 
   //cout<<" New Event: alluncorrevents= "<<alluncorrevents<<" allevents= "<<allevents<<endl;
   if (iorder==2) alluncorreventslo++; 
  }



  myevent->SetQ2( t.fac_scale * t.fac_scale ); // applgrid takes sqrt(fac_scale) when evaluating pdf and alphas

  //
  // fill incoming partons
  //
  double pxin=0.;
  double pyin=0.;
  double pzin=0.;
  double Ein=0.;

  for ( int ip=0 ; ip<np ; ip++ ) {
   pxin += t.px[ip];
   pyin += t.py[ip];
   pzin += t.pz[ip];
   Ein  += t.E[ip];
  }

  if (debug) cout<<" pxin= "<<pxin<<" pyin= "<<pyin<<" pzin= "<<pzin<<" Ein= "<<Ein<<endl;

  double shat = Ein*Ein - ( pxin*pxin + pyin*pyin + pzin*pzin );
  double s    = shat / ( t.x1*t.x2 );
  double sqrs=0.;
  if (s>0.) sqrs=sqrt(s);
  myevent->SetCMS(sqrs);

  //if (debug) cout<<" s= "<<sqrs<<endl;

  double ep = ( sqrs / 2.0 );

  int pid   = t.id1;
  if (pid==21) pid = 0; //conversion from sherpa gluon to appl_grid convention
  myevent->push_back( pxin, pyin, t.x1 * ep, t.x1 * ep, pid);

  pid=t.id2;
  if (pid==21 ) pid = 0; //conversion from sherpa gluon to appl_grid convention
  myevent->push_back( pxin, pyin, -t.x2 * ep, t.x2 * ep, pid );
  // cout<<" makegridfromsherpa::main: pid2: "<<pid<<endl;

  for (int ip=0; ip<np; ip++)
  {
   // Momentum components of the partons  kf  Parton PDG code
   pid=t.kf[ip];
   if (pid==21) pid=0; //conversion from sherpa gluon to appl_grid convention

   //if (abs(pid)==6) {
   myevent->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
   //}

   if (abs(pid)==11||abs(pid)==12) {
    myevent->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
    continue;
   } else {
    // need to put in Gavins code here for pid
    // for top that will only work in LO
    jetclus->push_back(t.px[ip],t.py[ip],t.pz[ip],t.E[ip],pid);
   }
  }

  //if (debug) myevent->Print();

  bool mycut[NGrid];
  //double obs[NGrid];
  double obs=0.; //for time being only one grid
  for ( int igrid=0 ; igrid<NGrid ; igrid++ ) {
   mycut[igrid]=myobs[igrid]->eventcuts(myevent);
   if (mycut[igrid]==false) {
    cout<<igrid<<" makegridfromsherpa event rejected !"<<endl;
    cout<<" makegridfromsherpa cut by event cut no filling"<<endl;
    continue;
   }
   obs=myobs[igrid]->GetObservable(myevent);
   if (debug) cout<<" makegridfromsherpa igrid= "<<igrid<<" obs= "<<obs<<endl;
  }


  //double obs = myevent->GetInvariantMassTop12(); //only needs to be computed once per events with the same id

  // Need subprocess ID to compute the number of pairs as part of a combo. Refer to Lumi_pdf implimentation
  int subProcID = mygrid->GetDecideSubProcess( t.id1==21? 0:t.id1, t.id2==21? 0:t.id2 );
  double npairs = mygrid->GetNSubProcessPairs( subProcID );

  if (debug) cout<<" makegridfromsherpa subProcID= "<<subProcID<<" npairs= "<<npairs<<endl;

  // weight multiplication factor needed by ApplGrid when setting Grid weights when filling mygrid
  wgt_fac = pow( ( 2.0 * PI ) / t.alphas, iorder );

  // needed for filling reference histos

  int nWgts=13;   //number of flavors tBar(0)...gluon(6)...t(13)
  // coefficients used in Stefan's code.
  double f1  [nWgts];
  double f2  [nWgts];
  double xf1p[nWgts];
  double xf2p[nWgts];
  for ( int iproc=0 ; iproc<nWgts ; iproc++ ) {
   xf1[iproc] = 0.0; //emptied
   xf2[iproc] = 0.0; //emptied
  }

  // evolve to get weights in xf1, xf2
  getpdf( t.x1, t.fac_scale, xf1 );
  getpdf( t.x2, t.fac_scale, xf2 );

  // Filling f1 and f2 now for easy access in Stefan's code
  for ( int iproc=0 ; iproc<nWgts ; iproc++ ) {
   f1[iproc] = xf1[iproc] / t.x1;
   f2[iproc] = xf2[iproc] / t.x2;
  }

  //computing f1 from xf1 and f2 from xf2 for later usage (and readability) conveyance
  double fa = f1[id1]; // (xf1/x1)[id1]
  double fb = f2[id2]; // (xf2/x2)[id2]
  double lr = 0.0;     //NOTE: NOT SET CURRENTLY!
  double lf = 0.0;     //NOTE: NOT SET CURRENTLY!

  if (debug) cout<<" makegridfromsherpa fa= "<<fa<<" fb= "<<fb<<endl;

  double wgtCheck = 0.;

  // fill weight depending on LO or NLO
  if (iorder==2) {
   //
   // iorder = 2 = LO
   //
   // compute weight to fill my event with when LO
   double wgt = (fillme_wgt * wgt_fac) / npairs;
   if (debug) cout<<" fillme_wgt: "<<fillme_wgt<<" fa= "<<fa<<" fb= "<<fb<<endl;
   wgtCheck = fillme_wgt * fa * fb;

   //if (debug) cout<<" LO wgtCheck: "<<wgtCheck<<", t.weight2: "<<t.weight2<<", iorder: "<<iorder<<endl;
   if (debug) cout<<" LO wgtCheck: "<<wgtCheck<<", t.weight2: "<<t.weight2<<", iorder: "<<iorder<<endl;
   if (debug) cout<<" LO wgt: "<<wgt<<" x1= "<<t.x1<<" x2= "<<t.x2<<endl;
   if (fabs(wgtCheck-fillweight) >1.e-3) cout<<" XXX LO problem weight  XXX "
                                           <<" wgtCheck= "<<wgtCheck<<" fillweight= "<<fillweight<<endl;
   //set event based on LO
   myevent->SetWeight(wgt); 
   myevent->SetX1(t.x1);
   myevent->SetX2(t.x2);

   if (nloworks) {
     mygrid->fill(myevent,obs);
    myevent->SetCountEvent(false); 
    if (debug2) cout<<" makegridfromsherpa filling LO obs= "<<obs<<" wgtCheck= "<<wgtCheck<<endl;
    for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill(obs,wgtCheck);                
   }
  } else { // 
   //
   // iorder = 3 = NLO
   //
   int currentSubProc = -1; //holder for result of deciding the sub process of two partons
   // w used for setting weights
   if (t.nuwgt==0) {
    wgtCheck = fillme_wgt * fa * fb;
    double wgt = (fillme_wgt * wgt_fac) / npairs;

    //if (NLOdebug) cout<<" NLO-V t.me_wgt2: "<<t.me_wgt2<<" fa= "<<fa<<" fb= "<<fb<<endl;
    if (NLOdebug) cout<<" NLO-V fillme_wgt: "<<fillme_wgt<<" fa= "<<fa<<" fb= "<<fb<<endl;
    if (NLOdebug) cout<<" NLO-V wgtCheck: "<<wgtCheck
                      <<" t.weight2: "<<t.weight2
                      <<" t.weight: " <<t.weight
                      <<" iorder: "<<iorder<<endl;
    if (NLOdebug) cout<<" NLO-V wgt: "<<wgt<<" x1= "<<t.x1<<" x2= "<<t.x2<<endl;
    if (fabs(wgtCheck-fillweight) >1.e-3) cout<<" XXX NLO-V problem weight  XXX "
                                            <<" wgtCheck= "<<wgtCheck<<" fillweight= "<<fillweight<<endl;

    myevent->SetWeight(wgt); //dealing entirely with weight2
    myevent->SetX1(t.x1);
    myevent->SetX2(t.x2);

    if (nloworks) {
      mygrid->fill(myevent,obs);
     myevent->SetCountEvent(false); 
      
     if (debug2) cout<<" makegridfromsherpa filling NLO nuwgt=0 obs= "<<obs<<" wgtCheck= "<<wgtCheck<<endl;
     for (int igrid=0; igrid<NGrid; igrid++) htest3[igrid]->Fill(obs,wgtCheck);
    }
   } else {

     //lr= lr=log(t.ren_scale/sqrt(ren_scale)), lf=log(muf2/sqr(p_vars->m_muf));
    
    double w[9];
    // here me_wgt is correct !
    w[0] = t.me_wgt + t.usr_wgts[0] * lr + t.usr_wgts[1] * lr * lr / 2.0; 
    if (debug) cout<<" t.me_wgt: "<<t.me_wgt<<" w[0]= "<<w[0]<<" nuwgt= "<<t.nuwgt<<endl;

    bool wnz=false;
    for ( int i=1 ; i<9 ; ++i ) {
     if (t.nuwgt<=2) 
      w[i] = 0.;
     else 
      w[i] = t.usr_wgts[i+1] + t.usr_wgts[i+9] * lf;
     if (w[i]==0) wnz=true;
    }

    wgtCheck=w[0]*fa*fb;

    double  wgt = 0.0; //weight for event
    wgt      = ( w[0] * wgt_fac ) / npairs;

    if (debug) cout<<" wgt: "<<wgt<<" w[0]= "<<w[0]<<endl;

    myevent->SetWeight( wgt );
    myevent->SetX1( t.x1 );
    myevent->SetX2( t.x2 );

    if (nloworks) { // XXX 1
      mygrid->fill( myevent,obs );
     myevent->SetCountEvent(false);  
     for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill(obs,wgtCheck);     
     if (debug2) cout<<" makegridfromsherpa filling NLO XXX1 w[0] obs= "<<obs<<" w= "<<wgtCheck<<endl;   
    }

    if (NLOdebug) cout<<"TEST: id1: "<<id1<<", id2: "<<id2<<", wnz: "<<wnz<<endl;

    if (wnz==true)
    {
     wgt = 0.0; //reset

     double  faq  = 0.0;
     double  faqx = 0.0;
     double  fag  = 0.0;
     double  fagx = 0.0;
     double  fbq  = 0.0;
     double  fbqx = 0.0;
     double  fbg  = 0.0;
     double  fbgx = 0.0;
     //
     // calculating coefficients to fill ref histos before filling/setting up mygrid
     // literally a copy of Stefan's code.
     //
     if ( id1 != GLUON ) {
      //Parton 1 = QUARK
      faq = fa;
      fag = f1[GLUON];

      getpdf( t.x1/t.x1p, t.fac_scale, xf1p );

      faqx = xf1p[id1] / t.x1;
      fagx = xf1p[6]   / t.x1;
     } else {
      // Parton 1 - GLU
      fag=fa;
      for ( int i=1 ; i<nWgts-1 ; ++i)
       if( i!=GLUON ) faq += f1[i];

      getpdf( t.x1/t.x1p, t.fac_scale, xf1p );

      fagx = (xf1p[id1] / t.x1);
      for ( int i=1 ; i<nWgts-1 ; ++i )
       if ( i!= GLUON ) faqx += (xf1p[i]/t.x1);
     }
     if ( id2 != GLUON ) {
      // Parton 2 - QUARK
      fbq = fb;
      fbg = f2[GLUON];

      getpdf( t.x2/t.x2p, t.fac_scale, xf2p );

      fbqx = xf2p[id2] / t.x2;
      fbgx = xf2p[GLUON] / t.x2;
     } else {
      // Parton 2 - GLU
      fbg = fb;
      for ( int i = 1 ; i<nWgts-1 ; ++i)
       if( i != GLUON ) fbq += f2[i];
 
      getpdf( t.x2/t.x2p, t.fac_scale , xf2p);

      fbgx = (xf2p[id2] / t.x2);
      for ( int i=1 ; i<nWgts-1 ; ++i )
       if ( i != GLUON ) fbqx += (xf2p[i] / t.x2);
     }

     if (NLOdebug) {
      cout<<" faq= "<<faq<<" faqx= "<<faqx<<" fag= "<<fag<<" fagx= "<<fagx<<endl;
      cout<<" fbq= "<<fbq<<" fbqx= "<<fbqx<<" fbg= "<<fbg<<" fbgx= "<<fbgx<<endl;
      for (int i=0; i<9; i++) cout<<" w["<<i<<"]= "<<w[i]<<endl;
     }

     if (NLOdebug) cout<<"NLO-I wgtCheck: "<<wgtCheck<<endl;

     wgtCheck += ( faq*w[1] + faqx*w[2] + fag*w[3] + fagx*w[4] ) * fb;
     wgtCheck += ( fbq*w[5] + fbqx*w[6] + fbg*w[7] + fbgx*w[8] ) * fa;

     if (NLOdebug) cout<<"NLO-I wgtCheck: "<<wgtCheck
                       <<" t.weight2: "<<t.weight2
                       <<" t.weight:  "<<t.weight
                       <<endl;
     if (fabs(wgtCheck-fillweight) >1.e-3) cout<<" XXX NLO-I problem weight  XXX "
                                               <<" wgtCheck= "<<wgtCheck<<" fillweight= "<<fillweight
					      <<" Ratio= "<<wgtCheck/fillweight
                                               <<endl;
     //if (NLOdebug) cout<<"TEST: jentry: "<<jentry<<", id1: "<<id1<<" , id2: "<<id2<<endl;

     //
     // Filling grid and ref histos from coeficients and 'w' just generated
     //
     if ( id1 != GLUON ) {
      //
      // parton 1 = QUARK
      //
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos
      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );
      wgt = ( w[1] * wgt_fac ) / npairs; //<--weights determined
      myevent->SetWeight( wgt );
      if (nloworks) {// XXX 2
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 
       
       //if (debug2) cout<<" makegridfromsherpa filling NLO XXX2 w[1] obs= "<<obs<<" fill= "<<fb*faq*w[1]<<endl;
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX2 w[1] obs= "<<obs<<" w= "<< fb*faq*w[1]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*faq*w[1] ); 
      }
			
      myevent->SetX1( t.x1 / t.x1p );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos
      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );
			
      wgt = ( w[2] * ( 1.0 / t.x1p ) * wgt_fac ) / npairs; //<--weights determined
      myevent->SetWeight( wgt );
      if (nloworks) { // XXX 3
	mygrid->fill( myevent, obs );
       myevent->SetCountEvent(false); 
      
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX3 w[2] obs= "<<obs<<" w= "<<fb*faqx*w[2]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*faqx*w[2]);
      }
			
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, 0 ); // term with first parton is gluon
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos
      currentSubProc = mygrid->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = ( w[3]  * wgt_fac ) / npairs; //<--weights determined
      myevent->SetWeight( wgt );			

      if (nloworks) {
	mygrid->fill( myevent, obs );
       myevent->SetCountEvent(false); 
       // 
       if (debug2) cout<<" makegridfromsherpa filling NLO w[3] obs= "<<obs<<" w= "<<fb*fag*w[3]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*fag*w[3] );
      }			

      myevent->SetX1( t.x1 / t.x1p );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, 0 );
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos

      currentSubProc = mygrid->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = ( w[4] * ( 1.0 / t.x1p ) * wgt_fac ) / npairs; //<--weights determined
      myevent->SetWeight( wgt );
      if (nloworks) { // XXX 5 
	mygrid->fill( myevent, obs );
       myevent->SetCountEvent(false);  
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX4 w[4] obs= "<<obs<<" w= "<<fb*fagx*w[4]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*fagx*w[4] );
      }			
     } else {
      //
      // parton 1 = GLU
      //
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
 
      // count the number of subprocs that make up
      for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav )
      {
       if( iflav != GLUON )
       {
        currentSubProc = mygrid->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );
        if (currentSubProc != -1) nSubCount[currentSubProc]++; //count the number of occurences of that parlticular subproc
       }
      }

      // filling grids for all subprocs
      //myevent->SetCountEvent(true); 
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav )
      {
       if ( iflav != GLUON )
       {
        currentSubProc = mygrid->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );

        if ( currentSubProc != -1 )
        {
         npairs = mygrid->GetNSubProcessPairs( currentSubProc );
         // should be changed according to 7
         wgt = w[1] * wgt_fac * ( 1.0 / npairs ) * ( 1.0 / nSubCount[currentSubProc] ); //<--weights determined
         myevent->SetID( 0, (iflav-6) );
         myevent->SetID( 1, (t.id2==21? 0:t.id2) );
         myevent->SetWeight( wgt );
         myevent->SetXSection(0.); //hack for not filling ref histos
         // XXX6
         if (nloworks) {
	   mygrid->fill( myevent, obs );
          myevent->SetCountEvent(false);
         } 
         //myevent->SetCountEvent(false); 
        }
       }
      }
      //myevent->SetCountEvent(true); 
      if (nloworks) {    
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX6 w[1] obs= "<<obs<<" w= "<<fb*faq*w[1]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*faq*w[1]);    
      }

      // set x1 and x2 for all of the upcoming subprocs event fills
      myevent->SetX1( t.x1 / t.x1p );
      myevent->SetX2( t.x2 );

      // count the number of subprocs that make up
      for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );
        if( currentSubProc != -1 ) nSubCount[currentSubProc]++;                                
       }
      }

      // filling grids for all subprocs
      //myevent->SetCountEvent(true); 
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (iflav-6), (t.id2==21? 0:t.id2) );
        if (NLOdebug) 
         cout<<"P2-22: flv1: "<<(iflav-6)<<", flv2: "<<(t.id2==21? 0:t.id2)<<", DecideSubProc: "<<currentSubProc<<endl;
        if ( currentSubProc != -1 ) {
         npairs = mygrid->GetNSubProcessPairs( currentSubProc );                                   
         // wgt = w[2] * wgt_fac * ( 1.0 / t.x1p ) * ( 1.0 / npairs) /* * ( 1.0 / nSubCount[currentSubProc])*/; //<--weight determine				    
         wgt = w[2] * wgt_fac * ( 1.0 / t.x1p ) * ( 1.0 / npairs)  /* * ( 1.0 / nSubCount[currentSubProc])*/; //<--weight determine
         if (NLOdebug)
          cout<<"XX7 : flv1: "<<(iflav-6)<<", flv2: "<<(t.id2==21? 0:t.id2)
              <<", DecideSubProc: "<<currentSubProc
              <<" wgt= "<<wgt<<" w[2]= "<< w[2]<<" wgt_fac= "<<wgt_fac<<" 1/x1p= "<<( 1 / t.x1p )
              <<" npair= "<<npairs<<" nsub= "<<nSubCount[currentSubProc]
              <<endl;                                    
         myevent->SetID( 0, (iflav-6) );
         myevent->SetID( 1, (t.id2==21? 0:t.id2) );
         myevent->SetWeight( wgt );
         myevent->SetXSection(0.); //hack for not filling ref histos
         // XXX 7 
         if (nloworks) {
	   mygrid->fill( myevent,obs );
          myevent->SetCountEvent(false); 
         }
        }
       }
      }

      //myevent->SetCountEvent(true); 
      if (nloworks){ 
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX7 w[2] obs= "<<obs<<" w= "<<fb*faqx*w[2]<<endl;
       if (NLOdebug) cout<<" htest3 wgt= "<<fb*faqx*w[2]<<" fb= "<<fb<<" faqx= "<<faqx<< " w[2]= "<<w[2]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*faqx*w[2] );          
      }
                        
      //glu components 
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, 0 ); //glu
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      currentSubProc = mygrid->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = w[3] * wgt_fac * ( 1.0 / npairs);//<--weight determined
      myevent->SetWeight( wgt );
			
      if (nloworks){ // xxx 8
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false);  
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX8 w[3] obs= "<<obs<<" w= "<<fb*fag*w[3]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*fag*w[3] );
      }			
                        
      myevent->SetX1( t.x1 / t.x1p );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, 0 ); //glu
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
                        
      currentSubProc = mygrid->GetDecideSubProcess( 0, (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = w[4] * wgt_fac * ( 1.0 / t.x1p ) * ( 1.0 / npairs);//<--weight determined
      myevent->SetWeight( wgt );

      if (nloworks){ // 9
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX9 w[4] obs= "<<obs<<" w= "<<fb*fagx*w[4]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fb*fagx*w[4] );
      }			
     }

     if ( id2 != GLUON ) {
     //
     // parton 2 = QUARK
     //
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos
      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = ( w[5] * wgt_fac ) / npairs;//<--weight determined

      myevent->SetWeight( wgt );
      if  (nloworks) {// XX 10
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 

       if (debug2) cout<<" makegridfromsherpa filling NLO XXX10 w[5] obs= "<<obs<<" fill= "<<fa*fbq*w[5]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbq*w[5] );
      }

      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 / t.x2p );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, (t.id2==21? 0:t.id2) );
      myevent->SetXSection( 0 ); //hack for not filling ref histos

      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (t.id2==21? 0:t.id2) );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = ( w[6] * ( 1.0 / t.x2p ) * wgt_fac ) / npairs;//<--weight determined

      myevent->SetWeight( wgt );
      if (nloworks) { // XXX 11
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX11 w[6] obs= "<<obs<<" w= "<<fa*fbqx*w[6]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbqx*w[6] );
      }	
		
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, 0 );
      myevent->SetXSection( 0 ); //hack for not filling ref histos

      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );
  
      wgt = ( w[7]  * wgt_fac ) / npairs;//<--weight determined
                       
      myevent->SetWeight( wgt );
			
      if (nloworks) { //
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 
       if (debug2) cout<<" makegridfromsherpa filling NLO w[7] obs= "<<obs<<" w= "<<fa*fbg*w[7]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbg*w[7] );
      }	
		
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 / t.x2p );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, 0 );
      myevent->SetXSection( 0 ); //hack for not filling ref histos

      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = ( w[8] * ( 1.0 / t.x2p ) * wgt_fac ) / npairs;//<--weight determined
      //
      myevent->SetWeight( wgt );
			
      if (nloworks) {// XXX 13
	mygrid->fill( myevent,obs );
       myevent->SetCountEvent(false); 
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX13 w[8] obs= "<<obs<<" w= "<<fa*fbgx*w[8]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbgx*w[8] );
      }
     } else {
      //
      // parton 2 = GLU
      //
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );

      //count the number of subprocs that make up
      for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );
        if ( currentSubProc != -1 ) nSubCount[currentSubProc]++;
       }
      }

      // filling grids for all subprocs
      //myevent->SetCountEvent(true); 
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );

        if ( currentSubProc != -1 ) {
         npairs = mygrid->GetNSubProcessPairs( currentSubProc );
         wgt = w[5] * wgt_fac * ( 1.0 / npairs ) * ( 1.0 / nSubCount[currentSubProc] );//<--weight determined

         myevent->SetID( 0, (t.id1==21? 0:t.id1) );
         myevent->SetID( 1, (iflav-6) );
         myevent->SetWeight( wgt );

         if (NLOdebug) cout<<"P2-21: filling ..."<<endl;

         myevent->SetXSection( 0 ); //hack for not filling ref histos
         // XXX 14
         if (nloworks) {
	   mygrid->fill( myevent, obs );
          myevent->SetCountEvent(false); 
         }
        }
       }
      }

      //myevent->SetCountEvent(true); 
      if (nloworks) {
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX14 w[5] obs= "<<obs<<" w= "<<fa*fbq*w[5]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbq*w[5] );
      }

      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 / t.x2p );
      // count the number of subprocs that make up
      for ( int i=0 ; i<nSubCount.size(); ++i ) nSubCount[i] = 0;//empty for next use
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );
        if ( currentSubProc != -1 ) nSubCount[currentSubProc]++;
       }
      }

      // filling grids for all subprocs
      // myevent->SetCountEvent(true); 
      for ( int iflav=1 ; iflav<nWgts-1 ; ++iflav ) {
       if ( iflav != GLUON ) {
        currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), (iflav-6) );
        if (NLOdebug) cout<<"P2-22: flv1: "<<(iflav-6)<<", flv2: "<<(t.id2==21? 0:t.id2)<<", DecideSubProc: "<<currentSubProc<<endl;
        if ( currentSubProc != -1 ) {
         npairs = mygrid->GetNSubProcessPairs( currentSubProc );                                    
         wgt = w[6] * wgt_fac * ( 1 / t.x2p )*( 1.0 / npairs)  /* ( 1.0 / nSubCount[currentSubProc])*/;//<--weight determined
                                   
         myevent->SetID( 0, (t.id1==21? 0:t.id1) );
         myevent->SetID( 1, (iflav-6) );
         myevent->SetWeight( wgt );
         if (NLOdebug) cout<<"P2-22: filling ..."<<endl;
         myevent->SetXSection( 0 ); //hack for not filling ref histos
         // XXX 15
         if (nloworks) {
	   mygrid->fill( myevent,obs );
          myevent->SetCountEvent(false); 
         }
        }
       }
      }
      //myevent->SetCountEvent(true); 
      if (nloworks) {
       if (debug2) cout<<" makegridfromsherpa filling NLO XXX15 w[6] obs= "<<obs<<" w= "<<fa*fbqx*w[6]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbqx*w[6] );   
      }

      // glu component
      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, 0 ); //glu

      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = w[7] * wgt_fac * ( 1.0 / npairs);//<--weight determined
      //
      myevent->SetWeight( wgt );
			
      if (nloworks){// XXX 16
	mygrid->fill( myevent, obs );
       myevent->SetCountEvent(false); 

       if (debug2) cout<<" makegridfromsherpa filling NLO XXX16 w[7] obs= "<<obs<<" w= "<<fa*fbg*w[7]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbg*w[7] );
      }			

      myevent->SetX1( t.x1 );
      myevent->SetX2( t.x2 / t.x2p );
      myevent->SetID( 0, (t.id1==21? 0:t.id1) );
      myevent->SetID( 1, 0 ); //glu

      currentSubProc = mygrid->GetDecideSubProcess( (t.id1==21? 0:t.id1), 0 );
      npairs = mygrid->GetNSubProcessPairs( currentSubProc );

      wgt = w[8] * wgt_fac * ( 1.0 / t.x2p ) * ( 1.0 / npairs);//<--weight determined
      //
      myevent->SetWeight( wgt );
      if (nloworks){	// XXX 17
	mygrid->fill( myevent, obs );
        myevent->SetCountEvent(false); 

       if (debug2) cout<<" makegridfromsherpa filling NLO XXX17 w[8] obs= "<<obs<<" w= "<<fa*fbgx*w[8]<<endl;
       for ( int igrid=0 ; igrid<NGrid ; igrid++ ) htest3[igrid]->Fill( obs, fa*fbgx*w[8] );
      }			
     }
     if (NLOdebug) cout<<"test: Current event number finishing "<<jentry<<endl;
    }
   }
  }
  //
  // fill each test histogram
  //

  //sigtot2sq=sigtot2*sigtot2;

  if (debug) {
   //cout<<" makegridfromsherpa fillweight = "<<fillweight<<" tweight= "<<sigtot2<<endl;
   cout<<" makegridfromsherpa sigtot2 = "<<sigtot2<<" sigtot= "<<sigtot<<endl;
   cout<<" makegridfromsherpa sigtot2/nall= "<<sigtot2/alluncorrevents
       <<" newevent= "<<newevent<<endl;
  }

  for ( int igrid=0 ; igrid<mygrid->GetNGrid() ; igrid++ )
  {

   hsigtot[igrid]->Fill(obs,fillweight,newevent);

   if (debug) cout<<" makegridfromsherpa fill htest event= "<<allevents<<" obs= "<<obs
		//<<" t.weight2= "<<t.weight2<<endl;
		  <<" t.weight= "<<t.weight<<endl;
   htest1[igrid]->Fill( obs, t.weight );            
   htest2[igrid]->Fill( obs, fillweight );             

   if (debug) {
    cout<<" makegridfromsherpa print htest1 event= "<<allevents<<endl;
    htest1[igrid]->Print("all");            

    cout<<" makegridfromsherpa print htest2 event= "<<alluncorrevents<<endl;
    htest2[igrid]->Print("all");            

    cout<<" print htest3 "<<endl;
    htest3[igrid]->Print("all");            

    //cout<<" print applgrid reference "<<endl;
    //TH1D *htmp=(TH1D*)mygrid->GetReference(0);
    //htmp->Print("all");
    
   }

   if (iorder==2) {
    htest2lo[igrid]->Fill( obs, fillweight );                
   }
  } //end loop over grid
 } //end loop over events


 
 if (debug) cout<<"\n makegridfromsherpa::main: Finished running over events!\n"<<endl;
 //if (debug) {
  cout<<" makegridfromsherpa events all= "<<allevents<<" LO= "<<alleventslo
      <<" alluncorr= "<<alluncorrevents<<" uncorrLO= "<<alluncorreventslo<<endl;

  double sigerr=sigtotsq/allevents-pow(sigtot/allevents,2);
  sigerr=sqrt(sigerr/(allevents-1.));
  cout<<" makegridfromsherpa sum over all weights (from all events): "<<endl;
  cout<<" makegridfromsherpa sigtot= "<<sigtot/allevents<<" sigerr= "<<sigerr<<endl;
  //}

 //
 // get and set up test external histograms
 //
 for (int igrid=0; igrid<mygrid->GetNGrid(); igrid++)
 {

  href[igrid]=(TH1D*)mygrid->GetReference(igrid);
  if (!href[igrid]) {
   cout<<" makegridfromsherpa::main: Reference from grid not found ! "<<endl;
   exit(0);
  } else {
   TString nameTitle="applgridreference"+version;
   href[igrid]->SetTitle(nameTitle);
   href[igrid]->SetName(nameTitle);
   href[igrid]->SetLineColor(7); 
  }

  hreference[igrid]=(TH1D*)mygrid->GetReferenceHisto(igrid);
  if (!hreference[igrid]) {  
   cout<<" makegridfromsherpa::main: Reference from MyGrid not found ! "<<endl;
  } else {
   TString nameTitle="hreference"+version;
   hreference[igrid]->SetTitle(nameTitle);
   hreference[igrid]->SetName(nameTitle);
   hreference[igrid]->SetLineColor(8); 
  }
 }
//} //end of loop over all histograms

 //write grid to scale the internal reference histos and save them to *.root files
 if (debug) cout<<" makegridfromsherpa Scaling and writing mygrid "<<endl;
 bool phasespacerun=false;
 cout<<" makegridfromsherpa call write_grid "<<endl;
 mygrid->write_grid();
 phasespacerun=   mygrid->IsPhaseSpaceRun(); 
 if (phasespacerun) {
  cout<<" makegridfromsherpa is phase space run return "<<endl;
  return 0;
 }

 if (debug) cout<<" makegridfromsherpa Grids written "<<endl;

 //
 // Scale htest histos
 // Add histograms of MyGrid for R-Like(i_R) to B-Like(i_B). Should be equal to both combined(i_RB)
 //
 if (debug) cout<<" makegridfromsherpa::main: normalise and adding R-Type and B-Type"<<endl;

 for (int igrid=0; igrid<NGrid; igrid++) {
  //scale each histograms by one devided by number of events depending on type
  cout<<" makegridfromsherpa: allevents= "<<allevents<<endl;
  if (allevents!=0) {
   htest1[igrid]->Scale( 1.0 / allevents ); 
  }
  cout<<" makegridfromsherpa: alluncorrevents= "<<alluncorrevents<<endl;
  if (alluncorrevents!=0) {
   htest2[igrid]->Scale( 1.0 / alluncorrevents ); 
   htest3[igrid]->Scale( 1.0 / alluncorrevents ); 
  }
  //if (alluncorreventslo!=0) htest2lo[igrid]->Scale(1.0 / alluncorreventslo); //
  if (alluncorrevents!=0) htest2lo[igrid]->Scale(1.0 / alluncorrevents); //

  //cout<<" sigtot from histo htest1= "<<htest1[igrid]->Integral()<<endl;
  //cout<<" sigtot from histo htest2= "<<htest2[igrid]->Integral()<<endl;

  //hsigtot[igrid]->GetHisto()->Scale( 1.0 / alluncorrevents ); 
  hsigtot[igrid]->AddLastEvent();
 
  cout<<" makegridfromsherpa hsigtot2 = "<<  hsigtot[igrid]->GetSum()
                                         <<" +- "<< hsigtot[igrid]->GetSumError()<<endl;
  
  hsigtot[igrid]->Print();

  //if (debug) {
  cout<<" makegridfromsherpa::main: htest1 printing..."<<endl;
  htest1[igrid]->Print("all");
  cout<<" makegridfromsherpa::main: htest2 printing..."<<endl;
  htest2[igrid]->Print("all");
  cout<<" makegridfromsherpa::main: htest3 printing..."<<endl;
  htest3[igrid]->Print("all");
  //}

  cout<<" makegridfromsherpa::main: htest1 printing after correction error..."<<endl;
  //mygrid->CorrectStatUncertainty(htest1[igrid]);
  htest1[igrid]->Print("all");

  //cout<<" makegridfromsherpa::main: hnorm printing (htest1/binswidth)... "<<endl;
  TString hnormname=TString("hnorm")+htest1[igrid]->GetName();
  TH1D *hnorm1=(TH1D*)htest1[igrid]->Clone(hnormname);
  mygrid->DivideByBinWidth(hnorm1);

  TString hnormname2=TString("hnorm2")+htest2[igrid]->GetName();
  TH1D *hnorm2=(TH1D*)htest2[igrid]->Clone(hnormname2);
  mygrid->DivideByBinWidth(hnorm2);

  TString hnormnamelo=TString("hnorm")+htest2lo[igrid]->GetName();
  TH1D *hnorm2lo=(TH1D*)htest2lo[igrid]->Clone(hnormnamelo);
  mygrid->DivideByBinWidth(hnorm2lo);
  //if (debug){
   cout<<" makegridfromsherpa::main: hnorm2 printing..."<<endl;
   hnorm2->Print("all");
   cout<<" makegridfromsherpa::main: href printing..."<<endl;
   href[igrid]->Print("all");
  //}
 }

 //
 // Convolute and print all external histograms
 //
 if (debug) cout<< "\n makegridfromsherpa::main: Performing convolute"<<endl;
 //NGrid = mygrid[0]->GetNGrid(); //NGrid will be the same for all grids, so grid[0] is arbirary
  TH1D* convGridHistos[NGrid];
  TH1D* LOconvGridHistos[NGrid];
  TH1D* subProcConvGridHistos[NGrid][121];
  TH1D* LOsubProcConvGridHistos[NGrid][121];

  int nLoops = 1;
  //creating an identity matrix needed in convolute function - temp printing to make sure it's correct
  std::vector<std::vector<double> > ckm2 = std::vector<std::vector<double> >(13, std::vector<double>(13,0));
  for(int x=0; x<ckm2.size(); x++) {
   for(int y=0; y<ckm2.at(x).size(); y++) {
    if (x==y) ckm2.at(x).at(y)=1;
   }
  }

  TFile *fout;

  for (int igrid=0; igrid<NGrid; igrid++)
  {

   MyData *mydata=mygrid->GetMyData(igrid);
   bool binw=mydata->DivideByBinWidth();

   string filename=mygrid->GetDataName(igrid);
   filename+=version+"_histos.root";
   fout= new TFile(filename.c_str(),"recreate");
   mygrid->GetGrid(igrid)->setckm(ckm2);

   //NLO convolute for all
   //convolute histograms do not need to be normalised !
   cout<<" makegridfromsherpa calling NLO convolute: "<<endl;
   convGridHistos[igrid] = (TH1D*)mygrid->GetGrid(igrid)->convolute( getpdf, alphaspdf_, nLoops );
   convGridHistos[igrid]->SetName((TString) ("NLOconvolute" + version));
   convGridHistos[igrid]->SetTitle((TString) ("NLOconvolute"+ version));
   convGridHistos[igrid]->SetLineColor(kBlue);
   convGridHistos[igrid]->Write();

   cout<<" makegridfromsherpa::main: Printing NLO convolute histo..."<<endl;
   convGridHistos[igrid]->Print("all");

   //LO convolute for all
   cout<<" makegridfromsherpa calling LO convolute: "<<endl;
   LOconvGridHistos[igrid] = (TH1D*)mygrid->GetGrid(igrid)->convolute( getpdf, alphaspdf_, 0 );
   LOconvGridHistos[igrid]->SetName((TString) ("LOconvolute" + version));
   LOconvGridHistos[igrid]->SetTitle((TString) ("LOconvolute"+ version));
   LOconvGridHistos[igrid]->SetLineColor(kBlue);
   LOconvGridHistos[igrid]->Write();

   cout<<" makegridfromsherpa::main: Printing LO convolute histo..."<<endl;
   LOconvGridHistos[igrid]->Print("all");

   for (int isubproc=0; isubproc<mygrid->GetNSubProcess(igrid);  isubproc++) {
    //NLO convolute for subprocs
    subProcConvGridHistos[igrid][isubproc] = (TH1D*)mygrid->GetGrid(igrid)->convolute_subproc(isubproc, getpdf, alphaspdf_, nLoops );
    string sub_proc_hist_name="subProc-"+to_string(isubproc)+"-convolute"+version;
    subProcConvGridHistos[igrid][isubproc]->SetName((TString) (sub_proc_hist_name));
    subProcConvGridHistos[igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name));
    subProcConvGridHistos[igrid][isubproc]->SetLineColor(kBlue);

    subProcConvGridHistos[igrid][isubproc]->Write();

    //LO convolute for subprocs
    LOsubProcConvGridHistos[igrid][isubproc] = (TH1D*)mygrid->GetGrid(igrid)->convolute_subproc(isubproc, getpdf, alphaspdf_, 0 );
    sub_proc_hist_name="LOsubProc-"+to_string(isubproc)+"-convolute_for"+version;
    LOsubProcConvGridHistos[igrid][isubproc]->SetName((TString) (sub_proc_hist_name));
    LOsubProcConvGridHistos[igrid][isubproc]->SetTitle((TString) (sub_proc_hist_name));
    LOsubProcConvGridHistos[igrid][isubproc]->SetLineColor(kBlue);

    LOsubProcConvGridHistos[igrid][isubproc]->Write();

   }

   href[igrid]->Write(); //scaled version

   // output origional htest1 
   htest1[igrid]->Write();
   // output origional htest2 
   htest2[igrid]->Write();

   //string name1 = "hnorm1"+ntup_names;
   //TH1D *htest1norm=(TH1D*)htest1[igrid]->Clone( name1.c_str() );
   //htest1norm->SetTitle( name1.c_str() );
   //htest1norm->SetName(  name1.c_str() );
   //if (binw) mygrid->DivideByBinWidth(htest1norm); //DivideByBinWidth version
   //cout<<" makegridfromsherpa::main: Printing htest2 normalised..."<<endl;
   //htest2norm->Print("all");
      
   string name2lo = "hnorm2lo"+version;
   TH1D *htest2lonorm=(TH1D*)htest2lo[igrid]->Clone( name2lo.c_str() );
   htest2lonorm->SetTitle( name2lo.c_str() );
   htest2lonorm->SetName(  name2lo.c_str() );
   if (binw) mygrid->DivideByBinWidth(htest2lonorm); //DivideByBinWidth version
   cout<<" makegridfromsherpa::main: Printing htest2lo normalised..."<<endl;
   htest2lonorm->Print("all");

   TH1D* ratioLO = divide( LOconvGridHistos[igrid], htest2lonorm); //LOconvolute/htestnorm
   if (ratioLO) {
    string ratioName = "ratio_LO-conv/htest2lonorm"+version;
    ratioLO->SetTitle( ratioName.c_str() );
    ratioLO->SetName ( ratioName.c_str() );
    ratioLO->Print("all");

    ratioLO->DrawCopy();
    gPad->Print("ratioLO.pdf");
   }

   // output htest3 & a DivideByBinWidth version
   htest3[igrid]->Print("all");
   htest3[igrid]->Write(); 

   string name3 = "hnorm3"+version;
   TH1D *htest3norm=(TH1D*)htest3[igrid]->Clone( name3.c_str() );
   htest3norm->SetTitle( name3.c_str() );
   htest3norm->SetName(  name3.c_str() );

   if (binw) mygrid->DivideByBinWidth(htest3norm); //DivideByBinWidth version
   cout<<" makegridfromsherpa::main: Printing htest3 normalised..."<<endl;
   htest3norm->Print("all");
   cout<<" makegridfromsherpa::main: Printing convoluted normalised..."<<endl;
   convGridHistos[igrid]->Print("all");

   TH1D* ratioNLO = divide( convGridHistos[igrid], htest3norm ); //NLOconvolute/htestnorm
   if ( ratioNLO ) {
    string ratioName = "ratio_NLO-conv/htest3norm"+version;
    ratioNLO->SetTitle( ratioName.c_str() );
    ratioNLO->SetName(  ratioName.c_str());
    ratioNLO->Print("all");

    ratioNLO->DrawCopy();
    gPad->Print("ratioNLO.pdf");
   }

   fout->Write(); //writing root file
   fout->Close();
  }

  if (debug) cout<<" makegridfromsherpa::main: Internal Reference histograms normalised!"<<endl;
  return 0;
}
