//
//   for sherpa
//
#include <iostream>
using namespace std;
#include <string>

#include "MyObservable.h"


/******************************************************************
 ** Method Implementations
 ******************************************************************/

MyObservable::MyObservable(MyData *data)
{

  debug=false;

  if (debug) cout<<" MyObservable called "<<endl;

  mydata=data;
  //myevent=event;

  incljets=mydata->isObsInclJets(); 
  observable=mydata->GetObservable(); 
  if (debug) cout<<" MyObservable observable= "<<observable<<endl;

  obsptjet=false;
  obsrapjet=false;
  if (observable.Contains("PTJET")) obsptjet=true;
  if (observable.Contains("RAPJET")) obsrapjet=true;
  if (obsptjet) cout<<" MyObservable observable is ptjet"<<endl;
  if (obsrapjet)cout<<" MyObservable observable is rapidity jet"<<endl;

  return;
}

//  std::vector<fastjet::PseudoJet*> myseljets=myevent->GetSelectedJets(ptmin,rapmin,rapmax);

double MyObservable::GetJetObservable(fastjet::PseudoJet *myjet){
 cout<<" MyObservable::GetJetObservableValue observable= "<<observable<<endl;
 if (!myjet)  cout<<" MyObservable::GetJetObservableValue myjet not found "<<endl;
 double obs=0.;
 if (obsptjet) obs=myjet->pt();
 if (obsrapjet) obs=myjet->rap();
 return obs;
}
double MyObservable::GetObservable(MyEvent *myevent) {
  double obs=0.;
  // need to find different solution here
 /*
  int id1=2;
  int id2=3;
  double p3[4];
  double p4[4];
  for (int  i = id1; i < id2; i++) {
    if ((myevent->Get(i))->user_index()!=6 ) cout<<" not top ? "<<myevent->Get(i)->user_index()<<endl;
  }
  p3[0]=(myevent->Get(id1))->e();
  p3[1]=(myevent->Get(id1))->px();
  p3[2]=(myevent->Get(id1))->py();
  p3[3]=(myevent->Get(id1))->pz();

  p4[0]=(myevent->Get(id2))->e();
  p4[1]=(myevent->Get(id2))->px();
  p4[2]=(myevent->Get(id2))->py();
  p4[3]=(myevent->Get(id2))->pz();
   
  double mass34 = 0;
  mass34 = std::sqrt( std::pow(p3[0] + p4[0],2) - std::pow(p3[1] + p4[1],2) - std::pow(p3[2] + p4[2],2) - std::pow(p3[3] + p4[3],2) );

  cout<<" MyObservable::GetObservable obs from mass34 = "<<mass34<<endl;
 */
  obs = myevent->GetInvariantMassTop12();
  if (debug) cout<<" MyObservable::GetObservable obs from GetInvariantMass = "<<obs<<endl;

  return obs;
}

bool  MyObservable::eventcuts(MyEvent *myevent) {
//
// analysis cuts defined in steering
//
    bool mydebug=false;

    bool accept=false;
    bool onejet=false;
    bool oklept=false;
    bool okneut=false;

    int idjetmax=-1;
    double ptold=-999.;
    int idlep=-1;
    int idneu=-1;
    if (!mydata) cout<<" MyObservable::eventcuts MyData object not found !"<<endl;

    int njets=0;
    if (mydebug)cout<<" MyObservable::eventcuts number of particle N= "<<myevent->GetN()<<" "<<endl;

    for (int  i = 0; i < myevent->GetN(); i++) {
        int id = myevent->GetID(i);
        if (mydebug) cout<<" MyObservable::eventcuts";
        double pt = (myevent->Get(i))->pt();
        double rap = fabs((myevent->Get(i))->rap());
        if (mydebug)
            printf(" %d pt= %6.2f rap= %6.2f  id= %d \n",i,pt, rap, id);

        bool isjet = myevent->IsJet(id);

        if (mydebug)
            if( mydata->JetCut()) cout<<" MyObservable::eventcuts JetCut on ! "<<endl;
            else                  cout<<" MyObservableGrid::eventcuts JetCut off !"<<endl;

        if (mydebug&&isjet) cout<<" MyObservable::eventcuts particle "<<i<<" is a jet ! "<<endl;
        if (mydata->JetCut()&&isjet) {
            if ((rap<mydata->GetJetMaxY()&&rap>mydata->GetJetMinY())
                    && pt>mydata->GetJetPtCut()) {
                onejet=true;
                if (mydebug) cout<<" MyObservable::eventcuts accept jet "<<i<<endl;
                njets++;
                if (pt>ptold) {
                    ptold=pt;
                    idjetmax=i;
                }
            };
        };

        //bool islepton=abs(id)==11||abs(id)==13;
        bool islepton=myevent->IsChargedLepton(id);
        if (mydebug&&islepton) cout<<" MyObservable::eventcuts is a lepton "<<i<<endl;
        if (mydata->LeptonCut()&&islepton) {
            idlep=i;
            if ((rap<mydata->GetLepMaxY()&&rap>mydata->GetLepMinY())
                    &&pt>mydata->GetLepPtCut()) oklept=true;
            if (mydebug)
                if (oklept) cout<<" MyObservable::eventcuts accept lepton cut "<<endl;
                else        cout<<" MyObservable::eventcuts lepton cut not passed "<<endl;
        };

        bool isneutrino=myevent->IsNeutrino(id);
        if (mydebug&&isneutrino) cout<<" MyObservable::eventcuts is a neutrino "<<i<<endl;
        if (mydata->NeutrinoCut()&&isneutrino) {
            idneu=i;
            if (pt>mydata->GetNeuPtCut()) okneut=true;
            if (mydebug) if (okneut) cout<<" MyObservable::eventcuts accept neutrino cut "<<endl;
        };
    };

    if (mydebug)
        if (okneut) cout<<" MyObservable::eventcuts neutrino cut passed "<<endl;
        else if (mydata->NeutrinoCut()) cout<<" MyObservable::eventcuts neutrino cut not passed "<<endl;

    if (mydebug)
        if (oklept) cout<<" MyObservable::eventcuts lepton cut passed "<<endl;
        else if (mydata->LeptonCut()) cout<<" MyObservable::eventcuts lepton cut not passed "<<endl;

    if (mydebug)
        if (onejet) cout<<" MyObservable::eventcuts jet cut passed "<<endl;
        else if (mydata->JetCut()) cout<<" MyObservable::eventcuts jet cut not passed "<<endl;

    if (mydata->JetCut()) accept=onejet;
    else accept=true;

    if (mydebug&&accept) cout<<" MyObservable::eventcuts event passed 1 "<<endl;

    if (accept) {
        if (mydata->NeutrinoCut()&&!okneut) {
            if (debug) cout<<" MyObservable::eventcut cut by neutrino "<<endl;
            return false;
        }
        if (mydata->LeptonCut()  &&!oklept) {
            if (debug) cout<<" MyObservable::eventcut cut by lepton "<<endl;
            return false;
        }
    }

    if (mydebug&&accept) cout<<" MyObservable::eventcuts event passed 2 "<<endl;

    if (accept&&mydata->MtwCut()) {
        if (mydebug) cout<<" MyObservable::eventcuts idlep= "<<idlep<<" idneu= "<<idneu<<endl;
        if (idlep!=-1&&idneu!=-1) {
            fastjet::PseudoJet pw=*myevent->Get(idlep)+*myevent->Get(idneu);;
            double mtw=pw.mt();
            if (mydebug) cout<<" MyObservable::eventcuts mtw= "<<mtw<<endl;
            if (mtw<mydata->GetMtWCut()) return false;
        };
    };

    if (mydebug&&accept) cout<<" MyObservable::eventcuts event passed 3 "<<endl;
    
    //std::cout<<"TEST: MyObservable::eventcuts: accept: "<<accept<<", mydata->GetNJetCut(): "<<mydata->GetNJetCut()<<std::endl;

    if (accept&&mydata->GetNJetCut()) {
      if (mydebug) cout<<" MyObservable::eventcuts njets= "<<njets<<" mydata->NJet= "<<mydata->GetNJet()<<endl;
        if(njets<=mydata->GetNJet()) return false;
        if (mydebug) if (!accept) cout<<" MyObservable::eventcuts event rejected by Njet cut"<<endl;
    }

    if (mydebug&&accept) cout<<" MyObservable::eventcuts event passed 4 "<<endl;

    for (int  i = 0; i < myevent->GetN(); i++) {
// jet are filled after the lepton
        int id=myevent->GetID(i);
        bool isjet=myevent->IsJet(id);
        if (onejet&&isjet&&oklept&&mydata->DeltaLepJetCut()) {
            double dR=999.;
            if (idjetmax<1) cout<<" MyObservable::eventcuts something is wrong idjetmax= "<<idjetmax<<endl;
            if (idlep<1)    cout<<" MyObservable::eventcuts something is wrong idlep= "<<idlep<<endl;
            if (idlep!=-1) dR=(myevent->Get(i))->delta_R(*myevent->Get(idlep));
            if (mydebug) cout<<i<<" MyObservable::eventcuts dR= "<<dR<<" cut= "<<mydata->GetDeltaRLepJet()<<endl;
            if (dR<mydata->GetDeltaRLepJet()) accept=false;
            if (mydebug) if (!accept) cout<<" MyObservable::eventcuts event cut by dR"<<endl;
        };
    };

    if (mydebug&&accept) cout<<" MyObservable::eventcuts event passed 5 "<<endl;

    if (mydebug||debug)
        if (accept) cout<<" MyObservable::eventcuts event accepted !"<<endl;
        else        cout<<" MyObservable::eventcuts event rejected ! "<<endl;

    return accept;
};
