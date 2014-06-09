#ifndef MyHisto_H
#define MyHisto_H
#include <iostream>
using namespace std;
#include <cmath>
#include <string>
#include "TH1D.h"
/******************************************************************
 ** Method Implementations
 ******************************************************************/

class MyHisto {
 private:

  bool debug; 
  bool newevent;

  int nev;
  int nevall;

  double wev;
  //double wevsq;
  double wsum;
  //double wsumsq;
  double wsum2;

  TH1D *hsum;
  TH1D *hsum2;
  TH1D *htmp;
  TH1D *htmp2;


public:
  MyHisto(string name, string htit, int nbin, double xmin, double xmax);
  MyHisto(TH1D *h1);
  
  void Initialise();
  void Fill(double obs, double w, bool newevent);
  void AddLastEvent();

  double GetSum(){return wsum/nev;};
  double GetSumError();
  double GetSumSquare(){ return wsum2/nev;};
  int GetNumberEvents(){ return nev;};
  int GetNewEvent(){ return newevent;};

  void Print();
  TH1D* GetHisto(){ return hsum;};

};


#endif

