#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <Rtypes.h>

#include <TMath.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TPaveText.h>
#include <TColor.h>

#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "HttStyles.h"

#endif
double delay1(double fC);
double delay2(double fC);

void getLandauFrac(Float_t tStart, Float_t tEnd, Float_t &sum);

void PulseFraction(Double_t fC, Double_t *TS46);
double Det2(double *b, double *c);
double Det3(double *a, double *b, double *c);

void LandauIntegration() {

  TF1* pulse_shape= new TF1("pulse_shape", "[0]*TMath::Landau((x+[1]),[2],[3])",-50,10000);
  pulse_shape->SetParameter(0,1.0);
  pulse_shape->SetParameter(1,0);
  pulse_shape->SetParameter(2,14.36);
  pulse_shape->SetParameter(3,4.46);

  cout <<  pulse_shape->Integral(0,10000) << endl;
  cout <<  pulse_shape->Integral(0,100) << endl;
  cout << 1/(pulse_shape->Integral(0,10000)/pulse_shape->Integral(0,100)) << endl;
    
  //TCanvas *c1 = MakeCanvas("c1", "c1", 800, 600);

  //TGraph *niFrac = new TGraph(200);
  //TGraph *parFrac = new TGraph(200); 

  Double_t charges[10]={-0.0373423, 1.25702, -0.321804, -1.19814, 11.9322, 3.4333, 2.96879, -0.107635, 2.13894, -0.919258};
  //Double_t charges[10]={-0.118811, -0.381584, 0.705285, -1.96923, 5.44974, 2.96178, -0.417049, -0.859477, -0.118811, -0.381584};
  //Double_t charges[10]={0.283782, 0.0382726, 0.43603, -1.27933, 2.48983, 1.14717, 0.43603, -1.27933, 0.283782, 1.14717};
  //cout << "{";
  //for (Int_t k=-25; k<100; k++) {
    //Float_t tsShift=delay2(k/2)-2;
    //Float_t tsShift=k;
    //cout << "-------" << endl;
    //cout << "charge: " << k << endl;
    //cout << "time slew corr: " << tsShift << endl;
    //Float_t i=0;
    //getLandauFrac(tsShift,25+tsShift,i);
    //cout << i << ", ";
    //Float_t n=0;
    //getLandauFrac(25+tsShift,50+tsShift,n);
    //cout << "intime fraction: " << i << endl;
    //cout << "next fraction:   " << n << endl;
    //cout << k/2 << " " << tsShift << " " << i << endl;
    //niFrac->SetPoint(k,k,i);

    //double TS35[3];
    //PulseFraction(k, TS35);
    //parFrac->SetPoint(k,k,TS35[0]);
  //}
  //cout << "}" << endl;

  /*  niFrac->Draw();
  niFrac->GetYaxis()->SetRangeUser(0.55,0.75);
  //niFrac->GetYaxis()->SetRangeUser(0.1,0.3);
  niFrac->GetXaxis()->SetTitle("Charge");
  niFrac->GetYaxis()->SetTitle("In-time fraction");
  parFrac->SetMarkerColor(kRed);
  parFrac->Draw("same p");*/

  Float_t tsShift3=delay1(charges[3]);
  Float_t tsShift4=delay1(charges[4]);
  Float_t tsShift5=delay1(charges[5]);

  Float_t i3=0;
  getLandauFrac(tsShift3,tsShift3+25,i3);
  Float_t n3=0;
  getLandauFrac(tsShift3+25,tsShift3+50,n3);

  Float_t i4=0;
  getLandauFrac(tsShift4,tsShift4+25,i4);
  Float_t n4=0;
  getLandauFrac(tsShift4+25,tsShift4+50,n4);

  Float_t i5=0;
  getLandauFrac(tsShift5,tsShift5+25,i5);
  Float_t n5=0;
  getLandauFrac(tsShift5+25,tsShift5+50,n5);

  Float_t ch3=charges[3]/i3;
  Float_t ch4=(i3*charges[4]-n3*charges[3])/(i3*i4);
  Float_t ch5=(n3*n4*charges[3]-i3*n4*charges[4]+i3*i4*charges[5])/(i3*i4*i5);

  double TS35[3];
  double TS46[3];
  double TS57[3];
  cout << "call on charge 3" << endl;
  PulseFraction(charges[3], TS35);
  cout << "call on charge 4" << endl;
  PulseFraction(charges[4], TS46);
  cout << "call on charge 5" << endl;
  PulseFraction(charges[5], TS57);
 
  cout << "TS 3 PULSE" << endl;
  cout << i3 << " " << TS35[0] << endl;
  cout << n3 << " " << TS35[1] << endl;
  cout << " 0 " << TS35[2] << endl;

  cout << "TS 4 PULSE" << endl;
  cout << i4 << " " << TS46[0] << endl;
  cout << n4 << " " << TS46[1] << endl;
  cout << " 0 " << TS46[2] << endl;
  
  cout << "TS 5 PULSE" << endl;
  cout << i5 << " " << TS57[0] << endl;
  cout << n5 << " " << TS57[1] << endl;
  cout << " 0 " << TS57[2] << endl;

  double a3[3] = {TS35[0], TS35[1], TS35[2]};
  double b3[3] = {0., TS46[0], TS46[1]};
  double c3[3] = {0., 0., TS57[0]};
  double d3[3] = {charges[3], charges[4], charges[5]};

  double deno3 = Det3(a3, b3, c3);

  double A3 = Det3(d3, b3, c3) / deno3;
  double A4 = Det3(a3, d3, c3) / deno3;
  double A5 = Det3(a3, b3, d3) / deno3;

  cout << "SOLUTION" << endl;
  cout << ch3 << " " << A3 << endl;
  cout << ch4 << " " << A4 << endl;
  cout << ch5 << " " << A5 << endl;

}

double delay1(double fC) {
  //double rawDelay=13.307784-1.556668*log(fC);
  //return (rawDelay<0)?(0):((rawDelay>10)?(10):(rawDelay));
  return 13.7436-1.57*log(fC+28);
}

double delay2(double fC) {
  return 9.453-1.948*log(fC+88.18);
}

// Landau function integrated in 1 ns intervals
//Landau pulse shape from https://indico.cern.ch/event/345283/contribution/3/material/slides/0.pdf
//Landau turn on by default at left edge of time slice 
// normalized to 1 on [0,10000]
void getLandauFrac(Float_t tStart, Float_t tEnd, Float_t &sum) {

  Float_t landauFrac[125] = {0, 7.6377e-05, 0.000418655, 0.00153692, 0.00436844, 0.0102076, 0.0204177, 0.0360559, 0.057596, 0.0848493, 0.117069, 0.153152, 0.191858, 0.23198, 0.272461, 0.312438, 0.351262, 0.388476, 0.423788, 0.457036, 0.488159, 0.517167, 0.54412, 0.569112, 0.592254, 0.613668, 0.633402, 0.651391, 0.667242, 0.680131, 0.688868, 0.692188, 0.689122, 0.67928, 0.662924, 0.64087, 0.614282, 0.584457, 0.552651, 0.51997, 0.487317, 0.455378, 0.424647, 0.395445, 0.367963, 0.342288, 0.318433, 0.29636, 0.275994, 0.257243, 0.24, 0.224155, 0.2096, 0.196227, 0.183937, 0.172635, 0.162232, 0.15265, 0.143813, 0.135656, 0.128117, 0.12114, 0.114677, 0.108681, 0.103113, 0.0979354, 0.0931145, 0.0886206, 0.0844264, 0.0805074, 0.0768411, 0.0734075, 0.0701881, 0.0671664, 0.0643271, 0.0616564, 0.0591418, 0.0567718, 0.054536, 0.0524247, 0.0504292, 0.0485414, 0.046754, 0.0450602, 0.0434538, 0.041929, 0.0404806, 0.0391037, 0.0377937, 0.0365465, 0.0353583, 0.0342255, 0.0331447, 0.032113, 0.0311274, 0.0301854, 0.0292843, 0.0284221, 0.0275964, 0.0268053, 0.0253052, 0.0238536, 0.0224483, 0.0210872, 0.0197684, 0.0184899, 0.01725, 0.0160471, 0.0148795, 0.0137457, 0.0126445, 0.0115743, 0.0105341, 0.00952249, 0.00853844, 0.00758086, 0.00664871, 0.00574103, 0.00485689, 0.00399541, 0.00315576, 0.00233713, 0.00153878, 0.000759962, 0 };

  // can be further optimized to reduce computational time
  if (abs(tStart-tEnd)!=25) {
    sum=0;
    return;
  }
  sum=landauFrac[int(ceil(tStart+25))];
  return;

}

void PulseFraction(Double_t fC, Double_t *TS46) {
  
  //static Double_t TS3par[3] = {0.44, -18.6, 5.136}; //Gaussian parameters: norm, mean, sigma for the TS3 fraction  
  static Double_t TS4par[3] = {0.71, -5.17, 12.23}; //Gaussian parameters: norm, mean, sigma for the TS4 fraction
  static Double_t TS5par[3] = {0.258, 0.0178, 4.786e-4}; // pol2 parameters for the TS5 fraction
  static Double_t TS6par[4] = {0.06391, 0.002737, 8.396e-05, 1.475e-06};// pol3 parameters for the TS6 fraction
  
  Double_t tslew = -delay1(fC);
  cout << "time slew: " << tslew << endl;
  
  TS46[0] = TS4par[0] * TMath::Gaus(tslew,TS4par[1],TS4par[2]); // fraction of pulse in the TS4
  cout << "intime fraction: " << TS46[0] << endl;
  TS46[1] = TS5par[0] + TS5par[1]*tslew + TS5par[2]*tslew*tslew; // fraction of pulse in the T5S
  cout << "next fraction:   " << TS46[1] << endl;
  TS46[2] = TS6par[0] + TS6par[1]*tslew + TS6par[2]*tslew*tslew + TS6par[3]*tslew*tslew*tslew; //fraction of pulse in the TS6

  cout << fC << " " << tslew << " " << TS46[0] << endl;
  
  return;
}

double Det2(double *b, double *c){
  return b[1]*c[2]-b[2]*c[1];
}

double Det3(double *a, double *b, double *c){
  return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]);
}
