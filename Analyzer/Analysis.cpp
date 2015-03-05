#define Analysis_cxx

#include "Analysis.h"
#include "readparameters/readparameters.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>
#include "TLine.h"
#include "TLegend.h"

using namespace std;

double round_nplaces(double value, int to){
  int places = 1, whole = value;
  for(int i = 0; i < to; i++) places *= 10;
  value -= whole; //leave decimals
  value *= places; //0.1234 -> 123.4
  value = round(value);//123.4 -> 123
  value /= places; //123 -> .123
  value += whole; //bring the whole value back
  return value;
}

string int2string(int i){
  stringstream ss;
  string ret;
  ss<<i;
  ss>>ret;
  return ret;
}

// Ugly hack to apply energy corrections to some HB- cells
double eCorr(int ieta, int iphi, double energy) {
// return energy correction factor for HBM channels 
// iphi=6 ieta=(-1,-15) and iphi=32 ieta=(-1,-7)
// I.Vodopianov 28 Feb. 2011
  static const float low32[7]  = {0.741,0.721,0.730,0.698,0.708,0.751,0.861};
  static const float high32[7] = {0.973,0.925,0.900,0.897,0.950,0.935,1};
  static const float low6[15]  = {0.635,0.623,0.670,0.633,0.644,0.648,0.600,
				  0.570,0.595,0.554,0.505,0.513,0.515,0.561,0.579};
  static const float high6[15] = {0.875,0.937,0.942,0.900,0.922,0.925,0.901,
				  0.850,0.852,0.818,0.731,0.717,0.782,0.853,0.778};
  
  double slope, mid, en;
  double corr = 1.0;

  if (!(iphi==6 && ieta<0 && ieta>-16) && !(iphi==32 && ieta<0 && ieta>-8)) 
    return corr;

  int jeta = -ieta-1;
  double xeta = (double) ieta;
  if (energy > 0.) en=energy;
  else en = 0.;

  if (iphi == 32) {
    slope = 0.2272;
    mid = 17.14 + 0.7147*xeta;
    if (en > 100.) corr = high32[jeta];
    else corr = low32[jeta]+(high32[jeta]-low32[jeta])/(1.0+exp(-(en-mid)*slope));
  }
  else if (iphi == 6) {
    slope = 0.1956;
    mid = 15.96 + 0.3075*xeta;
    if (en > 100.0) corr = high6[jeta];
    else corr = low6[jeta]+(high6[jeta]-low6[jeta])/(1.0+exp(-(en-mid)*slope));
  }

  return corr;
}

int main(int argc, char **argv)
{
  int ret=0;
  if (argc!=2) {
    cerr<<"Usage: ./Analysis <paramfile>"<<endl;
    ret=1;
  } else {

    readparameters rp(argv[1]);
    //TChain* ch = new TChain("HcalNoiseTree");
    TChain* ch = new TChain("ExportTree/HcalTree");

    string filelistname;

    filelistname=rp.get<string>((string("in_filelist")).c_str());
    string line;
    ifstream filelist(filelistname.c_str());
    if (filelist.fail()) { //catch
      cerr << "\nERROR: Could not open " << filelistname << endl;
      exit(1);
    }
    while (getline(filelist,line)) {
      ch->Add(line.c_str());
    }

  Analysis Ana25ns((TTree*) ch);

  Ana25ns.Init(argv[1]);
  Ana25ns.Process();
  Ana25ns.Finish();
  }
 return ret;
}

Analysis::~Analysis() {
}
Analysis::Analysis(TTree *tree):analysistree(tree){};

void Analysis::Init(char* paramfile)
{
  try {
    readparameters rp(paramfile);
    try {Output_File=rp.get<string>("Output_File");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Entries=rp.get<int>("Entries");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Plot_Dir=rp.get<string>("Plot_Dir");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Region=rp.get<int>("Region");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Condition=rp.get<int>("Condition");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Baseline=rp.get<int>("Baseline");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Time_Slew=rp.get<int>("Time_Slew");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Neg_Charges=rp.get<int>("Neg_Charges");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Threshold=rp.get<float>("Threshold");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
    try {Quantile=rp.get<float>("Quantile");}
    catch (exception& e) {cerr<<e.what()<<endl;} 
  } 
  catch (exception& e) {cerr<<e.what()<<endl;} 

  cout << "Running on ";
  if (Entries==-1) cout << "all events." << endl;
  else cout << "first " << Entries << " events. " << endl;

  cout << "Output ROOT file: " << Output_File << endl;

  cout << "Using channels in ";
  if (Region==All) { cout << "all HCAL regions. " << endl; }
  else if (Region==Barrel) { cout << "HCAL barrel. " << endl; }
  else if (Region==Endcap) { cout << "HCAL endcap. " << endl; }

  if (Condition==0) {
    cout << "With no PU. " << endl;
  }
  else if (Condition==50) {
    cout << "50 ns spacing, 20 PU." << endl;
  }
  else if (Condition==25) {
    cout << "25 ns spacing, 20 PU." << endl;
  }
  else {
    Condition=0;
    cout << "Unrecognized run condition, using no PU." << endl;
  }

  int check=mkdir(Plot_Dir.c_str(),755);
  if (!check) {
    cout << "Saving files to: " << Plot_Dir << endl;
  }
  else {
    cout << "Double check your plot directory exists! Something funny happened." << endl;
    //exit(1);
  }

  if (Baseline==PedestalSub::DoNothing) {
    cout << "Pedestal subtraction only." << endl;
  }
  else if (Baseline==PedestalSub::AvgWithThresh) {
    cout << "Pedestal subtraction, average baseline subtraction with threshold: " << Threshold << endl;
  }
  else if (Baseline==PedestalSub::AvgWithoutThresh) {
    cout << "Pedestal subtraction, average baseline subtraction with no threshold. " << endl;
  }
  else if (Baseline==PedestalSub::AvgWithThreshNoPedSub) {
    cout << "Average baseline+pedestal subtraction with threshold: " << Threshold << endl;
  }
  else if (Baseline==PedestalSub::Percentile) {
    cout << "Percentile-based pedestal subtraction ";
    if (Quantile<0 || Quantile>1) {
      cout << endl << "Quantile value out of range. Not running." << endl;
      exit(1);
    }
    else  {
      cout << "with quantile value: " << Quantile << endl;
    }
  }

  if (Time_Slew==HcalTimeSlew::TestStand) cout << "Using test stand medium WP time slew parameterization." << endl;
  else cout << "Sorry, I don't know which time slew correction you asked for." << endl;

  if (Neg_Charges==HLTv2::DoNothing) cout << "Not requiring positive charge outputs." << endl;
  else cout << "Requiring positive charge outputs." << endl;

  return; 
}

void Analysis::Process() {
  if (fChain == 0) return;

  if (Entries==-1) Entries=fChain->GetEntries();

  fout = new TFile(Output_File.c_str(), "RECREATE");

  //DeriveTimeslew();
  DoHlt();
}

void Analysis::DeriveTimeslew() {

  //Get Pulse Fractions

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  gStyle->SetOptStat(0);
  
  TF1* pulse_shape= new TF1("pulse_shape", "[0]*TMath::Landau((x+[1]),[2],[3])",-50,10000);
  pulse_shape->SetParameter(0,1.0/4.45795);
  pulse_shape->SetParameter(1,0);
  pulse_shape->SetParameter(2,14.36);
  pulse_shape->SetParameter(3,4.46);
  
  TGraph *grTS_P = new TGraph(150);
  TGraph *grTS_IT = new TGraph(150);
  TGraph *grTS_N = new TGraph(150);
  TGraph *grTS_NN = new TGraph(150);

  TGraph *grTS_NIT = new TGraph(150);
  TGraph *grNIT_TS = new TGraph(150);

  for (Int_t i=0; i<150; i++) {
    Float_t slew = Float_t(i)/10-5;
    grTS_P->SetPoint(i,slew,pulse_shape->Integral(-25-slew,-slew));
    grTS_IT->SetPoint(i,slew,pulse_shape->Integral(-slew,25-slew));
    grTS_N->SetPoint(i,slew,pulse_shape->Integral(25-slew,50-slew));
    grTS_NN->SetPoint(i,slew,pulse_shape->Integral(50-slew,75-slew));

    grTS_NIT->SetPoint(i,slew,pulse_shape->Integral(25-slew,50-slew)/pulse_shape->Integral(-slew,25-slew));
    grNIT_TS->SetPoint(i,pulse_shape->Integral(25-slew,50-slew)/pulse_shape->Integral(-slew,25-slew),slew);
  }

  grTS_P->GetXaxis()->SetTitle("TIME SLEW");
  grTS_P->GetYaxis()->SetTitle("fraction in PREVIOUS time slice");
  grTS_P->SetTitle("From Landau fxn");
  grTS_P->SetLineWidth(3);
  grTS_P->Draw("al");
  c1->SaveAs("ts_prev.png");

  grTS_IT->GetXaxis()->SetTitle("TIME SLEW");
  grTS_IT->GetYaxis()->SetTitle("fraction in INTIME time slice");
  grTS_IT->SetTitle("From Landau fxn");
  grTS_IT->SetLineWidth(3);
  grTS_IT->Draw("al");
  c1->SaveAs("ts_intime.png");

  grTS_N->GetXaxis()->SetTitle("TIME SLEW");
  grTS_N->GetYaxis()->SetTitle("fraction in NEXT time slice");
  grTS_N->SetTitle("From Landau fxn");
  grTS_N->SetLineWidth(3);
  grTS_N->Draw("al");
  c1->SaveAs("ts_next.png");

  grTS_NN->GetXaxis()->SetTitle("TIME SLEW");
  grTS_NN->GetYaxis()->SetTitle("fraction in NEXT-TO-NEXT time slice");
  grTS_NN->SetTitle("From Landau fxn");
  grTS_NN->SetLineWidth(3);
  grTS_NN->Draw("al");
  c1->SaveAs("ts_nexttonext.png");

  grTS_NIT->GetXaxis()->SetTitle("TIME SLEW");
  grTS_NIT->GetYaxis()->SetTitle("NEXT/INTIME");
  grTS_NIT->SetTitle("From Landau fxn");
  grTS_NIT->SetLineWidth(3);
  grTS_NIT->Draw("al");
  c1->SaveAs("ts_nit.png");

  grNIT_TS->GetXaxis()->SetTitle("NEXT/INTIME");
  grNIT_TS->GetYaxis()->SetTitle("TIME SLEW");
  grNIT_TS->SetTitle("From Landau fxn");
  grNIT_TS->SetLineWidth(3);
  grNIT_TS->Draw("al");
  c1->SaveAs("nit_ts.png");

  TProfile *pIT_NIT = new TProfile("pIT_NIT","",50,0,500);
  TProfile *pIT_TS = new TProfile("pIT_TS","",50,0,500);

  for (Int_t i=0; i<Entries; i++) {
    fChain->GetEntry(i);
    for (Int_t j=0; j<PulseCount;j++) {
      Float_t it=Charge[j][4]+Pedestal[j][4];
      Float_t nt=Charge[j][5]+Pedestal[j][5];
      if (it!=0) {
        pIT_NIT->Fill(it,nt/it);
        pIT_TS->Fill(it, grNIT_TS->Eval(nt/it));
      }
    }
  }

  pIT_NIT->GetXaxis()->SetTitle("INTIME charge");
  pIT_NIT->GetYaxis()->SetTitle("NEXT/INTIME");
  pIT_NIT->SetTitle("From MC");
  pIT_NIT->Draw("");
  c1->SaveAs("intime_nit.png");

  TF1 *fitTS = new TF1("fitTS", "[0]+[1]*TMath::Log(x)",10,1000);
  fitTS->SetLineColor(kRed);

  pIT_TS->Fit("fitTS");

  pIT_TS->GetXaxis()->SetTitle("INTIME charge");
  pIT_TS->GetYaxis()->SetTitle("TIME SLEW");
  pIT_TS->SetTitle("From MC + Landau fxn");
  pIT_TS->Draw("");
  fitTS->Draw("same");
  c1->SaveAs("intime_ts.png");

}

void Analysis::DoHlt() {

  //=========================================================================      
  // These are the values we should set for Method 2 config with the python files
  // Don't currently have a setup to input with python but we shouldn't have to
  // change these values for our tests
  
  // --------------------------------------------------------------------
  bool iPedestalConstraint = true;
  bool iTimeConstraint = true;
  bool iAddPulseJitter = false;
  bool iUnConstrainedFit = false;
  bool iApplyTimeSlew = true;
  double iTS4Min = 5.;
  double iTS4Max = 500.;
  double iPulseJitter = 1.;
  double iTimeMean = -5.5;
  double iTimeSig = 5.;
  double iPedMean = 0.;
  double iPedSig = 0.5;
  double iNoise = 1.;
  double iTMin = -18.;
  double iTMax = 7.;
  double its3Chi2 = 5.;
  double its4Chi2 = 15.;
  double its345Chi2 = 100.;
  double iChargeThreshold = 6.;
  int iFitTimes = 1;
    
  //========================================================================
  // Set the Method 2 Parameters here
  psFitOOTpuCorr_->setPUParams(iPedestalConstraint,iTimeConstraint,iAddPulseJitter,iUnConstrainedFit,
			       iApplyTimeSlew,iTS4Min, iTS4Max, iPulseJitter,iTimeMean,iTimeSig,
			       iPedMean,iPedSig,iNoise,iTMin,iTMax,its3Chi2,its4Chi2,its345Chi2,
			       iChargeThreshold,HcalTimeSlew::Medium, iFitTimes);
  
  // Now set the Pulse shape type
  psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(105));
  
  //Setup HLT pedestal/baseline subtraction module
  //pedSubFxn_->Init(((PedestalSub::Method)Baseline), Condition, Threshold, Quantile);
  pedSubFxn_->Init(((PedestalSub::Method)1), Condition, 2.7, 0.0);

  //Set HLT module
  //hltv2_->Init((HcalTimeSlew::ParaSource)Time_Slew, HcalTimeSlew::Medium, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  //hltv2_->Init((HcalTimeSlew::ParaSource)2, HcalTimeSlew::Medium, (HLTv2::NegStrategy)0, *pedSubFxn_); // no neg. correction
  //hltv2_->Init((HcalTimeSlew::ParaSource)2, HcalTimeSlew::Medium, (HLTv2::NegStrategy)1, *pedSubFxn_); // "KISS" correction
  hltv2_->Init((HcalTimeSlew::ParaSource)2, HcalTimeSlew::Medium, (HLTv2::NegStrategy)2, *pedSubFxn_); // Greg's correction

  //Setup plots for what we care about
  int xBins=200, xMin=0,xMax=200;

  TH1D *a3 = new TH1D("a3","", xBins,xMin,xMax);
  TH1D *a4 = new TH1D("a4","", xBins,xMin,xMax);
  TH1D *a5 = new TH1D("a5","", xBins,xMin,xMax);

  TH2D *a4v3 = new TH2D("a4v3","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a4v5 = new TH2D("a4v5","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a5v3 = new TH2D("a5v3","", xBins,xMin,xMax,xBins,xMin,xMax);

  TH2D* h45vHLT = new TH2D("h45vHLT", "", xBins,xMin,xMax,xBins,xMin,xMax);
  TProfile* p45vHLT = new TProfile("p45vHLT", "", xBins,xMin,xMax,-10,10);

  TH2D* hM2vHLT = new TH2D("hM2vHLT", "", xBins,xMin,xMax,xBins,xMin,xMax);
  TProfile *pM2vHLT = new TProfile("pM2vHLT", "", xBins,xMin,xMax,-10,10);

  //Loop over all events
  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
    cout << jentry << endl;
    for (int j = 0; j < (int)PulseCount; j++) {
      if (IEta[j]>16 && Region==Barrel) continue;
      if (IEta[j]<17 && Region==Endcap) continue;

      std::vector<double> inputCaloSample, inputPedestal, inputGain;
      std::vector<double> offlineAns, hltAns;

      for (int i=0; i<10; i++) {
	inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
	inputPedestal.push_back(Pedestal[j][i]);
	inputGain.push_back(Gain[j][i]);
      }
      
      // Begin Method 2
      psFitOOTpuCorr_->apply(inputCaloSample,inputPedestal,inputGain,offlineAns);
      
      // Begin Online
      hltv2_->apply(inputCaloSample,inputPedestal,hltAns);
      //hltv2_->applyXM(inputCaloSample,inputPedestal,hltAns);

      if (hltAns.size()>1) {

	//Fill Histograms
	a3->Fill(hltAns.at(0));
	a4->Fill(hltAns.at(1));
	a5->Fill(hltAns.at(2));
	
	a4v3->Fill(hltAns.at(1), hltAns.at(0));
	a4v5->Fill(hltAns.at(1), hltAns.at(2));	
	a5v3->Fill(hltAns.at(2), hltAns.at(0));

	h45vHLT->Fill( (Charge[j][4]), hltAns.at(1),1);
	p45vHLT->Fill((Charge[j][4]), -(hltAns.at(1)-(Charge[j][4]))/((Charge[j][4])),1);

	if (offlineAns.size()>1) {
	  hM2vHLT->Fill( offlineAns.at(0)/Gain[j][0], hltAns.at(1),1);
	  pM2vHLT->Fill( offlineAns.at(0)/Gain[j][0], -(hltAns.at(1)*Gain[j][0]-(offlineAns.at(0)))/((offlineAns.at(0))),1);
	}
      }

    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  gStyle->SetOptStat(0);

  a3->GetXaxis()->SetTitle("A3 [fC]");
  a3->GetXaxis()->SetTitleSize(0.05);
  a3->GetYaxis()->SetTitle("Counts");
  a3->GetYaxis()->SetTitleSize(0.05);
  a3->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a3.png");
  
  a4->GetXaxis()->SetTitle("A4 [fC]");
  a4->GetXaxis()->SetTitleSize(0.05);
  a4->GetYaxis()->SetTitle("Counts");
  a4->GetYaxis()->SetTitleSize(0.05);
  a4->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4.png");
  
  a5->GetXaxis()->SetTitle("A5 [fC]");
  a5->GetXaxis()->SetTitleSize(0.05);
  a5->GetYaxis()->SetTitle("Counts");
  a5->GetYaxis()->SetTitleSize(0.05);
  a5->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5.png");

  a4v3->GetXaxis()->SetTitle("A4 [fC]");
  a4v3->GetXaxis()->SetTitleSize(0.05);
  a4v3->GetYaxis()->SetTitle("A3 [fC]");
  a4v3->GetYaxis()->SetTitleSize(0.05);
  a4v3->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v3.png");

  a4v5->GetXaxis()->SetTitle("A4 [fC]");
  a4v5->GetXaxis()->SetTitleSize(0.05);
  a4v5->GetYaxis()->SetTitle("A5 [fC]");
  a4v5->GetYaxis()->SetTitleSize(0.05);
  a4v5->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v5.png");

  a5v3->GetXaxis()->SetTitle("A5 [fC]");
  a5v3->GetXaxis()->SetTitleSize(0.05);
  a5v3->GetYaxis()->SetTitle("A3 [fC]");
  a5v3->GetYaxis()->SetTitleSize(0.05);
  a5v3->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5v3.png");

  h45vHLT->GetXaxis()->SetTitle("Charge TS4 [fC]");
  h45vHLT->GetXaxis()->SetTitleSize(0.05);
  h45vHLT->GetYaxis()->SetTitle("Charge from HLT [fC]");
  h45vHLT->GetYaxis()->SetTitleSize(0.05);
  h45vHLT->Draw("");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/h4vHLT.png");
  
  p45vHLT->GetXaxis()->SetTitle("Charge in TS4 [fC]");
  p45vHLT->GetXaxis()->SetTitleSize(0.05);
  p45vHLT->GetYaxis()->SetTitle("(HLT - E4)/E4");
  p45vHLT->GetYaxis()->SetTitleSize(0.05);
  p45vHLT->GetYaxis()->SetRangeUser(-1,1);
  p45vHLT->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/p4vHLT.png");

  hM2vHLT->GetXaxis()->SetTitle("M2 Charge [fC]");
  hM2vHLT->GetXaxis()->SetTitleSize(0.05);
  hM2vHLT->GetYaxis()->SetTitle("Charge from HLT [fC]");
  hM2vHLT->GetYaxis()->SetTitleSize(0.05);
  hM2vHLT->Draw("");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hM2vHLT.png");

  pM2vHLT->GetXaxis()->SetTitle("M2 Charge [fC]");
  pM2vHLT->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetTitle("(HLT - M2)/(M2)");
  pM2vHLT->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetRangeUser(-1,1);
  pM2vHLT->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/pM2vHLT.png");

}

void Analysis::Finish()
{

  fout->cd();
  fout->Write();
  fout->Close();
}
