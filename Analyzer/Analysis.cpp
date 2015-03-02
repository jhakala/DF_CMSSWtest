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
  Ana25ns.DefineHistograms();
  Ana25ns.TSP();
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

  DoHlt();
  MakeTimeSlewPlots();
  //MakePedestalPlots();
  
}

void Analysis::DefineHistograms()
{
}
 
void Analysis::TSP()
{
  TimeSlewPulse_All = new TH2F("TimeSlewPulse_All","Time Slew [ns] vs Charge in TS4 [fC]",25,-14.5,10.5,100,0.,1000.);
  TimeSlewPulse_All->Sumw2();
  TimeSlewPulse_HB = new TH2F("TimeSlewPulse_HB","Time Slew [ns] vs Charge in TS4 [fC] in HB",25,-14.5,10.5,100,0.,1000.);
  TimeSlewPulse_HB->Sumw2();
  TimeSlewPulse_HE = new TH2F("TimeSlewPulse_HE","Time Slew [ns] vs Charge in TS4 [fC] in HE",25,-14.5,10.5,100,0.,1000.);
  TimeSlewPulse_HE->Sumw2();

  //Logrithmic
  timeslewFit = new TF1("timeslewFit", "[0]+[1]*TMath::Log(x+[2])",0.,1000.);
  timeslewFit->SetParameters(9.27591e+02,-9.71576e+01, 1.39366e+04);

  //Linear
//  timeslewFit = new TF1("timeslewFit", "[0]+[1]*x",0.,1000.);
//  timeslewFit->SetParameters(4.61732e-01,-6.77490e-03);

  slewFit = new TF1("slewFit","pol4*expo(5)",-10.,14.);
  slewFit->SetParameters(1.07618e-02,-4.19145e-06,2.70310e-05,-8.71584e-08,1.86597e-07,3.59216e+00,-1.02057e-01);

  fTSP = new TFile("TSP.root", "RECREATE");
  pedSubFxn_->Init(((PedestalSub::Method)Baseline), Condition, Threshold, Quantile);

  if (Entries==-1) Entries=fChain->GetEntries();
  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
    for (int j = 0; j < (int)PulseCount; j++) {
      std::vector<double> inputCaloSample, inputPedestal, inputGain;
      for (int i=0; i<10; i++) {
        inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
        inputPedestal.push_back(Pedestal[j][i]);
        inputGain.push_back(Gain[j][i]);
      }

    double RatioTS54 = 0, TimeSlew = 0., Pulse = 0.;
    TimeSlewParameters->getParameters(inputCaloSample, inputPedestal, RatioTS54, TimeSlew, Pulse, slewFit, *pedSubFxn_);
    TimeSlewPulse_All->Fill(TimeSlew, Pulse);
    if(IEta[j] < 16) TimeSlewPulse_HB->Fill(TimeSlew, Pulse);
    if(IEta[j] >= 16) TimeSlewPulse_HE->Fill(TimeSlew, Pulse);
    }
  }
  gStyle->SetOptFit(1);
  TimeSlewPulse_All->Draw("BOX");
  TimeSlewPulse_All->ProfileY("Fit Y Profile",1,-1,"dos")->Fit("timeslewFit");
  timeslewFit->Write();
  TimeSlewPulse_HB->ProfileY("HBFit Y Profile",1,-1,"dos")->Fit("timeslewFit");
  TimeSlewPulse_HE->ProfileY("HEFit Y Profile",1,-1,"dos")->Fit("timeslewFit");

  fTSP->cd();
  fTSP->Write();
  fTSP->Close();

}

void Analysis::MakeCutflow() 
{
}

void Analysis::FillHistograms()
{
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
  pedSubFxn_->Init(((PedestalSub::Method)Baseline), Condition, Threshold, Quantile);

  //Set HLT module
  hltv2_->Init((HcalTimeSlew::ParaSource)Time_Slew, HcalTimeSlew::Medium, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);

  //Setup plots for what we care about
  int xBins=200, xMin=-10,xMax=40;

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
  h45vHLT->Draw("colz");
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
  hM2vHLT->Draw("colz");
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hM2vHLT.png");

  pM2vHLT->GetXaxis()->SetTitle("M2 Charge [fC]");
  pM2vHLT->GetXaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetTitle("(HLT - M2)/(M2)");
  pM2vHLT->GetYaxis()->SetTitleSize(0.05);
  pM2vHLT->GetYaxis()->SetRangeUser(-1,1);
  pM2vHLT->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/pM2vHLT.png");

}

void Analysis::MakePedestalPlots() {
  /*
  double QUANTILE=0.25;

  TH1D *hist1 = new TH1D("hist1", "", 100,-5,45);
  TH1D *hist2 = new TH1D("hist2", "", 100,-5,45);

  double mean1=0, mean2=0, rms1=0, rms2=0;
  double q1=0, q2=0;

  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
    for (int j = 0; j < (int)PulseCount; j++) {
      std::vector<double> tenC;
      std::vector<double> sevenC;  
      for (int i=0; i<10; i++) {
	mean1+=Charge[j][i];
	rms1+=Charge[j][i]*Charge[j][i];
	tenC.push_back(Charge[j][i]);
	hist1->Fill(Charge[j][i]);
	if (i==4||i==5||i==6) continue;
	mean2+=Charge[j][i];
	rms2+=Charge[j][i]*Charge[j][i];
	sevenC.push_back(Charge[j][i]);
	hist2->Fill(Charge[j][i]);
      }
      //if (j==100) {
      //for (int i=0; i<10; i++) {
      //  cout << tenC[i] << ", ";
      //}
      //cout << "and 0.2 quantile is " << sampleQuantile<10>(&tenC[0],0.2) << endl;
      //}
      q1+=sampleQuantile<10>(&tenC[0],QUANTILE);
      q2+=sampleQuantile<7>(&sevenC[0],QUANTILE);
    }
  }

  mean1=mean1/(Entries*PulseCount*10);
  mean2=mean2/(Entries*PulseCount*7);
  rms1=TMath::Sqrt(rms1/(Entries*PulseCount*10)-mean1*mean1);
  rms2=TMath::Sqrt(rms2/(Entries*PulseCount*7)-mean2*mean2);
  q1=q1/(Entries*PulseCount);
  q2=q2/(Entries*PulseCount);

  cout << -inverseGaussCDF(QUANTILE) << endl;

  cout << "no time slices excluded: " << mean1 << " +/- " << rms1 << endl;
  cout << "excluding 4,5,6:         " << mean2 << " +/- " << rms2 << endl;
  cout << "20th quantile of 10 TS:  " << q1 << endl;
  cout << "20th quantile of 7 TS:   " << q2 << endl;
  cout << "Would subtract:          " << q2-inverseGaussCDF(QUANTILE)*rms2 << endl;
  cout << "or.... Would subtract:   " << q1-inverseGaussCDF(QUANTILE)*rms1 << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  hist1->SetLineWidth(3);
  hist2->SetLineWidth(3);
  hist1->GetXaxis()->SetTitle("Pedestal-Subtracted Charge [fC]");
  hist1->GetYaxis()->SetTitle("Counts");
  hist1->GetYaxis()->SetRangeUser(1,1e6);
  hist1->Draw();

  hist2->SetLineColor(kRed);
  hist2->Draw("same");
  
  TLine *line = new TLine(q2-inverseGaussCDF(QUANTILE)*rms2, 1, q2-inverseGaussCDF(QUANTILE)*rms2, 1e6);
  line->SetLineWidth(3);
  line->SetLineColor(kBlue);
  line->Draw();
  char fname[50];
  sprintf(fname, "quantile_%.2f_%i.png", QUANTILE,Condition);

  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(hist1, "All TS", "l");
  leg->AddEntry(hist2, "Not TS4/5/6", "l");
  leg->AddEntry(line, "Avg. BSE", "l");
  leg->Draw();

  c1->SaveAs(fname);
  */

  vector<PedestalSub*> vPedSub;
  vector<TH1D*> vCorrHist;
  
  char hname[50];

  int nMethod=8;

  for (int i=0; i<nMethod; i++) {
    vPedSub.push_back(new PedestalSub);
    sprintf(hname, "pedMethod_%i",i);
    vCorrHist.push_back(new TH1D(hname, "", 60, -2,4));
  }
  
  vPedSub[0]->Init(PedestalSub::Percentile, Condition, 0.0, 0.25);
  vCorrHist[0]->GetXaxis()->SetTitle("PedestalSub::Percentile 0.25 [fC]");

  vPedSub[1]->Init(PedestalSub::AvgWithThresh, Condition, 1.7, 0.0);
  vCorrHist[1]->GetXaxis()->SetTitle("PedestalSub::AvgWithThresh 1.7 [fC]");

  vPedSub[2]->Init(PedestalSub::AvgWithThresh, Condition, 2.7, 0.0);
  vCorrHist[2]->GetXaxis()->SetTitle("Baseline estimate <Min[2.7,Q]> [fC]");

  vPedSub[3]->Init(PedestalSub::AvgWithThreshNoPedSub, Condition, 5.0, 0.0);
  vCorrHist[3]->GetXaxis()->SetTitle("PedestalSub::AvgWithThreshNoPedSub 5.0 [fC]");

  vPedSub[4]->Init(PedestalSub::AvgWithThreshNoPedSub, Condition, 6.0, 0.0);
  vCorrHist[4]->GetXaxis()->SetTitle("PedestalSub::AvgwithThreshNoPedSub 6.0 [fC]");

  vPedSub[5]->Init(PedestalSub::Percentile, Condition, 0.0, 0.2);
  vCorrHist[5]->GetXaxis()->SetTitle("PedestalSub::Percentile 0.2 [fC]");

  vPedSub[6]->Init(PedestalSub::Percentile, Condition, 0.0, 0.3);
  vCorrHist[6]->GetXaxis()->SetTitle("PedestalSub::Percentile 0.3 [fC]");

  vPedSub[7]->Init(PedestalSub::AvgWithoutThresh, Condition, 0.0, 0.0);
  vCorrHist[7]->GetXaxis()->SetTitle("PedestalSub::AvgWithoutThresh [fC]");

  
  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
    for (int j = 0; j < (int)PulseCount; j++) {
      std::vector<double> inputCaloSample, inputPedestal;
      for (int i=0; i<10; i++) {
	inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
	inputPedestal.push_back(Pedestal[j][i]);
      }
      for (int k=0; k<nMethod; k++) {
	vCorrHist[k]->Fill(vPedSub[k]->GetCorrection(inputCaloSample, inputPedestal));
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(0);

  char fname[50];

  for (int i=0; i<nMethod; i++) {
    vCorrHist[i]->GetYaxis()->SetTitle("Counts");
    vCorrHist[i]->Draw("hist");
    cout << "Mean for distribution " << i << ": " << vCorrHist[i]->GetMean() << endl;
    cout << "RMS for distribution " << i << ": " << vCorrHist[i]->GetRMS() << endl;
    sprintf(fname, "ped_plots/subtract_%i_%i.png", i, Condition);
    c1->SaveAs(fname);
  }
  
}

void Analysis::MakeTimeSlewPlots() {
  
  //Setup HLT pedestal/baseline subtraction module
  //pedSubFxn_->Init(((PedestalSub::Method)Baseline), Condition, Threshold, Quantile);
  pedSubFxn_->Init(((PedestalSub::Method)1), Condition, 2.7, 0.0);

  //Set HLT module
  HLTv2* hlt_wNeg= new HLTv2;
  HLTv2* hlt_noNeg= new HLTv2;
  hlt_wNeg->Init(HcalTimeSlew::MC, HcalTimeSlew::Medium, (HLTv2::NegStrategy)0, *pedSubFxn_);
  hlt_noNeg->Init(HcalTimeSlew::MC, HcalTimeSlew::Medium, (HLTv2::NegStrategy)1, *pedSubFxn_);

  int xBins=100, xMin=-10,xMax=40;

  TH1D *a3j = new TH1D("a3j","", xBins,xMin,xMax);
  TH1D *a4j = new TH1D("a4j","", xBins,xMin,xMax);
  TH1D *a5j = new TH1D("a5j","", xBins,xMin,xMax);

  TH1D *a3x = new TH1D("a3x","", xBins,xMin,xMax);
  TH1D *a4x = new TH1D("a4x","", xBins,xMin,xMax);
  TH1D *a5x = new TH1D("a5x","", xBins,xMin,xMax);

  TH2D *a4v3j = new TH2D("a4v3j","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a4v5j = new TH2D("a4v5j","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a5v3j = new TH2D("a5v3j","", xBins,xMin,xMax,xBins,xMin,xMax);

  TH2D *a4v3x = new TH2D("a4v3x","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a4v5x = new TH2D("a4v5x","", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D *a5v3x = new TH2D("a5v3x","", xBins,xMin,xMax,xBins,xMin,xMax);

  TH2D* hJvX3 = new TH2D("hJvX3", "", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D* hJvX4 = new TH2D("hJvX4", "", xBins,xMin,xMax,xBins,xMin,xMax);
  TH2D* hJvX5 = new TH2D("hJvX5", "", xBins,xMin,xMax,xBins,xMin,xMax);

  //Loop over all events                                                                                                                             
  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
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

      std::vector<double> jmlAns;
      std::vector<double> xmnAns;
      
      // Begin Online
      hlt_wNeg->apply(inputCaloSample,inputPedestal,jmlAns);
      hlt_noNeg->apply(inputCaloSample,inputPedestal,xmnAns);

      if (jmlAns.size()>1) {

        a3j->Fill(jmlAns.at(0));
        a4j->Fill(jmlAns.at(1));
        a5j->Fill(jmlAns.at(2));
	
        a4v3j->Fill(jmlAns.at(1), jmlAns.at(0));
        a4v5j->Fill(jmlAns.at(1), jmlAns.at(2));
        a5v3j->Fill(jmlAns.at(2), jmlAns.at(0));
      }	

      if (xmnAns.size()>1) {

        a3x->Fill(xmnAns.at(0));
        a4x->Fill(xmnAns.at(1));
        a5x->Fill(xmnAns.at(2));
	
        a4v3x->Fill(xmnAns.at(1), xmnAns.at(0));
        a4v5x->Fill(xmnAns.at(1), xmnAns.at(2));
        a5v3x->Fill(xmnAns.at(2), xmnAns.at(0));
      }	

      if (jmlAns.size()>1 && xmnAns.size()>1) {
	hJvX3->Fill( jmlAns.at(0), xmnAns.at(0) );
	hJvX4->Fill( jmlAns.at(1), xmnAns.at(1) );
	hJvX5->Fill( jmlAns.at(2), xmnAns.at(2) );
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  gStyle->SetOptStat(0);

  hJvX3->GetXaxis()->SetTitle("No Neg. Corr. A3 [fC]");
  hJvX3->GetYaxis()->SetTitle("With Neg. Corr. A3 [fC]");
  hJvX3->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hJvX3.png");

  hJvX4->GetXaxis()->SetTitle("No Neg. Corr. A4 [fC]");
  hJvX4->GetYaxis()->SetTitle("With Neg. Corr. A4 [fC]");
  hJvX4->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hJvX4.png");

  hJvX5->GetXaxis()->SetTitle("No Neg. Corr. A5 [fC]");
  hJvX5->GetYaxis()->SetTitle("With Neg. Corr. A5 [fC]");
  hJvX5->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/hJvX5.png");

  a3j->GetXaxis()->SetTitle("No Neg. Corr. A3 [fC]");
  a3j->GetXaxis()->SetTitleSize(0.05);
  a3j->GetYaxis()->SetTitle("Counts");
  a3j->GetYaxis()->SetTitleSize(0.05);
  a3j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a3j.png");

  a4j->GetXaxis()->SetTitle("No Neg. Corr. A4 [fC]");
  a4j->GetXaxis()->SetTitleSize(0.05);
  a4j->GetYaxis()->SetTitle("Counts");
  a4j->GetYaxis()->SetTitleSize(0.05);
  a4j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4j.png");

  a5j->GetXaxis()->SetTitle("No Neg. Corr. A5 [fC]");
  a5j->GetXaxis()->SetTitleSize(0.05);
  a5j->GetYaxis()->SetTitle("Counts");
  a5j->GetYaxis()->SetTitleSize(0.05);
  a5j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5j.png");

  a4v3j->GetXaxis()->SetTitle("No Neg. Corr.  A4 [fC]");
  a4v3j->GetXaxis()->SetTitleSize(0.05);
  a4v3j->GetYaxis()->SetTitle("No Neg. Corr.  A3 [fC]");
  a4v3j->GetYaxis()->SetTitleSize(0.05);
  a4v3j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v3j.png");

  a4v5j->GetXaxis()->SetTitle("No Neg. Corr.  A4 [fC]");
  a4v5j->GetXaxis()->SetTitleSize(0.05);
  a4v5j->GetYaxis()->SetTitle("No Neg. Corr.  A5 [fC]");
  a4v5j->GetYaxis()->SetTitleSize(0.05);
  a4v5j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v5j.png");

  a5v3j->GetXaxis()->SetTitle("No Neg. Corr.  A5 [fC]");
  a5v3j->GetXaxis()->SetTitleSize(0.05);
  a5v3j->GetYaxis()->SetTitle("No Neg. Corr.  A3 [fC]");
  a5v3j->GetYaxis()->SetTitleSize(0.05);
  a5v3j->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5v3j.png");

  a3x->GetXaxis()->SetTitle("With Neg. Corr.  A3 [fC]");
  a3x->GetXaxis()->SetTitleSize(0.05);
  a3x->GetYaxis()->SetTitle("Counts");
  a3x->GetYaxis()->SetTitleSize(0.05);
  a3x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a3x.png");

  a4x->GetXaxis()->SetTitle("With Neg. Corr.  A4 [fC]");
  a4x->GetXaxis()->SetTitleSize(0.05);
  a4x->GetYaxis()->SetTitle("Counts");
  a4x->GetYaxis()->SetTitleSize(0.05);
  a4x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4x.png");

  a5x->GetXaxis()->SetTitle("With Neg. Corr.  A5 [fC]");
  a5x->GetXaxis()->SetTitleSize(0.05);
  a5x->GetYaxis()->SetTitle("Counts");
  a5x->GetYaxis()->SetTitleSize(0.05);
  a5x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5x.png");

  a4v3x->GetXaxis()->SetTitle("With Neg. Corr.  A4 [fC]");
  a4v3x->GetXaxis()->SetTitleSize(0.05);
  a4v3x->GetYaxis()->SetTitle("With Neg. Corr.  A3 [fC]");
  a4v3x->GetYaxis()->SetTitleSize(0.05);
  a4v3x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v3x.png");

  a4v5x->GetXaxis()->SetTitle("With Neg. Corr.  A4 [fC]");
  a4v5x->GetXaxis()->SetTitleSize(0.05);
  a4v5x->GetYaxis()->SetTitle("With Neg. Corr.  A5 [fC]");
  a4v5x->GetYaxis()->SetTitleSize(0.05);
  a4v5x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a4v5x.png");

  a5v3x->GetXaxis()->SetTitle("With Neg. Corr.  A5 [fC]");
  a5v3x->GetXaxis()->SetTitleSize(0.05);
  a5v3x->GetYaxis()->SetTitle("With Neg. Corr.  A3 [fC]");
  a5v3x->GetYaxis()->SetTitleSize(0.05);
  a5v3x->Draw();
  c1->SaveAs(TString(Plot_Dir.c_str())+"/a5v3x.png");  

  /*  int nMethod=5;
  
  std::vector<HLTv2*> vHltMethods_;
  std::vector<TH1D*> vHltA4Plots;
  std::vector<TH2D*> vHltvHlt;

  //I'm sorry :(
  std::vector<std::vector<double>> vHltAns(nMethod, std::vector<double>(10));

  int xBins=100, xMin=-5,xMax=15;

  for (int i=0; i<nMethod; i++) {
    vHltMethods_.push_back(new HLTv2);
  }
  char hname[50];

  sprintf(hname, "ts_testStand_slow");
  vHltMethods_[0]->Init(HcalTimeSlew::TestStand, HcalTimeSlew::Slow, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  vHltA4Plots.push_back(new TH1D(hname,"",xBins,xMin,xMax));
  vHltA4Plots[0]->GetXaxis()->SetTitle("A4 with Slow test stand TS param [fC]");

  sprintf(hname, "ts_testStand_med");
  vHltMethods_[1]->Init(HcalTimeSlew::TestStand, HcalTimeSlew::Medium, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  vHltA4Plots.push_back(new TH1D(hname,"",xBins,xMin,xMax));
  vHltA4Plots[1]->GetXaxis()->SetTitle("A4 with Med test stand TS param [fC]");

  sprintf(hname, "ts_testStand_fast");
  vHltMethods_[2]->Init(HcalTimeSlew::TestStand, HcalTimeSlew::Fast, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  vHltA4Plots.push_back(new TH1D(hname,"",xBins,xMin,xMax));
  vHltA4Plots[2]->GetXaxis()->SetTitle("A4 with Fast test stand TS param [fC]");

  sprintf(hname, "ts_data");
  vHltMethods_[3]->Init(HcalTimeSlew::Data, HcalTimeSlew::Medium, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  vHltA4Plots.push_back(new TH1D(hname,"",xBins,xMin,xMax));
  vHltA4Plots[3]->GetXaxis()->SetTitle("A4 with capped MC TS param [fC]");

  sprintf(hname, "ts_mc");
  vHltMethods_[4]->Init(HcalTimeSlew::MC, HcalTimeSlew::Medium, (HLTv2::NegStrategy)Neg_Charges, *pedSubFxn_);
  vHltA4Plots.push_back(new TH1D(hname,"",xBins,xMin,xMax));
  vHltA4Plots[4]->GetXaxis()->SetTitle("A4 with uncapped MC TS param [fC]");

  sprintf(hname, "mediumVslow");
  vHltvHlt.push_back(new TH2D(hname, "", xBins, xMin, xMax, xBins, xMin, xMax));
  vHltvHlt[0]->GetXaxis()->SetTitle("A4 with medium test stand param [fC]");
  vHltvHlt[0]->GetYaxis()->SetTitle("A4 with slow test stand param [fC]");

  sprintf(hname, "mediumVfast");
  vHltvHlt.push_back(new TH2D(hname, "", xBins, xMin, xMax, xBins, xMin, xMax));
  vHltvHlt[1]->GetXaxis()->SetTitle("A4 with medium test stand param [fC]");
  vHltvHlt[1]->GetYaxis()->SetTitle("A4 with fast test stand param [fC]");

  sprintf(hname, "testStandVdata");
  vHltvHlt.push_back(new TH2D(hname, "", xBins, xMin, xMax, xBins, xMin, xMax));
  vHltvHlt[2]->GetXaxis()->SetTitle("A4 with test stand TS param [fC]");
  vHltvHlt[2]->GetYaxis()->SetTitle("A4 with capped MC TS param [fC]");

  sprintf(hname, "testStandVmc");
  vHltvHlt.push_back(new TH2D(hname, "", xBins, xMin, xMax, xBins, xMin, xMax));
  vHltvHlt[3]->GetXaxis()->SetTitle("A4 with test stand TS param [fC]");
  vHltvHlt[3]->GetYaxis()->SetTitle("A4 with uncapped MC TS param [fC]");

  sprintf(hname, "mcVdata");
  vHltvHlt.push_back(new TH2D(hname, "", xBins, xMin, xMax, xBins, xMin, xMax));
  vHltvHlt[4]->GetXaxis()->SetTitle("A4 with uncapped MC TS param [fC]");
  vHltvHlt[4]->GetYaxis()->SetTitle("A4 with capped MC TS param [fC]");

  for (int jentry=0; jentry<Entries;jentry++) {
    fChain->GetEntry(jentry);
    for (int j = 0; j < (int)PulseCount; j++) {
      std::vector<double> inputCaloSample, inputPedestal;
      for (int i=0; i<10; i++) {
	inputCaloSample.push_back(Charge[j][i]+Pedestal[j][i]);
	inputPedestal.push_back(Pedestal[j][i]);
      }
      for (int k=0; k<nMethod; k++) {
	vHltMethods_[k]->apply(inputCaloSample,inputPedestal,vHltAns[k]);
	vHltA4Plots[k]->Fill(vHltAns[k].at(1));
      }
      vHltvHlt[0]->Fill(vHltAns[1].at(1),vHltAns[0].at(1));
      vHltvHlt[1]->Fill(vHltAns[1].at(1),vHltAns[2].at(1));
      vHltvHlt[2]->Fill(vHltAns[1].at(1),vHltAns[3].at(1));
      vHltvHlt[3]->Fill(vHltAns[1].at(1),vHltAns[4].at(1));
      vHltvHlt[4]->Fill(vHltAns[3].at(1),vHltAns[4].at(1));
    }
  }

  TCanvas *c1 = new TCanvas("c1");
  gStyle->SetOptStat(0);

  char fname[50];

  for (int i=0; i<nMethod; i++) {
    vHltA4Plots[i]->GetYaxis()->SetTitle("Counts");
    vHltA4Plots[i]->Draw("hist");
    sprintf(fname, "timeslew_plots/result_%i_%i.png", i, Condition);
    c1->SaveAs(fname);

    vHltvHlt[i]->Draw("colz");
    sprintf(fname, "timeslew_plots/comp_%i_%i.png", i, Condition);
    c1->SaveAs(fname);
  }
  */
}

void Analysis::Finish()
{

  fout->cd();
  fout->Write();
  fout->Close();
}
