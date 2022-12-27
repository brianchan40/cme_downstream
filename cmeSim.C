#include <iostream>

class toyParticle;
class toyEvent;
class doubleRatioAnalyzer;

TFile *outputFile;
using namespace std;
//-------------------
// parameter setting
//-------------------

const double primaryParticle_a1[5] = {0., 0.002, 0.0035, 0.005, 0.01};
// const double primaryParticle_a1[5]      = {0., 0.01,  0.02, 0.03,  0.04}; //old standard set
// const double primaryParticle_a1[5]    = {0., 0.2,   0.4,  0.6,   0.8}; //wider range
// const double primaryParticle_a1[8]    = {0., 0.02,   0.04,  0.1,   0.2, 0.4, 0.6, 0.8}; //widest range
// const double primaryParticle_v2[5]    = {0.3, 0.4,   0.6,  0.8,   1}; //wider range
const double primaryParticle_v2[5] = {0., 0.05, 0.1, 0.15, 0.2};
// const double primaryParticle_v2Slope[5] = {0., 0.025, 0.05, 0.075, 0.1};
const double primaryParticle_v2PtParA[5] = {0.075, 0.125, 0.175, 0.225, 0.275};
const double primaryParticle_v3v2Ratio[5] = {0., 1. / 5., 2. / 5., 3. / 5., 4. / 5.}; // v3 over v2 ratio
const double primaryParticle_v2_2ndGrp[5] = {0.01, 0.03, 0.05, 0.07, 0.09};
// const double primaryParticle_InOutTChangeInFrac[5]      = {0.0, 0.05,  0.1,  0.15,  0.2}; //too wide range
const double primaryParticle_InOutTChangeInFrac[5] = {0.0, 0.025, 0.05, 0.075, 0.1};

const double resonance_v2[5] = {0., 0.05, 0.1, 0.15, 0.2};
const double resonance_v3[5] = {0., 0.01, 0.05, 0.1, 0.2};
// const double       resonance_v3[5]      = {0., -0.01,  -0.05, -0.1 ,  -0.2}; //neg. values
// const double       resonance_v2Slope[5] = {0., 0.05,  0.1,  0.15,  0.2};
const double resonance_v2PtParA[5] = {0.075, 0.125, 0.175, 0.225, 0.275};
const double resonance_fixedPt[5] = {0.5, 1., 5., 10., 20.};
const double resonance_T_forTMtSpectra[5] = {0.8 * 0.317, 1. * 0.317, 1.2 * 0.317, 1.4 * 0.317, 1.6 * 0.317};
const double resonance_v3v2Ratio[5] = {0., -1. / 5., -2. / 5., -3. / 5., -4. / 5.}; // v2 over v3 ratio

const double resonance_temp[5] = {0.25, 0.317, 0.4, 1., 2.0};              // temperature
const double resonance_fixedPt_fqPaperFig3[5] = {0.5, 0.75, 1., 1.25, 2.}; // for case for fig3 of //https://arxiv.org/pdf/1803.02860.pdf
// const double  resonance_meanPt[5]       = {0.5,0.83,  1.5., 3.,    4.};
const double resonance_rho00[5] = {1. / 3. - 0.2,
                                   1. / 3. - 0.1,
                                   1. / 3.,
                                   1. / 3. + 0.1,
                                   1. / 3. + 0.2}; // rho_00 value for spin alignment.

const double resonance_rho00_closeUpView[5] = {1. / 3. - 0.2,
                                               1. / 3. - 0.15,
                                               1. / 3. - 0.1,
                                               1. / 3. - 0.05,
                                               1. / 3.}; // rho_00 value for spin alignment. Close-up view

const double etaAcceptance[5] = {-0.2, -0.4, -0.6, -0.8, -0.1};

//--------------------------
// reaction plane rotating
//--------------------------

double Chi(double res);
double res_keq1(double chi); // k==1
double res_keq2(double chi);
double generateDeltaPhi(double chi);
double PDF(double chi, double phi);

TRandom3 *rdm = new TRandom3();
rdm->SetSeed(0);

TH1D *hDeltaPhi = new TH1D("hDeltaPhi", "", 100, -4., 4.);

//------------
// main
//------------
void cmeSim(const int variableIdx, const TString jobId)
{

  /*
    gROOT->ProcessLine(".L toyParticle_cxx.so");
    gROOT->ProcessLine(".L toyEvent_cxx.so");
    gROOT->ProcessLine(".L eventPreparer_cxx.so");
    gROOT->ProcessLine(".L eventProcessor_cxx.so");
    gROOT->ProcessLine(".L BFEventProcessor_cxx.so");
    gROOT->ProcessLine(".L RCorrEventProcessor_cxx.so");
    gROOT->ProcessLine(".L doubleRatioAnalyzer_cxx.so");
    gROOT->ProcessLine(".L qaMonitor_cxx.so");
  */
  gSystem->Load("libtoy");
  //      double varyingValue =  resonance_v2[variableIdx];
  //      double varyingValue =  resonance_v2Slope[variableIdx];
  //      double varyingValue =  resonance_v3[variableIdx];
  //      double varyingValue =  resonance_fixedPt[variableIdx];
  //	    double varyingValue =  primaryParticle_v2PtParA[variableIdx];
  //	    double varyingValue =  resonance_v2PtParA[variableIdx];
  //	    double varyingValue =  primaryParticle_v3v2Ratio[variableIdx];
  //	    double varyingValue =  resonance_v3v2Ratio[variableIdx];

  //	    double varyingValue =  resonance_temp[variableIdx];
  //      double varyingValue =  resonance_fixedPt_fqPaperFig3[variableIdx];
  //      double varyingValue =  resonance_T_forTMtSpectra[variableIdx];
  //      double varyingValue =  resonance_meanPt[variableIdx];
  double varyingValue = resonance_rho00[variableIdx];
  //      double varyingValue =  resonance_rho00_closeUpView[variableIdx];
  //      double varyingValue =  primaryParticle_v2[variableIdx];
  //      double varyingValue =  primaryParticle_v2Slope[variableIdx];
  //      double varyingValue =  primaryParticle_a1[variableIdx];
  //      double varyingValue =  etaAcceptance[variableIdx];
  //      double varyingValue =  primaryParticle_v2_2ndGrp[variableIdx];
  //      double varyingValue =  primaryParticle_InOutTChangeInFrac[variableIdx];

  // default fixed value for particle V2 0.05, for resonance v2 0.15. for a1 0.03.

  eventPreparer *evtPr = new eventPreparer();
  int NEvent = 20000;                 // 10500 evts per job, 190 jobs per set. Total 10500x190 ~ 2M evts per set
  //int NEvent = 100000;
  int halfMult = 195;                 // 195; //refer to 30-40% in table II of https://journals.aps.org/prc/pdf/10.1103/PhysRevC.79.034909
  int resMult = int(halfMult * 0.17); // ratio of rho to NPositive or NNegative is 17% PhysRevLett.92.092301. for 40-80% centrality.
  // int resMult = 0;

  evtPr->setNPositive(halfMult - resMult); //"primary" positives
  evtPr->setNNegative(halfMult - resMult);
  evtPr->setNResonance(resMult);

  //  evtPr->setNPositive_2ndGrp(97); //a feature that is rarely used.
  //  evtPr->setNNegative_2ndGrp(97);

  //========== v2 for primary
  evtPr->setV2Primary(0.05); // 0.05 for fixed value //primary (non-resonance decayed) particle v2
  //  evtPr->setV2PtParametersPrimary(0.125,0.17,0.177,0.031); //default (0.125,0.17,0.177,0.031);
  // If called, override the constant v2 above and v3(pt) is set to be 1/5 of v2(pt)
  // default v2(pt) para. from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.95.051901   30-40% cent.
  // v3 can be set as ~1/5 of v2 values. This based on Qiye's SQM slides below.
  // https://indico.cern.ch/event/204432/contributions/1490921/attachments/309863/432561/SQM2013_qyshou_presentation_v2.pdf
  // A better reference for v3/v2 ~ 1/5 : Nucl.Phys. A931 (2014) 758-762   Qiyes' QM 2014 proceeeding.
  // evtPr->setV2PtParAPrimary(varyingValue);

  // v2 for second group of primaries. Rarelly used!
  // evtPr->setV2Primary_2ndGrp(varyingValue);

  //========== v2 for resonance

  evtPr->setV2Resonance(0.05); // 0.15 for fixed value (but 0.06 was used in Fuqiang's paper PhysRevC.98.034904)
  // evtPr->setV2PtParametersResonance(0.125,0.17,0.177,0.031);
  // If called, override the constant v2 above and v3(pt) (below) is set to be 1/5 of v2(pt)
  // evtPr->setV2PtParAResonance(varyingValue);

  //========== v3 for primary,.
  // evtPr->setV3Primary(0.0);   // if called, use value, otherwise set to be 1/5 of v2 by default or set by below
  // evtPr->setV3V2RatioPrimary(varyingValue); // default 1/5.

  //========== v3 for  resonance.
  // evtPr->setV3Resonance(0.0); // if called, use value, otherwise set to be 1/5 of v2 by default or set by below.
  // evtPr->setV3V2RatioResonance(varyingValue); //default 1/5.

  //========== rho00 (for Res.) and a1 (for primary)
  evtPr->setRho00Resonance(varyingValue); // 1./3. default value (no alignment)
                                          // evtPr->setRho00Resonance(1./3.-0.1); // 1./3. default value (no alignment)
  // evtPr->setA1(0.015);  //0.015 for fixed value a1

  //======== spectra for primary, resonance.
  // type 1: Bose-Einstein 2: simple exponential 3: pt exp, 4: pt gaussian 5: pt3 exp 6: mt exp 7: TMt exp 8: Boltzmann
  evtPr->setSpectraTypePion(1);
  evtPr->setSpectraTypeRho(7);

  //======== In and out T change in fraction for primordial pion. Can be used only together with fixed v2 value, not v2(pt). not often used.
  // evtPr->setInOutTchangeInFraction(varyingValue); //default (not called) is 0.

  // these two lines below are for setting temperatures, not oftenly used. Default T is good enough.
  // evtPr->setTPion(0.2125); // 0.195 default value. Temperature of pions. see Table I PhysRevC.95.051901 value for 30-40% cent.
  // evtPr->setTRho(varyingValue); // 0.317 default vaule. Temp.of rho., No effect if setFixedPtResonance() is called with a nonzero value.

  //======== stuff that are not oftenly used.
  // evtPr->setFixedPtResonance(varyingValue); //if called will override spectra setting.
  // evtPr->setEtaAcceptance(varyingValue); // a negative values means not cut on eta. default -999 if not called.
  bool doEP = false;       // add by dshen 07/18
  bool doRotation = false; // add by dshen 08/08;
  double resolution = 0.6;

  doubleRatioAnalyzer *BFAnalyzer = new doubleRatioAnalyzer(1, evtPr->eventMult(), doEP, doRotation);
  doubleRatioAnalyzer *RCorrAnalyzer = new doubleRatioAnalyzer(2, evtPr->eventMult(), doEP, doRotation);
  //   doubleRatioAnalyzer* RestFrameBFAnalyzer = new doubleRatioAnalyzer(3, evtPr->eventMult(), doEP, doRotation);
  deltaGamma *DGAnalyzer = new deltaGamma();

  qaMonitor *monitor = new qaMonitor(); // contains varying variables for drawGraph, must be included in this run.
  monitor->setVaryingValueForPlotting(varyingValue);

  for (int ievt = 0; ievt < NEvent; ievt++)
  {
    if (ievt % 100 == 0)
      cout << "ievt " << ievt << endl;
    toyEvent *evt = evtPr->prepareOneEvent();

    double chi = Chi(resolution);
    double deltaPhi = generateDeltaPhi(chi) / 2.;
    hDeltaPhi->Fill(deltaPhi);

    monitor->check(evt); // pls makesure cuts for EP used here same as BFAnalyzer
    BFAnalyzer->analyze(evt, deltaPhi);
    //	RestFrameBFAnalyzer->analyze(evt,deltaPhi);
    RCorrAnalyzer->analyze(evt, deltaPhi);
    DGAnalyzer->analyze(evt);

    if (evt)
      delete evt;
  }
  TString fn("cmeTest");
  cout << fn << endl;
  fn += variableIdx;
  cout << fn << endl;
  fn.Append(TString::Format("_job%s", jobId.Data()).Data());
  cout << fn << endl;
  // fn += jobId;
  fn.Append(".root");
  cout << fn << endl;
  outputFile = new TFile(fn.Data(), "RECREATE");
  BFAnalyzer->summarizeAndWrite(outputFile);
  //     RestFrameBFAnalyzer->summarizeAndWrite(outputFile);
  RCorrAnalyzer->summarizeAndWrite(outputFile);
  monitor->summarizeAndWrite(outputFile); // common QA histo
  DGAnalyzer->save(outputFile);
  outputFile->cd();
  hDeltaPhi->Write();
  outputFile->Close();

  cout << "I have saved file" << endl;
}

//------------------------------
double Chi(double res)
{
  double chi = 2.0;
  double delta = 1.0;

  for (int i = 0; i < 15; i++)
  {
    while (res_keq1(chi) < res)
    {
      chi += delta;
    }
    delta = delta / 2.;
    while (res_keq1(chi) > res)
    {
      chi -= delta;
    }
    delta = delta / 2.;
  }
  // cout<<"chi="<<chi<<endl;
  return chi;
}

double res_keq1(double chi)
{
  double A = sqrt(TMath::Pi() / 2.) / 2.;
  double dog = pow(chi, 2) / 4.;

  double res = A * chi * exp(-dog) * (TMath::BesselI0(dog) + TMath::BesselI1(dog));
  return res;
}

double res_keq2(double chi)
{
  double A = sqrt(TMath::Pi() / 2.) / 2.;
  double dog = pow(chi, 2) / 4.;
  double pi_half = TMath::Pi() / 2.;

  double bessel_1to2 = sqrt(dog / pi_half) * TMath::SinH(dog) / dog;
  double bessel_3to2 = sqrt(dog / pi_half) * (TMath::CosH(dog) / dog - TMath::SinH(dog) / pow(dog, 2));
  double res = A * chi * exp(-dog) * (bessel_1to2 + bessel_3to2);
  return res;
}

double generateDeltaPhi(double chi)
{
  double deltaPhi = 0;
  double dog = 0;
  do
  {
    deltaPhi = ((rdm->Rndm()) * 2. - 1.) * TMath::Pi();                                        // uniform distribution of (-pi,pi)
    dog = (rdm->Rndm()) * (1. / (2. * TMath::Pi())) * (1 + sqrt(TMath::Pi() / 2.) * 2. * chi); // the benifit of do while loop is that we do not need a exact maximum value of PDF, shortage is we take the risk of dead loop
  } while (dog > PDF(chi, deltaPhi));
  return deltaPhi;
}

double PDF(double chi, double phi)
{
  double PI = TMath::Pi();
  double value = (1. / (2. * PI)) * (exp(-pow(chi, 2) / 2.) + sqrt(PI / 2.) * chi * cos(phi) * exp(-pow(chi * sin(phi), 2) / 2.) * (1 + TMath::Erf(chi * cos(phi) / sqrt(2))));
  return value;
}
