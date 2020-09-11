////postAnalyzer.C
//For use with Ntuples made from ggNtuplizer
//Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
//
//To compile using rootcom to an executable named 'analyze':
//$ ./rootcom postAnalyzer analyze
//
//To run, assuming this is compiled to an executable named 'analyze':
//$ ./analyze /hdfs/store/user/jjbuch/LatestNtuples/ /afs/hep.wisc.edu/user/jjbuchanan/private/CMSSW_7_4_9/src/output.root -1 10000
//Runs over every event in the folder LatestNtuples, reporting progress every 10000 events
//and storing the resulting histograms in the file output.root.
//
//To plot, for example, single photon trigger efficiency as a function of photon pt:
//$ root -l
//root[0] TFile *f = new TFile("output.root");
//root[1] TGraphAsymmErrors *efficiency = new TGraphAsymmErrors((TH1F*)f->Get("Photon_Et_300_2"),(TH1F*)f->Get("Photon_Et_300_1"));
//root[2] efficiency->Draw("AP")
//root[3] efficiency->SetTitle("Single photon trigger efficiency")
//root[4] efficiency->GetXaxis()->SetTitle("Photon p_{T}")
//root[5] efficiency->Draw("AP")
//

#define postAnalyzer_cxx
#include "postAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
	Long64_t maxEvents = atof(argv[3]);
	if (maxEvents < -1LL)
	{
		std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
		return 1;
	}
	int reportEvery = atof(argv[4]);
	if (reportEvery < 1)
	{
		std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
		return 1;
	}
	postAnalyzer t(argv[1],argv[2]);
	t.Loop(maxEvents,reportEvery);
	return 0;
}

void postAnalyzer::Loop(Long64_t maxEvents, int reportEvery){
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  int nBackupTriggerEvents, nBTMediumEvents, nBTMediumHLTsinglePhoEvents, nEffPhoptden, nEffPhoptnum, nEffMETden, nEffMETnum;
  nBackupTriggerEvents = nBTMediumEvents = nBTMediumHLTsinglePhoEvents = nEffPhoptden = nEffPhoptnum = nEffMETden = nEffMETnum = 0;
  int nHLTPassed, nGoodPhotonPassed, nPhotonPtPassed, nMETPassed, nDPhiPassed, nqcdden,nqcdnum;
  nHLTPassed = nGoodPhotonPassed = nPhotonPtPassed = nMETPassed = nDPhiPassed = nqcdden= nqcdnum= 0;
  
  std::vector<int> phoCand;
  phoCand.clear();
  
  std::vector<int> phoCand1;
  phoCand1.clear();
  
  
  std::vector<int> qcdden;
  qcdden.clear();
  //mark
  std::vector<int> phoCand_test;
  phoCand_test.clear();

  bool debug=true;
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  //Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;
  
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++){
    
    event_.clear();
    event_info.clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //=1.0 for real data
    double event_weight=1.0;
  //mark
    phoCand_test = getPhoCand1(1.4442,1);  
//if(phoCand_test.size() > 0)fillHistos(3,event_weight,phoCand_test[0]);
    
    // Trigger Efficiencies
    if((HLTPho>>4&1 == 1)||(HLTPho>>5&1 == 1) || (HLTPho>>6&1 == 1)) //HLT_Photon75_v*||HLT_Photon90_v*||HLT_Photon120_v*
      {
	nBackupTriggerEvents++;
	phoCand   = getPhoCand1(1.4442,1);
//	if(phoCand.size() >0 )fillHistos(4,event_weight,phoCand[0]);//none
	if(phoCand.size() >0 && pfMET>40)  //140 -> 190 -> 160
	  {
	    nBTMediumEvents++;
	    fillHistos(1,event_weight,phoCand[0]);
	    if((HLTJet>>11&1 == 1)||(HLTJet>>12&1 == 1)|| (HLTJet>>13&1 == 1) || (HLTPho>>11&1 == 1))
	      {
		nBTMediumHLTsinglePhoEvents++;
		fillHistos(2,event_weight,phoCand[0]);
	      }
	    
	    if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
	      {
		nEffMETnum++;
		fillHistos(3,event_weight,phoCand[0]);
	      }
	    
	  }
	if(phoCand.size() >0 && phoCand[0]>175.)
	  {
	    fillHistos(4,event_weight,phoCand[0]);
	    if((HLTPho>>22&1 == 1))
	      {
		fillHistos(5,event_weight,phoCand[0]);
	      }	
	    
	  }
	
	
      }
    
    
    //Trigger Efficiency for low pt triggers
    if((HLTPho>>0&1 == 1)||(HLTPho>>1&1 == 1)) //HLT_Photon75_v*||HLT_Photon90_v*||HLT_Photon120_v*
      {
	phoCand   = getPhoCand1(1.4442,1);
	if(phoCand.size() >0 && pfMET>100.0)
	  {
	    fillHistos(6,event_weight,phoCand[0]);
	    if((HLTPho>>23&1 == 1)||(HLTPho>>24&1 == 1)|| (HLTPho>>25&1 == 1) || (HLTPho>>26&1 == 1) || (HLTPho>>28&1 == 1))
	      {
		fillHistos(7,event_weight,phoCand[0]);
	      }
	    
	  }
	
	if(phoCand.size() >0 && phoCand[0]>140)
	  {
	    fillHistos(8,event_weight,phoCand[0]);
	    if((HLTPho>>23&1 == 1)||(HLTPho>>24&1 == 1)|| (HLTPho>>25&1 == 1) || (HLTPho>>26&1 == 1) || (HLTPho>>28&1 == 1))
	      {
		fillHistos(9,event_weight,phoCand[0]);
	      }
	    
	  }
	
	if(phoCand.size() >0)
	  {
	    fillHistos(10,event_weight,phoCand[0]);
	    if((HLTPho>>27&1 == 1)||(HLTPho>>29&1 == 1))
	      {
		fillHistos(11,event_weight,phoCand[0]);
	      }	
	  }  
	
      }
    
    
    //Sequential cuts
    phoCand1   = getPhoCand(175,1.4442,1);
    if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
      {
	nHLTPassed++;
	if(phoCand1.size() >0)
	  {
	    nGoodPhotonPassed++;
	    fillHistos(12,event_weight,phoCand1[0]);
	    if(pfMET>140)
	      {
		nMETPassed++;
		fillHistos(13,event_weight,phoCand1[0]);
		if(DeltaPhi(phoPhi->at(phoCand1[0]),pfMETPhi)>2.0)
		  {
		    nDPhiPassed++;
		    fillHistos(14,event_weight,phoCand1[0]);
		    std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
		  }
	      }
	  }
      }
    
    //qcd sequential cuts
    qcdden   = getQcdden(175,1.4442,1);
    
    if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
      {
	if(phoCand1.size() >0 && pfMET<30)
	  {
	    nqcdnum++;
	    std::cout<<"qcd num : run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
	    fillHistos(15,event_weight,phoCand1[0]);
	    if(qcdden.size()>0 && pfMET<30)
	      {
		nqcdden++;
		fillHistos(16,event_weight,qcdden[0]);
		std::cout<<"qcd den : run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
	      }
	  }
      }
    
    tree->Fill();
    
    if (jentry%reportEvery == 0)
      {
	std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
  }
  
  if((nentriesToCheck-1)%reportEvery != 0)std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  std::cout<<"All events checked."<<std::endl;
  //Report
  std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;
  std::cout << std::endl;
  std::cout<<"Total number of events: "<<nTotal<<std::endl;
  std::cout << "Number of events inspected: " << nTotal << std::endl;
  std::cout<<std::endl;
  cout<<"nBackupTriggerEvents: "<<nBackupTriggerEvents<<endl;
  cout<<"nBTMediumEvents: "<<nBTMediumEvents<<endl;
  cout<<"nBTMediumHLTsinglePhoEvents: "<<nBTMediumHLTsinglePhoEvents<<endl;
  cout<<"nEffMETnum: "<<nEffMETnum<<endl;
  cout<<endl;
  cout<<"nHLTPassed: "<<nHLTPassed<<endl;
  cout<<"nGoodPhotonPassed: "<<nGoodPhotonPassed<<endl;
  cout<<"nMETPassed: "<<nMETPassed<<endl;
  cout<<"nDPhiPassed: "<<nDPhiPassed<<endl;
  cout<<"nqcdnum: "<<nqcdnum<<endl;
  cout<<"nqcdden: "<<nqcdden<<endl;
  cout<<" This is  the end of code "<< endl;
}

void postAnalyzer::BookHistos(const char* file2){
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();
  Float_t PtBins[34]={75.0, 100.0, 125.0, 145.0, 155.0, 165.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 425.0, 450.0, 475.0, 500.0, 550.0, 600.0, 700.0, 800.0, 1000.0, 9999.0};
//  Float_t PtBins[12]={75.0,100.0,125.0,145., 155.,165., 175.,190.,250., 400., 700.0,1000.0};
  Float_t MetBins[9]={130.0, 150.0, 170.0, 190.0, 250.0, 400.0, 500.0, 700.0, 1000.0};
//  Float_t MetBins[8]={130., 150., 170., 190., 250., 400., 700.0,1000.0};
  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<25; i++)
    {
      char ptbins[100];
      sprintf(ptbins, "_%d", i);
      std::string histname(ptbins);
      //h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
      h_photon_Et[i] = new TH1F(("Photon_Et"+histname).c_str(), "Photon_Et",100,0,1000);h_photon_Et[i]->Sumw2();
      h_photon_Et_range[i] = new TH1F(("Photon_Et_range"+histname).c_str(), "Photon_Et",33,PtBins);h_photon_Et_range[i]->Sumw2();
      h_photon_eta[i] = new TH1F(("Photon_eta"+histname).c_str(), "Photon_eta",40,-1.4442,1.4442);h_photon_eta[i]->Sumw2();
      h_photon_SCEta[i] = new TH1F(("Photon_SCeta"+histname).c_str(), "Photon_SCeta",40,-1.4442,1.4442);h_photon_SCEta[i]->Sumw2();
      //h_photon_phi[i] = new TH1F(("Photon_phi"+histname).c_str(), "Photon_phi", 64,0,3.2);h_photon_phi[i]->Sumw2();
      //h_photon_SCPhi[i] = new TH1F(("Photon_SCphi"+histname).c_str(), "Photon_SCphi", 64,0,3.2);h_photon_SCPhi[i]->Sumw2();
      //h_photon_IDbit[i] = new TH1F(("Photon_ID_bit"+histname).c_str(), "Photon_ID_bit",8,0,8);h_photon_IDbit[i]->Sumw2();
      h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",8,MetBins);h_pfMET[i]->Sumw2();
      h_pfMET_300[i] = new TH1F(("h_pfMET_300"+histname).c_str(), "pfMET",25,0,300);h_pfMET_300[i]->Sumw2();
      h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",40,0,3.2);h_dPhi[i]->Sumw2();
      h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",120,0,120);h_nJet[i]->Sumw2();
      h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();
      h_leadingJetPt_300[i] = new TH1F(("leadingJetPt_300"+histname).c_str(),"leadingJetPt_300",25,0,300);h_leadingJetPt_300[i]->Sumw2();
      h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();
    }
}

//Fill the sequential histos at a particular spot in the sequence
void postAnalyzer::fillHistos(int histoNumber, double event_weight,int index){
  
  h_photon_Et[histoNumber]->Fill(phoEt->at(index),event_weight);
  h_photon_Et_range[histoNumber]->Fill(phoEt->at(index),event_weight);
  h_photon_eta[histoNumber]->Fill(phoEta->at(index),event_weight);
  h_photon_SCEta[histoNumber]->Fill(phoSCEta->at(index),event_weight);
  //h_photon_phi[histoNumber]->Fill(phoPhi->at(0),event_weight);
  //h_photon_SCPhi[histoNumber]->Fill(phoSCPhi->at(0),event_weight);
  h_pfMET[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET_300[histoNumber]->Fill(pfMET,event_weight);
  double dPhi = DeltaPhi(phoPhi->at(index),pfMETPhi);
  h_dPhi[histoNumber]->Fill(dPhi,event_weight);
  h_nJet[histoNumber]->Fill(nJet,event_weight);
  //if(nJet>0){
  //	h_leadingJetPt[histoNumber]->Fill(jetPt->at(0),event_weight);
  //	h_leadingJetPt_300[histoNumber]->Fill(jetPt->at(0),event_weight);
  //	h_leadingJetEta[histoNumber]->Fill(jetEta->at(0),event_weight);
  //}
  
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzer::DeltaPhi(double phi1, double phi2){
  double pi = 3.14159265359;
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}


//---------------------------------------------------                                                                                                                                
// get a photon candiate based on pt eta and isolation                                                                                                                               
//----------------------------------------------------                                                                                                                               

std::vector<int> postAnalyzer::getPhoCand(double phoPtCut, double phoEtaCut, Short_t isoBit){
  std::vector<int> tmpCand;
  tmpCand.clear();
  
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++){
    bool kinematic = (*phoEt)[p] > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
    //Short_t IsoPass;
    //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
    bool photonId = (
		     ((*phoHoverE)[p]                <  0.0232   ) &&
		     ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.01040 ) &&
		     ((*phohasPixelSeed)[p]              ==  0      ) &&
         	     ( TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 0.584 ) &&        
         ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (0.321 + (0.0112 * (*phoEt)[p]) + (0.000028 * pow((*phoEt)[p], 2.0))) )  &&
         ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (2.141 + (0.0043 * (*phoEt)[p])) ) 
   ); 
    if(photonId && kinematic){
      tmpCand.push_back(p);
    }                                                                                                                                                              
  }                                                                                                                                                                
    return tmpCand;
}

std::vector<int> postAnalyzer::getPhoCand1(double phoEtaCut, Short_t isoBit){
  
  std::vector<int> tmpCand;
  tmpCand.clear();
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++) {
    bool kinematic = fabs((*phoSCEta)[p])<phoEtaCut;
    //Short_t IsoPass;
    //if(isoBit<=2)IsoPass= ((*phoIDbit)[p])>>isoBit&1;
    bool photonId = (
		     ((*phoHoverE)[p]                <  0.0232   ) &&
		     ((*phoSigmaIEtaIEtaFull5x5)[p]  <  0.01040 ) &&
		     ((*phohasPixelSeed)[p]              ==  0      ) &&
         	     ( TMath::Max( ( (*phoPFChIso)[p] - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 0.584 ) &&        
         ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (0.321 + (0.0112 * (*phoEt)[p]) + (0.000028 * pow((*phoEt)[p], 2.0))) )  &&
         ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (2.141 + (0.0043 * (*phoEt)[p])) ) 
   ); 
    
    if(photonId && kinematic){
      tmpCand.push_back(p);
    }                                                                                                                                                              
  }                                                                                                                                                               
  return tmpCand;
}


std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){
  std::vector<int> tmpCand;
  tmpCand.clear();
  
  //Loop over photons                                                                                                                                                             
  for(int p=0;p<nPho;p++){

//    Float_t phoEt = ((*phoSCRawE)[p]/TMath::CosH((*phoSCEta)[p]));
    //Fail loose iso
    bool passChIsoLoose = TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < 3.32;
    bool passNeuIsoLoose = TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (10.910 + (0.0148* (*phoEt)[p]) + (0.000028 * pow((*phoEt)[p], 2.0)));
    bool passPhoIsoLoose = TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (3.630 + (0.0053 * (*phoEt)[p]));

    if(!passChIsoLoose || !passNeuIsoLoose || !passPhoIsoLoose)
    {
      bool kinematic = (*phoEt)[p] > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
      //"Very loose" ID cuts with inverted shape cut
      bool photonID = (
        ((*phoSigmaIEtaIEtaFull5x5)[p]  >  0.01040 ) &&
        ((*phoHoverE)[p]                <  0.05   ) &&
        ((*phohasPixelSeed)[p]              ==  0      ) &&
        ( TMath::Max( ( (*phoPFChIso)[p]  - rho*EAcharged((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(3.32) , 0.20*(*phoEt)[p]) )  &&
        ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(10.910 + (0.0148 * (*phoEt)[p]) + (0.000028 * pow((*phoEt)[p], 2.0))) , 0.20*(*phoEt)[p]) )  &&
        ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < TMath::Min(5.0*(3.630+ (0.0053 * (*phoEt)[p])) , 0.20*(*phoEt)[p]) )
      );

    if(photonID && kinematic){
      tmpCand.push_back(p);
    }                        }                                                                                                                                      
  }                                                                                                                                                                
  return tmpCand;
}

// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5
Double_t postAnalyzer::EAcharged(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0456;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0500;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0340;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0383;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0339;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0303;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0240;
  
  return EffectiveArea;
}

Double_t postAnalyzer::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.0599;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.0819;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0696;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.0360;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.0462;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.0656;

  return EffectiveArea;
}

Double_t postAnalyzer::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.;
  if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.1271;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.1101;
  if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.0756;
  if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.1175;
  if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.1498;
  if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.1857;
  if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.2183;

  return EffectiveArea;
}

