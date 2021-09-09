#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TVector3.h"
#include "TRandom.h"
#include <Math/UnaryOperators.h>
#include <iostream>
#include <vector>
#include <map>
#include <bitset>
#include "TLorentzVector.h"
#include <numeric>
#include <algorithm> 
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<int> >+;
#endif
#ifdef __MAKECINT__
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class vector< vector < TLorentzVector > >+;
#endif

using namespace std;

//jet_comb("/squark1/kumari/new_sample/user.aknue.346344.PhPy8EG.DAOD_TOPQ1.e7148_s3126_r9364_p4346._TTH_PFlow_212171_V1_out_root/user.aknue.25595592._00000*.out.root","out_jetComb.root")

enum JetTruthMatchDefinition{JetMatchLeadingBHiggs=0,
			     JetMatchSubLeadingBHiggs=1,
			     JetMatchBLepTop=2,
			     JetMatchBHadTop=3,
			     JetMatchLeadingJetW=4,
			     JetMatchSubLeadingJetW=5};



void jet_comb(TString input_file,TString outfile){

  TChain* nominal_Loose = new TChain("nominal_Loose");//emtopo
  nominal_Loose->Add(input_file);
 
  TFile* f_new = TFile::Open(outfile);
  f_new = TFile::Open(outfile,"RECREATE");

  TH1F* hist_truthPt = new TH1F("h_truthPt","h_truthPt",100,0,1000000 );  
  TH1F* hist_truthPt_inv = new TH1F("h_truthPt_inv","h_truthPt_inv",100,0,1000000 );  

  int _nJets;
  int _nElectrons;
  int _nMuons;
  int _eventNumber;
  int _nBTags_DL1r_85;
  int _nBTags_DL1r_70;
  float _weight_mc;
  
  //study variable
  float _truth_higgs_pt;
  float _TTHReco_withH_best_Higgs_pt;
  float _TTHReco_best_Higgs_pt;

  //jet//
  vector<float> *_jet_pt;
  vector<float> *_jet_eta;
  vector<float> *_jet_phi;
  vector<float> *_jet_e;
  vector<int> *_jet_truthmatch;
  vector<int> *_jet_tagWeightBin;
  vector<char> *_jet_isbtagged_DL1r_85;
  vector<char> *_jet_isbtagged_DL1r_70;

  //lepton//

  vector<float> *_el_pt;
  vector<float> *_el_eta;
  vector<float> *_el_phi;
  vector<float> *_el_e;
  vector<float> *_el_charge;
  vector<float> *_mu_pt;
  vector<float> *_mu_eta;
  vector<float> *_mu_phi;
  vector<float> *_mu_e;
  vector<float> *_mu_charge;

  float _met_met;
  float _met_phi;

  nominal_Loose->SetBranchAddress("nJets",&_nJets);
  nominal_Loose->SetBranchAddress("nElectrons",&_nElectrons);
  nominal_Loose->SetBranchAddress("nMuons",&_nMuons);
  nominal_Loose->SetBranchAddress("nBTags_DL1r_85",&_nBTags_DL1r_85);
  nominal_Loose->SetBranchAddress("nBTags_DL1r_70",&_nBTags_DL1r_70);
  nominal_Loose->SetBranchAddress("eventNumber",&_eventNumber);
  nominal_Loose->SetBranchAddress("weight_mc",&_weight_mc);
  nominal_Loose->SetBranchAddress("jet_pt",&_jet_pt);
  nominal_Loose->SetBranchAddress("jet_eta",&_jet_eta);
  nominal_Loose->SetBranchAddress("jet_phi",&_jet_phi);
  nominal_Loose->SetBranchAddress("jet_e",&_jet_e);
  nominal_Loose->SetBranchAddress("jet_truthmatch",&_jet_truthmatch);
  nominal_Loose->SetBranchAddress("jet_tagWeightBin",&_jet_tagWeightBin);
  nominal_Loose->SetBranchAddress("jet_isbtagged_DL1r_85",&_jet_isbtagged_DL1r_85);
  nominal_Loose->SetBranchAddress("jet_isbtagged_DL1r_70",&_jet_isbtagged_DL1r_70);
  nominal_Loose->SetBranchAddress("el_pt",&_el_pt);
  nominal_Loose->SetBranchAddress("el_eta",&_el_eta);
  nominal_Loose->SetBranchAddress("el_phi",&_el_phi);
  nominal_Loose->SetBranchAddress("el_e",&_el_e);
  nominal_Loose->SetBranchAddress("el_charge",&_el_charge);
  nominal_Loose->SetBranchAddress("mu_pt",&_mu_pt);
  nominal_Loose->SetBranchAddress("mu_eta",&_mu_eta);
  nominal_Loose->SetBranchAddress("mu_phi",&_mu_phi);
  nominal_Loose->SetBranchAddress("mu_e",&_mu_e);
  nominal_Loose->SetBranchAddress("mu_charge",&_mu_charge);
  nominal_Loose->SetBranchAddress("met_met",&_met_met);
  nominal_Loose->SetBranchAddress("met_phi",&_met_phi);
  
  nominal_Loose->SetBranchAddress("truth_higgs_pt",&_truth_higgs_pt);
  nominal_Loose->SetBranchAddress("TTHReco_withH_best_Higgs_pt",&_TTHReco_withH_best_Higgs_pt);
  nominal_Loose->SetBranchAddress("TTHReco_best_Higgs_pt",&_TTHReco_best_Higgs_pt);

  int nentries = nominal_Loose->GetEntries();

  TTree* tree_new=new TTree("nominal_Loose_new","nominal_Loose_new");
  tree_new->SetDirectory(0);

  vector<float> _comb_higgs_mass;
  vector<float> _comb_higgs_pt;
  vector<float> _comb_hadW_mass;
  vector<float> _comb_hadtop_mass;
  vector<float> _comb_leptop_mass;
  vector<float> _comb_hadWblepTop_mass;
  vector<float> _comb_minbhadTopqhadW_dR;
  vector<float> _comb_hadWblepTop_dR;
  vector<float> _comb_blepTopbhadTop_dR;
  vector<float> _comb_bhadTopq2hadW_dR;
  vector<float> _comb_hadWbhadTop_dR;
  vector<float> _comb_diff_mindRbhadTopqhadW_dRlepblepTop;
  vector<float> _comb_Higgsq1hadW_mass;
  vector<float> _comb_bbHiggs_dR;
  vector<float> _comb_bhadTopq1hadW_dR;
  vector<float> _comb_qqhadW_dR;
  vector<float> _comb_lepbhadTop_dR;
  vector<float> _comb_lepb1Higgs_dR;
  vector<float> _comb_lepWbhadTop_mass;
  vector<float> _comb_lepblepTop_dR;

  vector<int>   _comb_bh1_index;
  vector<int>   _comb_bh2_index;
  vector<int>   _comb_bhadt_index;
  vector<int>   _comb_blept_index;
  vector<int>   _comb_qW1_index;
  vector<int>   _comb_qW2_index;
  vector<int>   _comb_truthMatchPattern;

  float _higgs_truthPt;
  float _best_higgs_Pt;
  float _best_higgs_Pt_withH;

  int _nBTags_DL1r85;
  int _nBTags_DL1r70;

  float _weight_mc_new;
  vector<float> _comb_weight_mc;
  vector<float> _comb_truthpt_rewg;

  float  _truthpt_weights;
  vector<float> _comb_higgs_truthPt;
  vector<float> _comb_best_higgs_Pt;
  vector<float> _comb_best_higgs_Pt_withH;

  tree_new->Branch("truthpt_weights",&_truthpt_weights);
  tree_new->Branch("weight_mc_new",&_weight_mc_new);
  tree_new->Branch("comb_weight_mc", &_comb_weight_mc);
  tree_new->Branch("comb_truthpt_rewg", &_comb_truthpt_rewg);
  tree_new->Branch("comb_higgs_mass", &_comb_higgs_mass);
  tree_new->Branch("comb_higgs_pt", &_comb_higgs_pt);
  tree_new->Branch("comb_hadW_mass", &_comb_hadW_mass);
  tree_new->Branch("comb_hadtop_mass", &_comb_hadtop_mass);
  tree_new->Branch("comb_leptop_mass", &_comb_leptop_mass);
  tree_new->Branch("comb_hadWblepTop_mass", &_comb_hadWblepTop_mass);
  tree_new->Branch("comb_minbhadTopqhadW_dR", &_comb_minbhadTopqhadW_dR);
  tree_new->Branch("comb_hadWblepTop_dR", &_comb_hadWblepTop_dR);
  tree_new->Branch("comb_blepTopbhadTop_dR", &_comb_blepTopbhadTop_dR);
  tree_new->Branch("comb_bhadTopq2hadW_dR", &_comb_bhadTopq2hadW_dR);
  tree_new->Branch("comb_hadWbhadTop_dR", &_comb_hadWbhadTop_dR);
  tree_new->Branch("comb_diff_mindRbhadTopqhadW_dRlepblepTop", &_comb_diff_mindRbhadTopqhadW_dRlepblepTop);
  tree_new->Branch("comb_Higgsq1hadW_mass", &_comb_Higgsq1hadW_mass);
  tree_new->Branch("comb_bbHiggs_dR", &_comb_bbHiggs_dR);
  tree_new->Branch("comb_bhadTopq1hadW_dR", &_comb_bhadTopq1hadW_dR);
  tree_new->Branch("comb_qqhadW_dR", &_comb_qqhadW_dR);

  tree_new->Branch("comb_lepbhadTop_dR", &_comb_lepbhadTop_dR);
  tree_new->Branch("comb_lepb1Higgs_dR", &_comb_lepb1Higgs_dR);
  tree_new->Branch("comb_lepWbhadTop_mass", &_comb_lepWbhadTop_mass);
  tree_new->Branch("comb_lepblepTop_dR", &_comb_lepblepTop_dR);

  tree_new->Branch("comb_bh1_index", &_comb_bh1_index);
  tree_new->Branch("comb_bh2_index", &_comb_bh2_index);
  tree_new->Branch("comb_bhadt_index", &_comb_bhadt_index);
  tree_new->Branch("comb_blept_index", &_comb_blept_index);
  tree_new->Branch("comb_qW1_index", &_comb_qW1_index);
  tree_new->Branch("comb_qW2_index", &_comb_qW2_index);
  tree_new->Branch("comb_truthMatchPattern", &_comb_truthMatchPattern);
 
  //
  tree_new->Branch("comb_higgs_truthPt", &_comb_higgs_truthPt);
  tree_new->Branch("comb_best_higgs_Pt", &_comb_best_higgs_Pt);
  tree_new->Branch("comb_best_higgs_Pt_withH", &_comb_best_higgs_Pt_withH);

  tree_new->Branch("higgs_truthPt", &_higgs_truthPt);
  tree_new->Branch("best_higgs_Pt", &_best_higgs_Pt);
  tree_new->Branch("best_higgs_Pt_withH", &_best_higgs_Pt_withH);

  tree_new->Branch("nBTags_DL1r85",&_nBTags_DL1r85);
  tree_new->Branch("nBTags_DL1r70",&_nBTags_DL1r70);

  

  //4c1*3c1*2c2
  //hard-coded combination now:
  
  // combinations:
  
  //[bhadt,blept,h,h]
  //[blept,bhadt,h,h]
  //[h,h,bhadt,blept]
  //[h,h,blept,bhadt]
  //[h,bhadt,h,blept]
  //[h,blept,h,bhadt]
  //[h,bhadt,blept,h]
  //[h,blept,bhadt,h]
  //[bhadt,h,h,blept]
  //[blept,h,h,bhadt]
  //[bhadt,h,blept,h]
  //[blept,h,bhadt,h]
  
  vector<vector<int> > comb{ {0,1,2,3,4,5},
                             {1,0,2,3,4,5},
	                     {2,3,0,1,4,5},
	                     {3,2,0,1,4,5},
	                     {1,3,0,2,4,5},
	                     {3,1,0,2,4,5},
	             	     {1,2,0,3,4,5},
		             {2,1,0,3,4,5},
		             {0,3,1,2,4,5},
		             {3,0,1,2,4,5},
			     {0,2,1,3,4,5},
			     {2,0,1,3,4,5} };

  static const int bhadt = 0;
  static const int blept = 1;
  static const int bh1 = 2;
  static const int bh2 = 3;
  static const int qW1 = 4;
  static const int qW2 = 5;

  struct Jet{
    TLorentzVector tlv;
    int tagWeightBin;
    int index;
  };

  vector<Jet> jets;
  vector<Jet> jet_comb;
  //vector<int> jets;
  //vector<int> jet_comb;

  TLorentzVector lep;
  vector<TLorentzVector> nusp4;
   
  for(int i=0;i<nentries;i++){
  //for(int i=0;i<100;i++){
    if(i%1000==0)
       cout<<"i="<<i<<endl;
    
     _nJets = -999;
     _nElectrons = -999;
     _nMuons = -999;
     _nBTags_DL1r_85 = -999;
     _nBTags_DL1r_70 = -999;
     _eventNumber=-999;
     _weight_mc = -999;
     _met_met = -999 ;
     _met_phi = -999 ;
     _jet_pt = 0 ;
     _jet_eta = 0 ;
     _jet_phi = 0 ;
     _jet_e = 0 ;
     _el_pt = 0 ;
     _el_eta = 0 ;
     _el_phi = 0 ;
     _el_e = 0 ;
     _el_charge = 0 ;
     _mu_pt = 0 ;
     _mu_eta = 0 ;
     _mu_phi = 0 ;
     _mu_e = 0 ;
     _mu_charge = 0 ;
     _jet_truthmatch = 0;
     _jet_tagWeightBin = 0;
     _jet_isbtagged_DL1r_85=0;
     _jet_isbtagged_DL1r_70=0;
     
     _truth_higgs_pt = -999;
     _TTHReco_withH_best_Higgs_pt = -999;
     _TTHReco_best_Higgs_pt = -999;
     
     ///new variables//

     _comb_higgs_mass.clear();
     _comb_higgs_pt.clear();
     _comb_hadW_mass.clear();
     _comb_hadtop_mass.clear();
     _comb_leptop_mass.clear();
     _comb_hadWblepTop_mass.clear();
     _comb_minbhadTopqhadW_dR.clear();
     _comb_hadWblepTop_dR.clear();
     _comb_blepTopbhadTop_dR.clear();
     _comb_bhadTopq2hadW_dR.clear();
     _comb_hadWbhadTop_dR.clear();
     _comb_diff_mindRbhadTopqhadW_dRlepblepTop.clear();
     _comb_Higgsq1hadW_mass.clear();
     _comb_bbHiggs_dR.clear();
     _comb_bhadTopq1hadW_dR.clear();
     _comb_qqhadW_dR.clear();
     _comb_lepbhadTop_dR.clear();
     _comb_lepb1Higgs_dR.clear();
     _comb_lepWbhadTop_mass.clear();
     _comb_lepblepTop_dR.clear();

     _comb_bh1_index.clear();
     _comb_bh2_index.clear();
     _comb_bhadt_index.clear();
     _comb_blept_index.clear();
     _comb_qW1_index.clear();
     _comb_qW2_index.clear();
     _comb_truthMatchPattern.clear();
         
     _comb_higgs_truthPt.clear();	 
     _comb_best_higgs_Pt.clear();	 
     _comb_best_higgs_Pt_withH.clear();

     _higgs_truthPt = -999;
     _best_higgs_Pt = -999;
     _best_higgs_Pt_withH = -999;
     
     _nBTags_DL1r85 = -999;
     _nBTags_DL1r70 = -999;

     _weight_mc_new = -999;
     _comb_weight_mc.clear();
     _comb_truthpt_rewg.clear();
     _truthpt_weights = -999;

     nominal_Loose->GetEntry(i);

     jets.clear();
     jet_comb.clear();
     nusp4.clear();
     lep.SetPtEtaPhiE(0,0,0,0);
     
     int number_nu_sol = 0;

     //To check if the two cases are really exclusive
     if(_nElectrons>0) lep.SetPtEtaPhiE((*_el_pt)[0],(*_el_eta)[0],(*_el_phi)[0],(*_el_e)[0]);
     else if(_nMuons>0) lep.SetPtEtaPhiE((*_mu_pt)[0],(*_mu_eta)[0],(*_mu_phi)[0],(*_mu_e)[0]);
     
     //neutrino
     float Wmass = 80385.;
     float neut_pt  = _met_met;
     float neut_phi = _met_phi;
     float neut_px  = neut_pt*TMath::Cos(neut_phi);
     float neut_py  = neut_pt*TMath::Sin(neut_phi);

     float Beta = Wmass*Wmass - lep.M()*lep.M() +2*lep.Px()*neut_px + 2*lep.Py()*neut_py;
     float Delta = lep.E()*lep.E()*(Beta*Beta + (2*lep.Pz()*neut_pt)*(2*lep.Pz()*neut_pt)-(2*lep.E()*neut_pt)*(2*lep.E()*neut_pt));

     if(Delta<=0){
       float neut_pz = 0.5*(lep.Pz()*Beta)/(lep.E()*lep.E() - lep.Pz()*lep.Pz());
       TLorentzVector nup4;
       nup4.SetXYZM(neut_px, neut_py, neut_pz, 0.);
       nusp4.push_back(nup4);
       number_nu_sol = 1;
     }
     else{
       float neut_pz1 = 0.5*(lep.Pz()*Beta + TMath::Sqrt(Delta))/(lep.E()*lep.E() - lep.Pz()*lep.Pz());
       float neut_pz2 = 0.5*(lep.Pz()*Beta - TMath::Sqrt(Delta))/(lep.E()*lep.E() - lep.Pz()*lep.Pz());
       TLorentzVector nup41;
       TLorentzVector nup42;
       nup41.SetXYZM(neut_px, neut_py, neut_pz1, 0);
       nup42.SetXYZM(neut_px, neut_py, neut_pz2, 0);
       nusp4.push_back(nup41);
       nusp4.push_back(nup42);
       number_nu_sol = 2;
     }

      _higgs_truthPt = _truth_higgs_pt;
      _best_higgs_Pt = _TTHReco_best_Higgs_pt;
      _best_higgs_Pt_withH = _TTHReco_withH_best_Higgs_pt;
     
      _nBTags_DL1r85 = _nBTags_DL1r_85;
      _nBTags_DL1r70 = _nBTags_DL1r_70;
      _weight_mc_new = _weight_mc;
  
      //_truthpt_weights = 1/_truth_higgs_pt;
      hist_truthPt->Fill(_truth_higgs_pt);
      //cout<<hist_truthPt->GetBinContent(i)<<"\n";
      //hist_truthPt_inv->Fill(1);
      //double c1 = hist_truthPt->GetBinContent(i);
      //if (c1 > 0){
      //hist_truthPt->SetBinContent(i, _truth_higgs_pt);
      hist_truthPt_inv->SetBinContent(i, 1);
      //cout<<hist_truthPt_inv->GetBinContent(i)<<"\n";
      //}
      //hist_truthPt_inv->Fill(1);
      //hist_truthPt_inv->SetBinContent(1);
      //TH1F* h_ratio  = (TH1F*)hist_truthPt->Clone("h_ratio");
      //h_ratio->Divide(hist_truthPt_inv);
     
      //int bin_wg = h_ratio->FindBin(_truth_higgs_pt);
      //double c2 = hist_truthPt->GetBinContent(i);
      //_truthpt_weights = 1/c2;

     for(unsigned int i_jet=0; i_jet<(*_jet_pt).size(); i_jet++){

       Jet jet;
       jet.tlv.SetPtEtaPhiE((*_jet_pt)[i_jet],(*_jet_eta)[i_jet],(*_jet_phi)[i_jet],(*_jet_e)[i_jet]);
       jet.index = i_jet;
       jet.tagWeightBin = (*_jet_tagWeightBin)[i_jet];

       jets.push_back(jet);

     }

     std::sort(jets.begin(), jets.end(),
	       [&](Jet jet1, Jet jet2){
		 int a = jet1.tagWeightBin;
		 int b = jet2.tagWeightBin;
		 if( a == b ) return jet1.tlv.Pt() > jet2.tlv.Pt();
		 return a>b;
	       }
	       );
     
     float higgs_mass = 0.;
     float higgs_pt = 0.;
     float hadW_mass = 0;
     float hadtop_mass = 0;
     float leptop_mass = 0;
     float hadWblepTop_mass = 0;
     float minbhadTopqhadW_dR = 0;
     float hadWblepTop_dR = 0;
     float blepTopbhadTop_dR = 0;
     float bhadTopq2hadW_dR = 0;
     float hadWbhadTop_dR = 0;
     float diff_mindRbhadTopqhadW_dRlepblepTop = 0;
     float Higgsq1hadW_mass = 0;
     float bbHiggs_dR = 0;
     float bhadTopq1hadW_dR = 0;
     float qqhadW_dR = 0;
     float lepbhadTop_dR = 0;
     float lepb1Higgs_dR = 0;
     float lepWbhadTop_mass = 0;
     float lepblepTop_dR = 0;

     int bh1_index = -1;
     int bh2_index = -1;
     int bhadt_index = -1;
     int blept_index = -1;
     int qW1_index = -1;
     int qW2_index = -1;

     int truthMatchPattern = 0;

     bool is_new_jetcomb = true;
     
     //if(jets.size()>=6){
     if(jets.size()>=6 && (*_jet_isbtagged_DL1r_85)[jets[3].index]){

       // 4 leading jets ordered by tagWeightBin as b-jets + 2 next jets as light-jets
       jet_comb = jets;
       jet_comb.resize(6);

       for(unsigned int q=0;q<comb.size();q++){ //combination index: q

	 if(is_new_jetcomb){

	   TLorentzVector bh1_jet = jet_comb[comb[q][bh1]].tlv;
	   TLorentzVector bh2_jet = jet_comb[comb[q][bh2]].tlv;
	   TLorentzVector bhadt_jet = jet_comb[comb[q][bhadt]].tlv;
	   TLorentzVector blept_jet = jet_comb[comb[q][blept]].tlv;
	   TLorentzVector qW1_jet = jet_comb[comb[q][qW1]].tlv;
	   TLorentzVector qW2_jet = jet_comb[comb[q][qW2]].tlv;

	   TLorentzVector Higgs = bh1_jet+bh2_jet;
	   TLorentzVector HadW = qW1_jet+qW2_jet;
           TLorentzVector HadTop = bhadt_jet+HadW;
	   

	   higgs_mass = Higgs.M();
	   higgs_pt = Higgs.Pt();
	   hadtop_mass = HadTop.M();
	   hadW_mass = HadW.M();
	   hadWblepTop_mass = (blept_jet+HadW).M();
	   minbhadTopqhadW_dR = std::min(bhadt_jet.DeltaR(qW1_jet),
				 bhadt_jet.DeltaR(qW2_jet));
	   hadWblepTop_dR = blept_jet.DeltaR(HadW);
	   blepTopbhadTop_dR = blept_jet.DeltaR(bhadt_jet);
	   bhadTopq2hadW_dR = bhadt_jet.DeltaR(qW2_jet);
	   hadWbhadTop_dR = bhadt_jet.DeltaR(HadW);
	   Higgsq1hadW_mass = (Higgs+qW1_jet).M();
	   bbHiggs_dR = bh1_jet.DeltaR(bh2_jet);
	   bhadTopq1hadW_dR = bhadt_jet.DeltaR(qW1_jet);
	   qqhadW_dR = qW1_jet.DeltaR(qW2_jet);

	   lepblepTop_dR = lep.DeltaR(blept_jet);
	   lepbhadTop_dR = lep.DeltaR(bhadt_jet);
	   lepb1Higgs_dR = lep.DeltaR(bh1_jet);
	   diff_mindRbhadTopqhadW_dRlepblepTop = std::min(bhadt_jet.DeltaR(qW1_jet),
						  bhadt_jet.DeltaR(qW2_jet))
	                                         - lep.DeltaR(blept_jet);

	   TLorentzVector LepTop = blept_jet+lep+nusp4[0];
	   leptop_mass = LepTop.M();
	   lepWbhadTop_mass = (lep+bhadt_jet+nusp4[0]).M();

	   bh1_index = jet_comb[comb[q][bh1]].index;
	   bh2_index = jet_comb[comb[q][bh2]].index;
	   bhadt_index = jet_comb[comb[q][bhadt]].index;
	   blept_index = jet_comb[comb[q][blept]].index;
	   qW1_index = jet_comb[comb[q][qW1]].index;
	   qW2_index = jet_comb[comb[q][qW2]].index;

	   truthMatchPattern = 0;

	   //Also consider good matching if Higgs jets are swapped
	   if((*_jet_truthmatch)[bh1_index] & (1<<JetMatchLeadingBHiggs)) truthMatchPattern |= (1<<JetMatchLeadingBHiggs);
	   else if((*_jet_truthmatch)[bh1_index] & (1<<JetMatchSubLeadingBHiggs)) truthMatchPattern |= (1<<JetMatchSubLeadingBHiggs);
	   if((*_jet_truthmatch)[bh2_index] & (1<<JetMatchSubLeadingBHiggs)) truthMatchPattern |= (1<<JetMatchSubLeadingBHiggs);
	   else if((*_jet_truthmatch)[bh2_index] & (1<<JetMatchLeadingBHiggs)) truthMatchPattern |= (1<<JetMatchLeadingBHiggs);

	   if((*_jet_truthmatch)[bhadt_index] & (1<<JetMatchBHadTop)) truthMatchPattern |= (1<<JetMatchBHadTop);
	   if((*_jet_truthmatch)[blept_index] & (1<<JetMatchBLepTop)) truthMatchPattern |= (1<<JetMatchBLepTop);

	   //Also consider good matching if W jets are swapped
	   if((*_jet_truthmatch)[qW1_index] & (1<<JetMatchLeadingJetW)) truthMatchPattern |= (1<<JetMatchLeadingJetW);
	   else if((*_jet_truthmatch)[qW1_index] & (1<<JetMatchSubLeadingJetW)) truthMatchPattern |= (1<<JetMatchSubLeadingJetW);
	   if((*_jet_truthmatch)[qW2_index] & (1<<JetMatchSubLeadingJetW)) truthMatchPattern |= (1<<JetMatchSubLeadingJetW);
	   else if((*_jet_truthmatch)[qW2_index] & (1<<JetMatchLeadingJetW)) truthMatchPattern |= (1<<JetMatchLeadingJetW);

	   if(number_nu_sol==2){
	     is_new_jetcomb = false;
	     q--;
	   }

	 }

	 else{
	   TLorentzVector blept_jet = jet_comb[comb[q][blept]].tlv;
	   TLorentzVector bhadt_jet = jet_comb[comb[q][bhadt]].tlv;
	   TLorentzVector LepTop = blept_jet+lep+nusp4[1];
	   leptop_mass = LepTop.M();
	   lepWbhadTop_mass = (lep+bhadt_jet+nusp4[1]).M();

	   is_new_jetcomb = true;
	 }

	 _comb_higgs_mass.push_back(higgs_mass);
	 _comb_higgs_pt.push_back(higgs_pt);
	 _comb_hadtop_mass.push_back(hadtop_mass);
	 _comb_hadW_mass.push_back(hadW_mass);
	 _comb_hadWblepTop_mass.push_back(hadWblepTop_mass);
	 _comb_minbhadTopqhadW_dR.push_back(minbhadTopqhadW_dR);
	 _comb_hadWblepTop_dR.push_back(hadWblepTop_dR);
	 _comb_blepTopbhadTop_dR.push_back(blepTopbhadTop_dR);
	 _comb_bhadTopq2hadW_dR.push_back(bhadTopq2hadW_dR);
	 _comb_hadWbhadTop_dR.push_back(hadWbhadTop_dR);
	 _comb_Higgsq1hadW_mass.push_back(Higgsq1hadW_mass);
	 _comb_bbHiggs_dR.push_back(bbHiggs_dR);
	 _comb_bhadTopq1hadW_dR.push_back(bhadTopq1hadW_dR);
	 _comb_qqhadW_dR.push_back(qqhadW_dR);
	 _comb_leptop_mass.push_back(leptop_mass);
	 _comb_lepWbhadTop_mass.push_back(lepWbhadTop_mass);
	 _comb_lepbhadTop_dR.push_back(lepbhadTop_dR);
	 _comb_lepb1Higgs_dR.push_back(lepb1Higgs_dR);
	 _comb_lepblepTop_dR.push_back(lepblepTop_dR);
	 _comb_diff_mindRbhadTopqhadW_dRlepblepTop.push_back(diff_mindRbhadTopqhadW_dRlepblepTop);

	 _comb_bh1_index.push_back(bh1_index);
	 _comb_bh2_index.push_back(bh2_index);
	 _comb_bhadt_index.push_back(bhadt_index);
	 _comb_blept_index.push_back(blept_index);
	 _comb_qW1_index.push_back(qW1_index);
	 _comb_qW2_index.push_back(qW2_index);
	 _comb_truthMatchPattern.push_back(truthMatchPattern);

         //
	 _comb_higgs_truthPt.push_back(_higgs_truthPt);
	 _comb_best_higgs_Pt.push_back(_TTHReco_best_Higgs_pt);
         _comb_best_higgs_Pt_withH.push_back(_TTHReco_withH_best_Higgs_pt);
         _comb_weight_mc.push_back(_weight_mc);
         _comb_truthpt_rewg.push_back(1/_truth_higgs_pt);
       }

     }

     //if(i<1e6){
     tree_new->Fill();
     
    
  }

   
  //for(int j=0 ; j <= hist_truthPt->GetNbinsX() ; j++){
    //double c2 = hist_truthPt->GetBinContent(j);
    //_truthpt_weights = 1/c2;
	  
  //hist_truthPt_inv->SetBinContent(j, 1);
  //}

  //tree_new->AddFriend(nominal_Loose);
  f_new->cd();
  tree_new->Write();
  hist_truthPt->Write();
  hist_truthPt_inv->Write();
  f_new->Close();

  return;

}


