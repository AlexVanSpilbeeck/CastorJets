////////////////////////////////////////
//////// New CMSSW_4_2_X version ///////
////////////////////////////////////////

/* Updated by Alex Van Spilbeeck.

 This code looks at:
	-> Gen Jets (-6.1 < eta < -5.7), E above Ethreshold, pT > 1 GeV
	-> Castor Jets, E > 0 GeV, pT > 1 GeV
	-> Gen jets are matched to closest detector level jet in \Delta phi. 
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "JetAnalyzer_radii.h"
#include "HistoRetriever.h"

//STANDARD ROOT INCLUDES

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TBranch.h>
#include <TRandom3.h>
#include <TThread.h>
#include <TStopwatch.h>
#include <TBranchElement.h>

//STANDARD C++ INCLUDES
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>

// own classes includes
#include "../src/MyCastorRecHit.h"
#include "../src/MyCastorDigi.h"
#include "../src/MyCastorTower.h"
#include "../src/MyCastorJet.h"
#include "../src/MyEvtId.h"
#include "../src/MyGenKin.h"
#include "../src/MyDiJet.h"
#include "../src/MyVertex.h"
#include "../src/MyHLTrig.h"
#include "../src/MyL1Trig.h"
#include "../src/MyJet.h"
#include "../src/MyBeamSpot.h"
#include "../src/MyGenPart.h"
#include "../src/MyCaloTower.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

//Fastjet
//#include "fastjet/PseudoJet.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/PseudoJet.hh"
//using namespace fastjet;

#define jetPtThreshold 35.
#define jetEThreshold_gen 0. 
#define jetEThreshold_det 0. 
#define EbinWidth 5.
#define EbinWidth_rel 1.4

#define jet_distance 10
#define jet_distance_string "JER_allPairs__JetSorted_ak5"
//#define gen_radius_ "ak7"
//#define det_radius_ "ak7"
#define GenJetContained 0.5
//#define totalEvents_ 10000
#define PI 3.14159265359


#include "../../../RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../../../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#endif



TFile *JetAnalyzer_radii::currentStaticTFile_ = new TFile();

JetAnalyzer_radii::JetAnalyzer_radii(TString inputdir, TObjArray* filelist, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date) {
    
	std::cout << "constructing JetAnalyzer_radii class..." << std::endl;
	
    inputdir_ = inputdir;
    filelist_ = filelist;
    isData_ = isData;
    outputname_ = outputname;
    gen_radius_ = gen_radius;
    det_radius_ = det_radius;
    totalEvents_ = totalEvents;
    date_ = date;	
	std::cout << "initializing basic variables..." << std::endl;
    
    // initialize basic variables

    
    LoopOutputFile_ = "";

	currentfile_ = "";
	currentTFile_ = new TFile();
	
	std::cout << "all initialisations done, class constructed!" << std::endl;
    
}

JetAnalyzer_radii::~JetAnalyzer_radii() { }

void JetAnalyzer_radii::Loop() {

#ifdef __CINT__
  gSystem->Load("../../../RooUnfold-1.1.1/libRooUnfold.so");
#endif
	
	
	std::cout << " JetAnalyzer_radii Loop function is started " << std::endl;
	
	TString tstring = outputname_;
	std::cout << " TString outputname_ = " << tstring << std::endl;
	TString string_tag = jet_distance_string;
	TString string_det_radius = det_radius_;
        TString string_gen_radius;	if (!isData_){ string_gen_radius = gen_radius_; }
 
	
    // reweight the MC in this case
    bool reweightMC = false;
    if (tstring.Contains("Reweighted")) {
        std::cout << "We will reweight the MC now !!!" << std::endl;
        reweightMC = true;
    }
    
	
	using namespace std;
	int it = 0;
	int totalevents = 0;
	
	/////////////////////////////////////
	// Define all histograms
	/////////////////////////////////////


	int Ebins_fix = 100.;
	double Emin_fix = -100.;
	double Emax_fix = 3000.;
	
        // AVS - number of bins.
	/* We get binwidth from peak JER. */
	double Emin = jetEThreshold_det + 100., Emax = 3000.;
	int Ebins = static_cast<int> ( (Emax - Emin)/EbinWidth );


	/* If we want variable Ebins. */
	double current_lowE = 200.;
	int count_bincenters = 3;

	vector<double> binEdges;
	binEdges.push_back( -50.);
	binEdges.push_back( 0. );
	binEdges.push_back( 50. );
	binEdges.push_back( jetEThreshold_det );
	binEdges.push_back( current_lowE );
	while( current_lowE/EbinWidth_rel < Emax ){

	  binEdges.push_back( current_lowE * EbinWidth_rel );

 	  cout << "\tLower edge\t" << count_bincenters << "\t" << current_lowE << "\tBinwidth\t" << current_lowE *0.2 << endl;

	  current_lowE += current_lowE*0.2;
	  count_bincenters++;
	}

	Ebins = binEdges.size();
        float *Ebins_var = new float[Ebins];
//float Ebins_var[Ebins];
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_var[i] = binEdges[i];

	  cout << "\t\tBin " << i << "\tout of\t" << sizeof(Ebins_var) << " or " << binEdges.size() << " or " << Ebins << "\t" << Ebins_var[i] << "\t" << binEdges[i] << endl;
	}

	// We need an extended bin range for our matrix + hits & misses
	// Let us take bins 0 - 50, 50 - 100
	float *Ebins_ext = new float[Ebins+2];
	Ebins_ext[0] = 0.;
	Ebins_ext[1] = 50.;
	for(int i = 0; i < binEdges.size(); i++){
	  Ebins_ext[i+2] = binEdges[i];
	}

        /* End variable ebins. */

cout << "Variable bins done" << endl;
	// detector level histograms
	
	char name [100];
	char title [100];
		
	TH1D *hCASTORTowerMulti = new TH1D("hCASTORTowerMulti","CASTOR Tower Multiplicity (N above threshold) distribution",17,0,17);
	hCASTORTowerMulti->Sumw2();
	
	// default energy flow histos - using 5 modules
	TH1D *hCASTOReflow = new TH1D("hCASTOReflow","Total CASTOR energy flow in first 5 modules",252,-30,3750);
	hCASTOReflow->Sumw2();
	
	TH2D *h2CASTOReflow_grid = new TH2D("h2CASTOReflow_grid","CASTOR energy weighted module vs sector distribution",16,1,17,14,1,15);
	
	TH1D *hCASTOReflow_channel[224];
	for (int i=0;i<224;i++) {
		sprintf(name,"hCASTOReflow_channel_%d",i+1);
		sprintf(title,"CASTOR Energy distribution for channel %d",i+1);
		hCASTOReflow_channel[i] = new TH1D(name,title,100,-10,1800);
		hCASTOReflow_channel[i]->Sumw2();
	}
	
	// Castor jet histograms
//	TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",Ebins, Emin, Emax);
        TH1D *hCastorJet_energy = new TH1D("hCastorJet_energy","CastorJet energy distribution",200.,100.,3000.);
        TH1D *hGenJet_energy = new TH1D("hGenJet_energy","GenJet energy distribution",Ebins, Emin, Emax);
	TH1D *hCastorJet_pt = new TH1D("hCastorJet_pt","CastorJet pt distribution",30,0,30);
	TH1D *hCastorJet_em = new TH1D("hCastorJet_em","CastorJet EM energy distribution",150,0,1500);
	TH1D *hCastorJet_had = new TH1D("hCastorJet_had","CastorJet HAD energy distribution",150,0,1500);
	TH1D *hCastorJet_fem = new TH1D("hCastorJet_fem","CastorJet EM/(EM+HAD) distribution",60,-0.1,1.1);
	TH1D *hCastorJet_fhot = new TH1D("hCastorJet_fhot","CastorJet Fhot distribution",60,-0.1,1.1);
	TH1D *hCastorJet_width = new TH1D("hCastorJet_width","CastorJet width distribution",100,0,1);
	TH1D *hCastorJet_depth = new TH1D("hCastorJet_depth","CastorJet depth distribution",100,-16000,-14000);
	TH1D *hCastorJet_sigmaz = new TH1D("hCastorJet_sigmaz","CastorJet sigmaz distribution",100,0,500);
	TH1D *hCastorJet_ntower = new TH1D("hCastorJet_ntower","CastorJet ntower distribution",16,1,17);
	TH1D *hCastorJet_eta = new TH1D("hCastorJet_eta","CastorJet eta distribution",14,-6.6,-5.2);
	TH1D *hCastorJet_phi = new TH1D("hCastorJet_phi","CastorJet phi distribution",16,-M_PI,+M_PI);
	TH1D *hCastorJet_multi = new TH1D("hCastorJet_multi","CastorJet multiplicity distribution",17,0,17);

	// Central-forward configuration response matrix.
        TH1D *hCastorJet_cf_energy = new TH1D("hCastorJet_cf_energy", "CastorJet energy in with leading central jet",Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_gen = new TH1D("hCastorJet_cf_energy_gen", "CastorJet energy in with leading central jet",Ebins+10,Emin-10*EbinWidth,Emax);
	TH2D *hCastorJet_cf_energy_response = new TH2D("hCastorJet_cf_energy_response", "CastorJet energy response matrix for leading central jet",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_fakes = new TH1D("hCastorJet_cf_energy_fakes", "CastorJet energy in with leading central jet - fakes",Ebins+10,Emin-10*EbinWidth,Emax);
        TH1D *hCastorJet_cf_energy_misses = new TH1D("hCastorJet_cf_energy_misses", "CastorJet energy in with leading central jet - misses",Ebins+10,Emin-10*EbinWidth,Emax);

	// All Casto jets response matrix.
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution - Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);	// ,Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_fakes 	 = new TH2D("hCastorJet_energy_fakes",	 "CastorJet energy distribution - Fakes", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_energy_misses	 = new TH2D("hCastorJet_energy_misses",	 "CastorJet energy distribution - Misses", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);	// ,Ebins+10,Emin-10*EbinWidth,Emax);
	
	TH2D *hCastorJet_energy_complete_response = new TH2D("hCastorJet_energy_complete_response","CastorJet energy distribution - Complete Response", Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);

	// Response with variable binning.
/*
        TH2D *hCastorJet_energy_response = new TH2D("hCastorJet_energy_response","CastorJet energy distribution  - Response",Ebins-1, Ebins_var,Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_fakes = new TH1D("hCastorJet_energy_fakes","CastorJet energy distribution - Fakes", Ebins-1, Ebins_var);
        TH1D *hCastorJet_energy_misses = new TH1D("hCastorJet_energy_misses","CastorJet energy distribution - Misses",Ebins-1, Ebins_var);
*/
//	TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins-1, Ebins_var);
        TH2D *hCastorJet_energy_ratio = new TH2D("hCastorJet_energy_ratio", "Ratio of generator to detector", 1000, 0., 10.,Ebins, Emin, Emax);

	TH2D *hCastorJet_cf_Matrix = new TH2D("hCastorJet_cf_responseMatrix", "Central-forward jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
//        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins+10,Emin-10*EbinWidth,Emax,Ebins+10,Emin-10*EbinWidth,Emax);
        TH2D *hCastorJet_Matrix = new TH2D("hCastorJet_responseMatrix", "Castor jet Response Matrix",Ebins-1, Ebins_var,Ebins-1, Ebins_var);

	  TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)"	,Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);
          TH2D *hCastorJet_Matrix_had_e = new TH2D("hCastorJet_responseMatrix_had_e", 	"Castor jet Response Matrix (had-em)"	,Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);
          TH2D *hCastorJet_Matrix_em_pi = new TH2D("hCastorJet_responseMatrix_em_pi", 	"Castor jet Response Matrix (em-had)"	,Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);
          TH2D *hCastorJet_Matrix_em_e 	= new TH2D("hCastorJet_responseMatrix_em_e", 	"Castor jet Response Matrix (em-em)"	,Ebins_fix, Emin_fix, Emax_fix, Ebins_fix, Emin_fix, Emax_fix);
//          TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)",Ebins-1, Ebins_var,Ebins-1, Ebins_var);
//          TH2D *hCastorJet_Matrix_had_pi = new TH2D("hCastorJet_responseMatrix_had_pi", "Castor jet Response Matrix (had-had)",Ebins-1, Ebins_var,Ebins-1, Ebins_var);


	// Number of trackjets.
	TH2D *hTrackjets_2D_number = new TH2D("hTrackjets_2D_number", "Trackjets vs. Charged Gen Jet", 17,0.,17., 17,0.,17.);
        TH2D *hTrackjets_2D_pt = new TH2D("hTrackjets_2D_pt", "Trackjets vs. Charged Gen Jet", 30,0.,30., 30,0.,30.);

        TH1D *hCastorJet_energy_gen = new TH1D("hCastorJet_energy_gen","CastorJet energy distribution",Ebins-1, Ebins_var);
        TH1D *hCastorJet_pt_gen = new TH1D("hCastorJet_pt_gen","CastorJet pt distribution",500,0,10.);

	// Efficiency control.
	TH1D *hJER = new TH1D("hJER", "Jet Energy Resolution", 200, -5, 5);
        TH2D *hJER_per_energy = new TH2D("hJER_per_energy", "#DeltaE/E for fixed energies;E_{gen};JER", 	200.,100.,3000., 200, -5, 5);
	     TH2D *hJER_per_energy_had_pi = new TH2D("hJER_per_energy_had_pi", "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_energy_had_e  = new TH2D("hJER_per_energy_had_e",  "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_energy_em_pi  = new TH2D("hJER_per_energy_em_pi",  "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_energy_em_e	  = new TH2D("hJER_per_energy_em_e",   "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);

	     TH2D *hJER_per_energy_had_det= new TH2D("hJER_per_energy_had_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_energy_em_det= new TH2D("hJER_per_energy_em_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_energy_none_det= new TH2D("hJER_per_energy_none_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);

        TH2D *hJER_per_eDet = new TH2D("hJER_per_eDet", "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_had_pi = new TH2D("hJER_per_eDet_had_pi", "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_had_e  = new TH2D("hJER_per_eDet_had_e",  "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_em_pi  = new TH2D("hJER_per_eDet_em_pi",  "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_em_e   = new TH2D("hJER_per_eDet_em_e",   "#DeltaE/E for fixed energies;E_{gen};JER",     200.,100.,3000., 200, -5, 5);

             TH2D *hJER_per_eDet_had_det= new TH2D("hJER_per_eDet_had_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_em_det= new TH2D("hJER_per_eDet_em_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);
             TH2D *hJER_per_eDet_none_det= new TH2D("hJER_per_eDet_none_det", "#DeltaE/E for fixed energies;E_{gen};JER",   200.,100.,3000., 200, -5, 5);



        TH2D *hJER_per_distance = new TH2D("hJER_per_distance", "#DeltaE/E for distance;#DeltaR;JER", 	200, 0., 6.5, 
															200, -5, 5);
        TH2D *hJER_per_eta = new TH2D("hJER_per_eta", "#DeltaE/E for distance;#DeltaR;JER",               14,-6.6,-5.2,
                                                                                                                        200, -5, 5);
        TH2D *hEnergy_per_eta = new TH2D("hEnergy_per_eta", "E_{gen} vs. #eta of leading jet", 200, 100., 3000.,  14,-6.6,-5.2);

        // Matches.
        TH1D *hMatched = new TH1D("hMatched", "Castor and Gen jets match", 10, -0.5, 9.5);
        TH1D *hUnmatched = new TH1D("hUnmatched", "Castor and Gen jets don't match", Ebins, Emin, Emax);
        TH1D *hJRE = new TH1D("hJRE", "Jet Reconstruction Efficiency", Ebins, Emin, Emax);

        // JER
        TH1D *hJER_all = new TH1D("hJER_all", TString::Format("Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH2D *hJER_per_energy_all = new TH2D("hJER_per_energy_all", "#DeltaE/E for fixed energies;E_{gen};JER",	200.,100.,3000., 	200, -5, 5);
        TH2D *hJER_per_distance_all = new TH2D("hJER_per_distance_all", "#DeltaE/E for distance;#DeltaR;JER", 	200, 0., 6.5, 		200, -5, 5);
	TH2D *hJER_per_eta_all = new TH2D("hJER_per_eta_all", "#DeltaE/E for distance;#DeltaR;JER", 		14,	-6.6,-5.2, 	200, -5, 5);
	TH2D *hEnergy_per_eta_all = new TH2D("hEnergy_per_eta_all", "E_{gen} vs. #eta of leading jet", 				200, 100., 3000.,  	14,-6.6,-5.2);
	
        TH1D *hJER_had_pi = new TH1D("hJER_had_pi", TString::Format("Hadronic (gen and det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_had_e  = new TH1D("hJER_had_e", TString::Format("Hadronic (gen) and em (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_em_pi = new TH1D("hJER_em_pi", TString::Format("EM (gen) and hadronic (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_em_e  = new TH1D("hJER_em_e", TString::Format("EM (gen and det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_both_pi = new TH1D("hJER_both_pi", TString::Format("Hybrid (gen) and hadronic (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);
        TH1D *hJER_both_e  = new TH1D("hJER_both_e", TString::Format("Hybrid (gen) and em (det) Jet Energy Resolution, " + string_gen_radius + " " + string_det_radius ), 200, -5, 5);


//        TMap *JER_map  JER_map;
/*	JER_map["had_had"]  = *hJER_had_pi;
        JER_map["had_em"]   = *hJER_had_e;
        JER_map["em_had"]   = *hJER_em_pi;
        JER_map["em_em"]    = *hJER_em_e;
        JER_map["both_had"] = *hJER_both_pi;
        JER_map["both_em"]  = *hJER_both_e;	
*/
	// RooUnfold.

//	RooUnfoldResponse response (Ebins, Emin, Emax);
	RooUnfoldResponse response 		(hCastorJet_energy, hCastorJet_energy, "response");
	RooUnfoldResponse response_all		(hCastorJet_energy, hCastorJet_energy, "response_all");

//        RooUnfoldResponse response_cf (Ebins, Emin, Emax, "response_cf");
//        RooUnfoldResponse response_cf_onlyMatches (Ebins, Emin, Emax, "response_cf_onlyMatches");

	// Jet distance distribution.
	TH1D *hDistance 	= new TH1D("hJetDistance",	"Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
	TH1D *hPhiDiff 		= new TH1D("hPhiDiff", 		"Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
        TH1D *hEtaDiff 		= new TH1D("hEtaDiff", 		"Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, 0., 0.8);
	TH2D *hEtaPhiDiff 	= new TH2D("hEtaPhiDiff", 	"Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
        TH2D *hEtaRDiff 	= new TH2D("hEtaRDiff", 	"Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
        TH2D *hPhiRDiff 	= new TH2D("hPhiRDiff", 	"Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

        TH1D *hDistance_all 	= new TH1D("hJetDistance_all", 	"Distance between matched jets;#DeltaR;dN/d#DeltaR", 200, 0., 2.);
        TH1D *hPhiDiff_all 	= new TH1D("hPhiDiff_all", 	"Distance in #varphi;#Delta#varphi;dN/d#Delta#varphi", 200, 0., 3.15);
        TH1D *hEtaDiff_all 	= new TH1D("hEtaDiff_all", 	"Distance in #eta;#Delta#eta;dN/d#Delta#eta", 200, 0., 0.8);
        TH2D *hEtaPhiDiff_all 	= new TH2D("hEtaPhiDiff_all", 	"Distance in #eta and #varphi;#eta;#varphi",200,0.,0.8,200,0.,3.15);
        TH2D *hEtaRDiff_all 	= new TH2D("hEtaRDiff_all", 	"Distance in #eta and R;#Delta#eta;#DeltaR",200,0.,0.8,200,0.,3.3);
        TH2D *hPhiRDiff_all 	= new TH2D("hPhiRDiff_all", 	"Distance in #eta and #varphi;#Delta#varphi;#DeltaR",200,0.,3.3,200,0.,3.3);

	// Pion to electron ratio.
	TH1D *hElectron_energy	= new TH1D("hElectron_energy",	"Energy of electron jets;E_{e} (GeV)",	20, 0., 1000.);
        TH1D *hPion_energy  	= new TH1D("hPion_energy",  	"Energy of Pion jets;E_{#pi} (GeV)",   	20, 0., 1000.);
	TH1D *hPi_e_ratio	= new TH1D("hPi_e_ratio",	"Ratio of pions to electrons;E (GeV)",	20, 0., 1000.);
	
	// Study the jets' energy versus the number of jets.
        TH2D *hNjet_vs_Ejets_gen = new TH2D("hNjet_vs_Ejets_gen", "Number of jets versus E leading jet (gen)",        10, -0.5, 9.5, 50, 100., 3500.);
	TH2D *hNjet_vs_Ejets_det = new TH2D("hNjet_vs_Ejets_det", "Number of jets versus E leading jet (det)", 	10, -0.5, 9.5, 50, 100., 3500.);								
	// Valid gen jets versus Castor jets.
	TH2D *hNumber_of_match_jets = new TH2D("hNumber_of_match_jets", "Castor jet versus Gen jets;N_{Castor};N_{GEN}", 11, -0.5, 10.5, 11, -0.5, 10.5);

	hCastorJet_energy->Sumw2();
	hCastorJet_pt->Sumw2();
	hCastorJet_em->Sumw2();
	hCastorJet_had->Sumw2();
	hCastorJet_fem->Sumw2();
	hCastorJet_fhot->Sumw2();
	hCastorJet_width->Sumw2();
	hCastorJet_depth->Sumw2();
	hCastorJet_sigmaz->Sumw2();
	hCastorJet_ntower->Sumw2();
	hCastorJet_eta->Sumw2();
	hCastorJet_phi->Sumw2();
	hCastorJet_multi->Sumw2();
	
        hCastorJet_energy_response->Sumw2();
        hCastorJet_energy_fakes->Sumw2();
        hCastorJet_energy_misses->Sumw2();
	hCastorJet_Matrix->Sumw2();

	hCastorJet_cf_energy->Sumw2();
        hCastorJet_cf_energy_gen->Sumw2();
	hCastorJet_cf_energy_response->Sumw2();
	hCastorJet_cf_energy_fakes->Sumw2();
	hCastorJet_cf_energy_misses->Sumw2();
	hCastorJet_Matrix->Sumw2();

	hJRE->Sumw2();
	hJER->Sumw2();
	hJER_per_energy->Sumw2();
	hJER_per_distance->Sumw2();
	hMatched->Sumw2();
	hUnmatched->Sumw2();
	
	TIter       next(filelist_); 
	TObjString* fn = 0;
	
	bool isMC = false;
    
	int counter_jer = 0;
	int counter_events = 0;
	
	std::cout << "start looping over files" << std::endl;
	
	// start file loop
	while((fn = (TObjString*)next()) && counter_events < totalEvents_) { 
//      while((fn = (TObjString*)next()) ) {
		
		currentfile_.Clear();
		currentTFile_->Clear();
		
		currentfile_ = fn->GetString();
		
		std::cout << "opening file " << currentfile_ << " ... " << std::endl;
		
		TStopwatch *timer = new TStopwatch();
		TThread *fileopener;
		fileopener = new TThread("fileopener",(void(*) (void *))&OpenROOTFile,(JetAnalyzer_radii*) this);
		TThread::SetCancelAsynchronous();
		TThread::SetCancelOn();
		fileopener->Run();
		
		bool fileopen = false;
		bool timeout = true;
		while (timer->RealTime() < 30) {
			if (fileopener->GetState() != TThread::kRunningState) {
				// file open thread is not running anymore, was it successful?
				// set the static TFile to the instance one
				setCurrentTFile();
				timeout = false;
				break;
			} 
			
			timer->Continue();
			gSystem->Sleep(1000);
		}
		
		// after the time, stop & delete the open thread
		TThread::Kill("fileopener");
		TThread::Delete(fileopener);
		delete fileopener;
		delete timer;
		
		if (timeout) {
			std::cout << "file open thread has timed out, go to next file in the list" << std::endl;
			continue;
		} else {
			// is this instance open?
			if (currentTFile_->IsOpen()) {
				std::cout << "the currentTFile_ is correctly opened" << std::endl;
				fileopen = true;
			}
		}
		
		if (!fileopen) {
			std::cout << "no timeout, but the file could not be opened correctly: going to the next file in the list" << std::endl;
			continue;
		}
		
		
		//////////////////////////////////////////////////
		// Get tree from the files and define all branches
		//////////////////////////////////////////////////
		
		std::cout << "We are working with Det " << det_radius_ << " and Gen " << gen_radius_ << endl;

		// get tree from file
		TTree *tree = new TTree("CastorTree","");
		currentTFile_->GetObject("castortree/CastorTree",tree);
		
		// define objects and branches
		MyEvtId *evtid = NULL;
        MyGenKin *evtkin = NULL;
		MyBeamSpot *BeamSpot = NULL;
		MyHLTrig *HLTrig = NULL;
		MyL1Trig *L1Trig = NULL;
		std::vector<MyVertex> *Vertices = NULL;
		std::vector<MyCastorRecHit> *CastorRecHits = NULL;
		std::vector<MyCastorTower> *CastorTowers = NULL;
		std::vector<MyCastorJet> *CastorJets = NULL;
		std::vector<MyJet> *PFJets = NULL;
		std::vector<MyGenPart> *genParts = NULL;
		std::vector<MyCaloTower> *caloTowers = NULL;
		std::vector<MyGenJet> *genJets = NULL;
		std::vector<MyGenJet> *chargedGenJets = NULL;
		std::vector<MyTrackJet> *trackJets = NULL;
		  // Radius-energy study.
		std::vector<MyCastorJet> *CastorJets_ak3 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak5 = NULL;
                std::vector<MyCastorJet> *CastorJets_ak7 = NULL;
                std::vector<MyGenJet> *genJets_ak3 = NULL;
                std::vector<MyGenJet> *genJets_ak5 = NULL;
                std::vector<MyGenJet> *genJets_ak7 = NULL;

		
		TBranch *b_evtid = tree->GetBranch("EvtId");
        TBranch *b_evtkin = NULL;
        if (!isData_) b_evtkin = tree->GetBranch("GenKin");
		TBranch *b_BeamSpot = tree->GetBranch("beamSpot");
		TBranch *b_HLTrig = tree->GetBranch("HLTrig");
		TBranch *b_L1Trig = tree->GetBranch("L1Trig");
		TBranch *b_vertices = tree->GetBranch("primaryVertex");
		TBranch *b_castorrechits = tree->GetBranch("castorRecHit");
		TBranch *b_castortowers = tree->GetBranch("castorTower");
		TBranch *b_castorjets = tree->GetBranch("castorJet");

		TBranch *b_castorjets_ak3 = tree->GetBranch("ak3castorJet");	b_castorjets_ak3->SetAddress(&CastorJets_ak3);
                TBranch *b_castorjets_ak5 = tree->GetBranch("ak5castorJet");    b_castorjets_ak5->SetAddress(&CastorJets_ak5);
                TBranch *b_castorjets_ak7 = tree->GetBranch("ak7castorJet");    b_castorjets_ak7->SetAddress(&CastorJets_ak7);

		if( string_det_radius.Contains("ak3") ){ b_castorjets = tree->GetBranch("ak3castorJet"); }
                if( string_det_radius.Contains("ak5") ){ b_castorjets = tree->GetBranch("ak5castorJet"); }
                if( string_det_radius.Contains("ak7") ){ b_castorjets = tree->GetBranch("ak7castorJet"); }
//		else{ b_castorjets = tree->GetBranch("castorJet"); }
		
//                TBranch *b_castorjets = tree->GetBranch("ak7castorJet");
		TBranch *b_PFJets = tree->GetBranch("pfJet");
		TBranch *b_genParts = NULL;
		if (!isData_) b_genParts = tree->GetBranch("GenPart");
		TBranch *b_caloTowers = tree->GetBranch("caloTower");
		TBranch *b_genJets = NULL;
		TBranch *b_chargedGenJets = NULL;
//		if (!isData_) b_genJets = tree->GetBranch("GenJet");
                if (!isData_ && string_gen_radius.Contains("ak3")) b_genJets = tree->GetBranch("ak3GenJet");
                if (!isData_ && string_gen_radius.Contains("ak5")) b_genJets = tree->GetBranch("ak5GenJet");
                if (!isData_ && string_gen_radius.Contains("ak7")) b_genJets = tree->GetBranch("ak7GenJet");		
		if (!isData_) b_chargedGenJets = tree->GetBranch("ChargedGenJet");
		TBranch *b_trackJets = tree->GetBranch("trackJet");

                TBranch *b_genJets_ak3 = NULL;
                TBranch *b_genJets_ak5 = NULL;
                TBranch *b_genJets_ak7 = NULL;

		if(!isData_){
		  b_genJets_ak3 = tree->GetBranch("ak3GenJet");
                  b_genJets_ak5 = tree->GetBranch("ak5GenJet");
                  b_genJets_ak7 = tree->GetBranch("ak7GenJet");
		}
		
		b_evtid->SetAddress(&evtid);
        if (!isData_) b_evtkin->SetAddress(&evtkin);
		b_BeamSpot->SetAddress(&BeamSpot);
		b_HLTrig->SetAddress(&HLTrig);
		b_L1Trig->SetAddress(&L1Trig);
		b_vertices->SetAddress(&Vertices);
		b_castorrechits->SetAddress(&CastorRecHits);
		b_castortowers->SetAddress(&CastorTowers);
		cout << "Before" << endl;
		b_castorjets->SetAddress(&CastorJets);
		cout << "After" << endl;
		b_PFJets->SetAddress(&PFJets);
		if (!isData_) b_genParts->SetAddress(&genParts);
		b_caloTowers->SetAddress(&caloTowers);
		if (!isData_) b_genJets->SetAddress(&genJets);
		if (!isData_) b_chargedGenJets->SetAddress(&chargedGenJets);
		b_trackJets->SetAddress(&trackJets);
		
		int Nevents = tree->GetEntriesFast();
		std::cout << "file opened, events in this file = " << Nevents << std::endl;
		totalevents += Nevents;
		
		// start event loop
		for (int i=0; i<Nevents && i < totalEvents_;i++) {
			counter_events++;
			if( counter_events%1000== 0){ cout << "\t" << counter_events << "\tpassed" << endl; }
            
			bool passedHadronCuts = false;
			bool passedDetectorCuts = false;
			
			double xix = 10;
			double xiy = 10;
			double xi = 10;
			double xidd = 10e10;
			double ymax = -1;
            
			
			/////////////////////////////////////////
			// Hadron level code
			/////////////////////////////////////////
			
			if (!isData_) {
				b_genParts->GetEntry(i);
				
				// calculate xi of the event
				
				// sort genParticles in y, from y_min to y_max
				std::vector<MyGenPart> myTempParticles;
				std::vector<MyGenPart> myRapiditySortedParticles;
				// copy only final stable particles with realistic Rapidity in tempvector
				for (unsigned int ipart=0;ipart<genParts->size();ipart++) {
					if ((*genParts)[ipart].status == 1) 
						myTempParticles.push_back((*genParts)[ipart]);
				}
				// do actual sorting
				while (myTempParticles.size() != 0) {
					double min_y = 100000;
					int min_y_pos = -1;
					for (unsigned int ipart = 0;ipart<myTempParticles.size();ipart++) {
						if (myTempParticles[ipart].Rapidity() < min_y) {
							min_y = myTempParticles[ipart].Rapidity();
							min_y_pos = ipart;
						}
					}
					myRapiditySortedParticles.push_back(myTempParticles[min_y_pos]);
					myTempParticles.erase(myTempParticles.begin()+min_y_pos);
				}
				
				// find deltaymax
				double deltaymax = 0;
				int deltaymax_pos = -1;
				for (unsigned int ipart=0;ipart<myRapiditySortedParticles.size()-1;ipart++) {
					double deltay = myRapiditySortedParticles[ipart+1].Rapidity() - myRapiditySortedParticles[ipart].Rapidity();
					if (deltay > deltaymax) {
						deltaymax = deltay;
						deltaymax_pos = ipart;
					}
				}
				ymax = deltaymax;
				
				// calculate Mx2 and My2
				long double XEtot = 0;
				long double XPxtot = 0;
				long double XPytot = 0;
				long double XPztot = 0;
				long double YEtot = 0;
				long double YPxtot = 0;
				long double YPytot = 0;
				long double YPztot = 0;
				
				for (int ipart=0;ipart<=deltaymax_pos;ipart++) {
					XEtot += myRapiditySortedParticles[ipart].E();
					XPxtot += myRapiditySortedParticles[ipart].Px();
					XPytot += myRapiditySortedParticles[ipart].Py();
					XPztot += myRapiditySortedParticles[ipart].Pz();
				}
				long double Mx2 = -1.;
				Mx2 = XEtot*XEtot - XPxtot*XPxtot - XPytot*XPytot - XPztot*XPztot;
				
				for (unsigned int ipart=deltaymax_pos+1;ipart<myRapiditySortedParticles.size();ipart++) {
					YEtot += myRapiditySortedParticles[ipart].E();
					YPxtot += myRapiditySortedParticles[ipart].Px();
					YPytot += myRapiditySortedParticles[ipart].Py();
					YPztot += myRapiditySortedParticles[ipart].Pz();
				}
				long double My2 = YEtot*YEtot - YPxtot*YPxtot - YPytot*YPytot - YPztot*YPztot;
				
				// calculate xix and xiy
				xix = Mx2/(7000*7000);
				xiy = My2/(7000*7000);
				
				// xi of event is max
				xi = std::max(xix,xiy);
				xidd=xix*xiy*7000*7000/(0.938*0.938);
                
                    
                // combine selection
        		// 7 TeV Xi cuts
				if (xix>0.04 || xiy>0.1 || xidd>0.5) passedHadronCuts = true;

				
				// if event passed hadron cut, execute analysis code
				if (passedHadronCuts) {
					
					
					
				} // end if passedHadronCuts
				
			} // end if not data
             
			
			
			/////////////////////////////////////////
			// Do stuff before filters
			/////////////////////////////////////////
						
			b_evtid->GetEntry(i);
			b_HLTrig->GetEntry(i);
			b_L1Trig->GetEntry(i);
			b_vertices->GetEntry(i);
			b_caloTowers->GetEntry(i);
			b_castorrechits->GetEntry(i);
			b_castortowers->GetEntry(i);
			
			// only process a certain run
			//if (evtid->Run != 135521) continue; // go to next event
			
			/////////////////////////////////////////
			// Filter the results
			/////////////////////////////////////////
			
			if (isData_) {
			
				// filter results
				// filter on phys declared bit
				bool physDeclresult = HLTrig->HLTmap["physDeclpath"];
				
				// filter on castor invalid data
				bool castorInvalidDataFilterresult = HLTrig->HLTmap["castorInvalidDataFilterpath"];
				
				// filter out scraping events
				bool noscrapingresult = HLTrig->HLTmap["noscrapingpath"];
				
				bool gooddata = physDeclresult && castorInvalidDataFilterresult && noscrapingresult;
				
				// L1 filter

				bool L1_BX = L1Trig->fTechDecisionBefore[0];
				bool L1_Veto = !L1Trig->fTechDecisionBefore[36] && !L1Trig->fTechDecisionBefore[37] && !L1Trig->fTechDecisionBefore[38] && !L1Trig->fTechDecisionBefore[39];
				
				bool L1_BSC = L1Trig->fTechDecisionBefore[40] || L1Trig->fTechDecisionBefore[41]; // default value for 900 GeV or 7000 GeV 
				
				bool HLT_BSC = true; // default value for 900 GeV or 7000 GeV (do we need an HLT here?)
				
				bool TriggerSelection = false;
                
				TriggerSelection = L1_BX && L1_Veto && L1_BSC && HLT_BSC;
				
				// ask for activity in HF
				bool HF_Activity = true;
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
				
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (gooddata && TriggerSelection && HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
				
			} else {
				
				// MC
				
				// ask for activity in HF
				bool HF_Activity = true; 
				bool HFplus = false;
				bool HFminus = false;
				for (unsigned int itow=0;itow<caloTowers->size();itow++) {
					MyCaloTower mytow = (*caloTowers)[itow];
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == 1 && mytow.Eta() > 3.23 && mytow.Eta() < 4.65) HFplus = true;
					if (mytow.hasHF && mytow.Energy() > 4. && mytow.zside == -1 && mytow.Eta() < -3.23 && mytow.Eta() > -4.65) HFminus = true;
				}
				if (!HFplus || !HFminus) HF_Activity = false; // real HF condition... 
				
				// ask for activity in CASTOR
				bool CASTOR_Activity = false;
				// if at least 1 tower is above noise threshold, set CASTOR activity to true
				if (CastorTowers->size() > 0) CASTOR_Activity = true;
                
				// get vertex info
				// do the oneGoodVertexFilter
				bool wehaveGoodVertex = false;
				for (unsigned int iVert=0;iVert<Vertices->size();iVert++) {
					MyVertex vertex = (*Vertices)[iVert];
					if (vertex.isGoodVertex) wehaveGoodVertex = true;
				}
				
				if (HF_Activity && CASTOR_Activity && wehaveGoodVertex) passedDetectorCuts = true;
                
			}
				
			/////////////////////////////////////////
			// Do stuff after filters
			/////////////////////////////////////////
			
			// get all the remaining branch entries
			b_trackJets->GetEntry(i);
                        if(!isData_){
			  b_genJets->GetEntry(i);
                          b_chargedGenJets->GetEntry(i);
                        }
				
			if (passedDetectorCuts) {
				
				
				/////////////////////////////////////////
				// Start Nvertex == 1 part of the code 
				/////////////////////////////////////////
				
				// only fill the histograms when there's 1 vertex (filter out pile-up)
				if (Vertices->size() == 1) {
					
					// get event id stuff
					if( ((i+1) % 10000) == 0) cout << " run " << evtid->Run << " isData = " << evtid->IsData << " lumiblock " << 
						evtid->LumiBlock << " event " << evtid->Evt << endl; 
					if (!evtid->IsData) isMC = true;
					
					// calculate energy flow in CASTOR
					double CASTOReflow = 0;
					
					for (unsigned int j=0;j<CastorRecHits->size();j++) {
						MyCastorRecHit rechit = (*CastorRecHits)[j];
						//if (rechit.bad) std::cout << " found bad rechit: channel " << rechit.cha << std::endl;
						if (!rechit.bad) {
							hCASTOReflow_channel[rechit.cha-1]->Fill(rechit.energy);
							if (rechit.cha <= 80) CASTOReflow += rechit.energy;
							h2CASTOReflow_grid->Fill(rechit.sec,rechit.mod,rechit.energy);
						}
					}
				
					// tower multiplicity code
					hCASTORTowerMulti->Fill(CastorTowers->size());
					
					// fill all minbias eflow histos
					hCASTOReflow->Fill(CASTOReflow);
										
					// ------------------------
					// No central jet required.		

					vector<MyCastorJet> good_castorJets;
					vector<MyGenJet> good_genJets;

					// analyse castor jets
					int NCastorJets = 0;
					vector<double> det_casjet, gen_casjet_energy;
					b_castorjets->GetEntry(i);
					for (unsigned int j=0;j<CastorJets->size();j++) {
						MyCastorJet casjet = (*CastorJets)[j];
						if (casjet.energy > jetEThreshold_det) {
							NCastorJets++;
							hCastorJet_energy->Fill(casjet.energy);
								det_casjet.push_back(casjet.energy); 
								good_castorJets.push_back( casjet );

							hCastorJet_pt->Fill(casjet.energy*sin(2*atan(exp(5.9))));
							hCastorJet_em->Fill(casjet.eem);
							hCastorJet_had->Fill(casjet.ehad);
							hCastorJet_eta->Fill(casjet.eta);
							hCastorJet_phi->Fill(casjet.phi);
							hCastorJet_fem->Fill(casjet.fem);
							hCastorJet_fhot->Fill(casjet.fhot);
							hCastorJet_width->Fill(casjet.width);
							hCastorJet_depth->Fill(casjet.depth);
							hCastorJet_sigmaz->Fill(casjet.sigmaz);
							hCastorJet_ntower->Fill(casjet.ntower);
						}
					}
					hCastorJet_multi->Fill(NCastorJets);

					// -------------------------
					// AVS - no Central jet required					
					if(!isData_){
					  for(int jet = 0; jet < genJets->size(); jet++){
					    MyGenJet currentJet = (*genJets)[jet];
					    //if( counter_events%1000 == 0) { cout << "\tGenjet\t" << jet << "\twith energy\t" << currentJet.Energy() << "\tand pT\t" << currentJet.Pt() << endl; }
					    hGenJet_energy->Fill( currentJet.Energy() );
					    if( currentJet.Eta() < (-5.2 - GenJetContained) && currentJet.Eta() > (-6.6 + GenJetContained) && currentJet.Energy() > jetEThreshold_gen){
					      gen_casjet_energy.push_back( currentJet.Energy() );
					      hCastorJet_energy_gen->Fill( currentJet.Energy() );
                                              hCastorJet_pt_gen->Fill( currentJet.Pt() );
					      good_genJets.push_back( currentJet );
					    }
					  }
					}
	
					/**********************
					 * ********************
					 * *******************/

					// Get the number of jets and the leading jet energy.
				        if( good_genJets.size() > 0) { hNjet_vs_Ejets_gen->Fill( good_genJets.size(), 	(good_genJets[ 0 ]).Energy() ); }
                                        else{ hNjet_vs_Ejets_gen->Fill(0.,0.); }
                                        if( good_castorJets.size() > 0) { hNjet_vs_Ejets_det->Fill( good_castorJets.size(), (good_castorJets[ 0 ]).energy ); }
					else{ hNjet_vs_Ejets_det->Fill(0., 0.); }

					// Match DET and GEN jets. 
					hNumber_of_match_jets->Fill( good_castorJets.size(), good_genJets.size());                             
					int matched_pairs = 0;


					while( good_castorJets.size() > 0 && good_genJets.size() > 0){ // All jets need to be matched, or at least tried to be.
										 
					/*
					// Pick a Castorjet and find the appropriate Gen Jet.
    					  MyCastorJet castorjet = good_castorJets[ 0 ];
    					  double eta_det = castorjet.eta;
    					  double phi_det = castorjet.phi;
    					  double det_energy = castorjet.energy;
					*/			
					// Pick a Gen jet and find the appropriate Castor Jet
					  MyGenJet genjet_castor = (good_genJets)[0];
					
                                          double eta_gen = genjet_castor.Eta();
                                          double phi_gen = genjet_castor.Phi();
                                          double gen_energy = genjet_castor.Energy();

					  int i_det = 0;
					  bool matched = false;
					  double lowest_distance = 10., lowest_phidiff = 10.;					
					  int match_det;
				  		  
					//
					// Dissect the jet to get hadronic/electromagnetic nature.
					//
					  int ipart = 0;
					  bool isHad = false, isEM = false, isBoth = false;
					  double eHad = 0., eEM = 0.;					 
					  TString genjettype;
 
					  for(vector<MyGenPart>::iterator it = (genjet_castor.JetPart).begin() ;it !=  (genjet_castor.JetPart).end() ;it++,++ipart) {
					    // cout << "\tipart\t" << ipart << "\t" << (*it).pdgId << endl;					   
					    if(  (*it).pdgId  == 22 || ( fabs( (*it).pdgId ) >= 11 && fabs( (*it).pdgId ) <= 18) ){ isEM = true; eEM += (*it).Energy();}
					    else if( fabs( (*it).pdgId ) > 100 ){ isHad = true; eHad += (*it).Energy();}
					  }
				
					  // In case of hybrid jets
					  if( isEM && isHad){
					    
					    if( eHad > eEM){
					      if( eHad >= 2.*eEM ){ isEM = false; }
					      else		 { isBoth = true; }
					    }
					    else{ 
                                              if( eHad < 2.*eEM ){ isHad = false; }
                                              else               { isBoth = true; }
					    }
					  } // End hybrids.
					  
					  if     ( isHad ){ genjettype = "had"; }
					  else if( isEM  ){ genjettype = "em"; }
					  else if( isBoth){ genjettype = "both";}


					  //
					  // Matching of jets.
					  //

  					  for( ; i_det < good_castorJets.size(); i_det++){
					   
                                            MyCastorJet castorjet = good_castorJets[ i_det ];
                                            double eta_det = castorjet.eta;
                                            double phi_det = castorjet.phi;

					    double phidiff = fabs(phi_det-phi_gen);	if( phidiff > PI ){ phidiff = 2.*PI - phidiff; }
					    double etadiff = fabs(eta_det-eta_gen);
					    
					    // Matching in phi.
				            if( phidiff < lowest_phidiff ){ 
					      lowest_phidiff = phidiff;
					      match_det = i_det;
					    }
					  } // Loop over Castor jets.
					
					  /* Fill in fake or match. */

					  MyCastorJet castorjet = good_castorJets[ match_det ];
					  double eta_det = castorjet.eta;
					  double phi_det = castorjet.phi;
                                          double phidiff = fabs(phi_det-phi_gen);     if( phidiff > PI ){ phidiff = 2.*PI - phidiff; }
                                          double etadiff = fabs(eta_det-eta_gen);
					  lowest_distance = sqrt( etadiff*etadiff + phidiff*phidiff );
                                          double det_energy = castorjet.energy;


					  double depth_jet = castorjet.depth;	// depth
					  double fhot_jet = castorjet.fhot;		// fhot
					  double fem_jet = castorjet.fem;
					  double sigmaz_jet = castorjet.sigmaz;
					  double width_jet = castorjet.width;

					  TString detjettype = "other";
					  // Count pions.
					  if( ! (depth_jet > -14450. && det_energy < 175.) ){ // Most likely not a pion.
					    if( ! (depth_jet > -14460. && det_energy > 175.) ){ // Most likely not a pion.
					      if( ! (fem_jet > 0.95) ){ // Most likely no pion.
					      	
						hPion_energy->Fill( det_energy );
						detjettype = "had";
					      } // Not a pion.
					    } // Not a pion.
				          } // Not a pion.
				
					  // Count electrons.
					  if( ! (fhot_jet < 0.45) ){
					    if( ! (fem_jet < 0.9) ){
					      if( ! ( sigmaz_jet > 30. && det_energy < 75.) ){
						if( ! ( sigmaz_jet > 40. && det_energy > 75. ) ){
						  if( ! ( width_jet > 11.5) ){
						    if( ! ( depth_jet < -14450. && det_energy < 125.) ){
						      if( ! ( depth_jet < -14460. && det_energy > 125.) ){						    
							hElectron_energy->Fill( det_energy );
							detjettype = "em";
						      }
						    }					
						  }
						}
					      }
					    }
					  } // Not an electron.

					/*
					  cout < endl << endl << endl 
					  << "\tEnergy\t" << det_energy << endl 
					  << "\tDepth\t" << depth_jet << endl
					  << "\tfhot\t" << fhot_jet << endl
					  << "\tfem\t" << fem_jet << endl
					  << "\tsigma_z\t" << sigmaz_jet << endl
					  << "\twidth\t" << width_jet << endl;
					*/
					  // -- looking into distance information.
					  // Leading jets.
					  if( matched_pairs == 0){
					    response.Fill( det_energy, gen_energy); 
					    
					    hCastorJet_energy_response->Fill( det_energy, gen_energy);
					 
                                            hDistance->Fill( lowest_distance );
                                            hPhiDiff -> Fill( phidiff );
                                            hEtaDiff -> Fill( etadiff );
                                            hEtaPhiDiff -> Fill( etadiff, phidiff );
                                            hEtaRDiff -> Fill( etadiff, lowest_distance );
                                            hPhiRDiff -> Fill( phidiff, lowest_distance );
					  } 
	
					  response_all.Fill( det_energy, gen_energy);
					  hDistance_all	->Fill( lowest_distance );
                                          hPhiDiff_all 	-> Fill( phidiff );
                                          hEtaDiff_all 	-> Fill( etadiff );
                                          hEtaPhiDiff_all -> Fill( etadiff, phidiff );
                                          hEtaRDiff_all -> Fill( etadiff, lowest_distance );
                                          hPhiRDiff_all -> Fill( phidiff, lowest_distance );

					  hCastorJet_energy_ratio->Fill( gen_energy/det_energy, det_energy);
					  /* Remove det and gen jet from vector. */
					  good_castorJets.erase( good_castorJets.begin() + i_det );
					  good_genJets.erase( good_genJets.begin() + 0 );
					  // matched = true;
		
					  // -- looking into JER information.
					  // 
                                          double JER = -(gen_energy - det_energy)/gen_energy;
					  double JER_eDet = (gen_energy - det_energy)/det_energy;					

                                          hJER_all		->Fill( JER );
                                          hJER_per_energy_all	->Fill( gen_energy, JER);
                                          hJER_per_distance_all	->Fill( sqrt(pow(phi_det - phi_gen, 2.) + pow(eta_det - eta_gen, 2.)), JER );
                                          hJER_per_eta_all	->Fill( eta_gen, JER);
                                          hEnergy_per_eta_all	->Fill( gen_energy, eta_gen );
					  
					  if( matched_pairs == 0){
                                            hJER		->Fill( JER );
                                            hJER_per_energy	->Fill( gen_energy, JER);
					    hJER_per_eDet	->Fill( det_energy, JER_eDet);
                                            hJER_per_distance	->Fill( sqrt(pow(phi_det - phi_gen, 2.) + pow(eta_det - eta_gen, 2.)), JER );
					    hJER_per_eta	->Fill( eta_gen, JER);
				            hEnergy_per_eta	->Fill( gen_energy, eta_gen );

					    // Look at jet types and fill in DeltaE/E..
					    if( detjettype == "had" ){
					      if( genjettype == "had") { 
						hJER_had_pi 		->Fill( JER ); 
						hCastorJet_Matrix_had_pi->Fill(det_energy, gen_energy);
						hJER_per_energy_had_pi  ->Fill( gen_energy, JER);
                                                hJER_per_eDet_had_pi  	->Fill( det_energy, JER_eDet);
					      }
                                              if( genjettype == "em") {  
						hJER_em_pi  		->Fill( JER ); 
						hCastorJet_Matrix_em_pi	->Fill(det_energy, gen_energy); 
						hJER_per_energy_em_pi  	->Fill( gen_energy, JER);
						hJER_per_eDet_em_pi	->Fill( det_energy, JER_eDet);
					      }
                                              if( genjettype == "both"){ 
						hJER_both_pi->Fill( JER ); 
					      }
					    }
                                            else if( detjettype == "em" ){
                                              if( genjettype == "had") { 
						hJER_had_e 		->Fill( JER ); 
						hCastorJet_Matrix_had_e	->Fill(det_energy, gen_energy); 
                                                hJER_per_energy_had_e  	->Fill( gen_energy, JER);
						hJER_per_eDet_had_e	->Fill( det_energy, JER_eDet);
					      }
                                              if( genjettype == "em") {  
						hJER_em_e  		->Fill( JER ); 
						hCastorJet_Matrix_em_e	->Fill(det_energy, gen_energy); 
                                                hJER_per_energy_em_e   	->Fill( gen_energy, JER);
						hJER_per_eDet_em_e	->Fill( det_energy, JER_eDet);
					      }
                                              if( genjettype == "both"){ 
						hJER_both_e->Fill( JER ); 
					      }
                                            }
					    else{
					      hJER_per_energy_none_det	->Fill( gen_energy, JER );
					      hJER_per_eDet_none_det	->Fill( det_energy, JER_eDet);
					    }

					/*
					    if( detjettype == "had" || detjettype == "em"){
					      (JER_map[genjettype + "_" + detjettype]).Fill( JER );
					*/
					    
					  }
					  matched_pairs++;

					  break; // End loop over gen jets.
					  //}// Loop gen jets.
					  if( !matched ){ // Det jet is a fake.
                                            response.Fake( det_energy );
                                            hCastorJet_energy_fakes->Fill( det_energy);
                                            good_castorJets.erase( good_castorJets.begin() +  0 );
					  }
					} // While loop.

					hMatched->Fill( matched_pairs );
					
					/* Fakes */
					for( int i_det = 0; i_det < good_castorJets.size(); i_det++){
					  MyCastorJet castorjet = (good_castorJets)[i_det];
					  double det_energy = castorjet.energy;
					  response.Fake( det_energy );
                                          hCastorJet_energy_fakes->Fill( det_energy, -50.);
					}

                                        /* Misses */
                                        for( int i_gen = 0; i_gen < good_genJets.size(); i_gen++){
					  MyGenJet genjet_castor = (good_genJets)[i_gen];
                                          double gen_energy = genjet_castor.Energy();
                                          response.Miss( gen_energy );
                                          hCastorJet_energy_misses->Fill( -50., gen_energy);
                                        }

//					cout << endl << endl;
						
					// end of event, print status
					if( ((i+1) % 10000) == 0) std::cout << i+1 <<"events done in file " << it << std::endl;
					totalevents++;
					
				} // end if statement for 1 vertex
								
			} // end if passed detector cuts
			
			// combined hadron and detector level code
			
			if (passedHadronCuts && passedDetectorCuts) {
			}
			
		} // end event loop
		
		delete tree;
		currentTFile_->Close();
		it++;
	} // end file loop
	
	std::cout << "file loop has ended" << std::endl;
    
    currentTFile_->Close();
	
	
    // check all histo's for overflow
	
	hf_.checkFlow(hCASTOReflow);
	hf_.checkFlow(hCASTORTowerMulti);
    
    // check distribution of each channel on under or overflow
	for (int icha=0;icha<224;icha++) {
		hf_.checkFlow(hCASTOReflow_channel[icha]);
	}
	
	hf_.checkFlow(hCastorJet_energy);
	hf_.checkFlow(hCastorJet_pt);
	hf_.checkFlow(hCastorJet_em);
	hf_.checkFlow(hCastorJet_had);
	hf_.checkFlow(hCastorJet_fem);
	hf_.checkFlow(hCastorJet_fhot);
	hf_.checkFlow(hCastorJet_width);
	hf_.checkFlow(hCastorJet_depth);
	hf_.checkFlow(hCastorJet_sigmaz);
	hf_.checkFlow(hCastorJet_ntower);
	hf_.checkFlow(hCastorJet_eta);
	hf_.checkFlow(hCastorJet_phi);
	hf_.checkFlow(hCastorJet_multi);
    

	// Create response matrix including misses and fakes.
	hCastorJet_energy_complete_response->Add( hCastorJet_energy_response);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_misses);
        hCastorJet_energy_complete_response->Add( hCastorJet_energy_fakes);
   
    // write all histo's to file
    
	std::cout << "total number of events = " << totalevents << " from " << it << " file(s)" << endl;
	
	
	// create output root file
	Char_t filename[200];
	const char* part = currentfile_.Data();
	std::string first(outputname_);
        float etamargin = GenJetContained;         
	//TString datestring = date;

	first += part;
	//strcat(temp,part);
//	sprintf(filename,"output_JetAnalyzer_radii_%s",first.c_str());
        if( !isData_) { sprintf(filename, date_ +"_Output_JetAnalyzer_radii_GEN_" + string_gen_radius + "_DET_" + string_det_radius + "_margin_%f_%i_%s", etamargin, counter_events, first.c_str()); }
        else{ sprintf(filename, date_ + "_Output_JetAnalyzer_radii_Data_" + string_det_radius + "_margin_%f_%i_%s", etamargin, counter_events, first.c_str()); }
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
		
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
	
	// detector level histograms
    
	// eflow histos
	hCASTORTowerMulti->Write();
	hCASTOReflow->Write();
	h2CASTOReflow_grid->Write();
	
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channel[icha]->Write();
	}

	hJRE->Divide( hMatched, hUnmatched );
	
	hCastorJet_energy->Write();
	hCastorJet_pt->Write();
	hCastorJet_em->Write();
	hCastorJet_had->Write();
	hCastorJet_fem->Write();
	hCastorJet_fhot->Write();
	hCastorJet_width->Write();
	hCastorJet_depth->Write();
	hCastorJet_sigmaz->Write();
	hCastorJet_ntower->Write();
	hCastorJet_eta->Write();
	hCastorJet_phi->Write();
	hCastorJet_multi->Write();
	hGenJet_energy->Write();

        hCastorJet_energy_gen->Write();
	hCastorJet_pt_gen->Write();

        hCastorJet_energy_response->Write();
        hCastorJet_energy_fakes->Write();
        hCastorJet_energy_misses->Write();
        hCastorJet_energy_complete_response->Write();

	hCastorJet_Matrix->Write();

	hCastorJet_energy_ratio->Write();

        hCastorJet_cf_energy->Write();
        hCastorJet_cf_energy_gen->Write();
        hCastorJet_cf_energy_response->Write();
        hCastorJet_cf_energy_fakes->Write();
        hCastorJet_cf_energy_misses->Write();
        hCastorJet_cf_Matrix->Write();

	hDistance->Write();
	hPhiDiff->Write();
        hEtaDiff->Write();
	hEtaPhiDiff->Write();
        hEtaRDiff->Write();
        hPhiRDiff->Write();

        hDistance_all->Write();
        hPhiDiff_all->Write();
        hEtaDiff_all->Write();
        hEtaPhiDiff_all->Write();
        hEtaRDiff_all->Write();
        hPhiRDiff_all->Write();

	hElectron_energy	->Write();
	hPion_energy		->Write();
	
	  hPi_e_ratio->Divide( hPion_energy, hElectron_energy );
	hPi_e_ratio		->Write();

	hNumber_of_match_jets	->Write();

	hTrackjets_2D_number	->Write();
	hTrackjets_2D_pt	->Write();

	hJER			->Write();
	hJER_per_energy		->Write();
        hJER_per_distance	->Write();
        hJER_per_eta		->Write();
        hEnergy_per_eta		->Write();

        hJER_all		->Write();
        hJER_per_energy_all	->Write();
        hJER_per_distance_all	->Write();
	hJER_per_eta_all	->Write();
	hEnergy_per_eta_all	->Write();

	/* Sort jets per type. */
	hJER_had_pi ->Write();
        hJER_had_e  ->Write();
        hJER_em_pi  ->Write();
        hJER_em_e   ->Write();
        hJER_both_pi->Write();
        hJER_both_e ->Write();

	hCastorJet_Matrix_had_pi->Write();
        hCastorJet_Matrix_had_e ->Write();
        hCastorJet_Matrix_em_pi ->Write();
        hCastorJet_Matrix_em_e	->Write();

	/* DeltaE/E as function of generator energy. */

	  hJER_per_energy_had_det->Add(hJER_per_energy_had_pi, hJER_per_energy_em_pi);
	hJER_per_energy_had_det->Write();

          hJER_per_energy_em_det->Add(hJER_per_energy_had_e, hJER_per_energy_em_e);
        hJER_per_energy_em_det->Write();
        hJER_per_energy_had_pi->Write();
        hJER_per_energy_had_e->Write();
        hJER_per_energy_em_pi->Write();
	hJER_per_energy_em_e->Write();

	hJER_per_energy_none_det->Write();

	/* DeltaE/E as function of detector energy. */

          hJER_per_eDet_had_det->Add(hJER_per_eDet_had_pi, hJER_per_eDet_em_pi);
        hJER_per_eDet_had_det->Write();

          hJER_per_eDet_em_det->Add(hJER_per_eDet_had_e, hJER_per_eDet_em_e);
        hJER_per_eDet_em_det->Write();
        hJER_per_eDet_had_pi->Write();
        hJER_per_eDet_had_e->Write();
        hJER_per_eDet_em_pi->Write();
        hJER_per_eDet_em_e->Write();

        hJER_per_energy_none_det->Write();

	hJER_per_eDet->Write();
	
	hMatched->Write();
	hUnmatched->Write();
	hJRE->Write();		

	hNjet_vs_Ejets_gen->Write();
        hNjet_vs_Ejets_det->Write();

	// Study of jet radius-energy
/*
	hNjet_3_gen->Write();
	hNjet_5_gen->Write();
	hNjet_7_gen->Write();
	
	hJet_37_gen->Write();
	hDistance_37_gen->Write();
	hPhiDiff_37_gen->Write(); 
	hEtaDiff_37_gen->Write(); 

	hJet_57_gen->Write(); 
	hDistance_57_gen->Write();
	hPhiDiff_57_gen->Write(); 
	hEtaDiff_57_gen->Write(); 

	hNjet_3_det->Write(); 
	hNjet_5_det->Write(); 
	hNjet_7_det->Write(); 

	hJet_37_det->Write(); 
	hDistance_37_det->Write();
	hPhiDiff_37_det->Write(); 
	hEtaDiff_37_det->Write(); 

	hJet_57_det->Write(); 
	hDistance_57_det->Write();
	hPhiDiff_57_det->Write();
	hEtaDiff_57_det->Write(); 
*/

	response.Write();
	response_all.Write();
//        response_cf.Write();
//	response_cf_onlyMatches.Write();

	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
 LoopOutputFile_ = filename;
        cout << totalevents << "\tevents" << endl;
	cout << counter_jer << "\tJER events" << endl;

}
	
void JetAnalyzer_radii::AfterLoopCalculations(TString file) {

	///////////////////////////////////////////////////////
	// perform calculations after event by event filling //
	///////////////////////////////////////////////////////
    
    // get the histograms from the file
    // get all the histograms
    HistoRetriever histogetter;
    std::vector<TH1D*> histovector = histogetter.getHistos(file);
    //std::vector<TH2D*> h2Dhistovector = histogetter_.get2DHistos(inputdir+file);
    
    char name [100];
    
    TH1D *hCASTOReflow_channel[224];
	for (int i=0;i<224;i++) {
		sprintf(name,"hCASTOReflow_channel_%d",i+1);
		hCASTOReflow_channel[i] = new TH1D(hf_.getHistoByName(histovector,name)); 
	}
    
    TH1D *hCASTOReflow = new TH1D(hf_.getHistoByName(histovector,"hCASTOReflow"));
	
	std::cout << "starting AfterLoop calculations" << std::endl;
	
	// get mean and error's from all channels and put it in one histo
    TH1D *hCASTOReflow_channels = new TH1D("hCASTOReflow_channels","average energy in used channels",224,1,225);
	for (int icha=0;icha<224;icha++) {
		hCASTOReflow_channels->SetBinContent(icha+1,hCASTOReflow_channel[icha]->GetMean());
		hCASTOReflow_channels->SetBinError(icha+1,hCASTOReflow_channel[icha]->GetMeanError());
	}
  	
	std::cout << "Mean CASTOR energy flow in first 5 modules = " << hCASTOReflow->GetMean() << " +/- " << hCASTOReflow->GetMeanError() << std::endl;
	
	//////////////////////////////////////////////////
	
	// create output root file
	Char_t filename[200];
	sprintf(filename,"AfterLoop_%s",file.Data());
	TFile* output = new TFile(filename,"RECREATE");
	output->cd();
	
	//////////////////////////////////////////
	// Save all your histograms in a root file
	//////////////////////////////////////////
        
  	hCASTOReflow_channels->Write();
		
	output->Close();
	std::cout << "file " << filename << " created." << std::endl;
	
}

TString JetAnalyzer_radii::getOutputFile() {
    return LoopOutputFile_;
}

TString JetAnalyzer_radii::getInputDir() {
	return inputdir_;
}

TString JetAnalyzer_radii::getCurrentFile() {
	return currentfile_;
}

void JetAnalyzer_radii::setCurrentTFile() {
	currentTFile_ = currentStaticTFile_;
}

void* JetAnalyzer_radii::OpenROOTFile(JetAnalyzer_radii* arg) {
	currentStaticTFile_ = TFile::Open(arg->getInputDir()+arg->getCurrentFile(),"READ");
	return 0;
}

int JetAnalyzer_radii::posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut) {
	
	// search for leading jets
	int posLeadingChargedGenJet = -1;

	double tempptchargedgenjet = 0;
	for (unsigned int ijet=0;ijet<JetVector.size();ijet++) {
		if (JetVector[ijet].Pt() > minptcut && fabs(JetVector[ijet].Eta()) < etacut) {
			if (JetVector[ijet].Pt() > tempptchargedgenjet) {
				tempptchargedgenjet = JetVector[ijet].Pt();
				posLeadingChargedGenJet = ijet;
			}
		}
	}
	
	return posLeadingChargedGenJet;
}

int JetAnalyzer_radii::posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut) {
	
	// search for leading jets
	int posLeadingTrackJetresult = -1;
	double temppttrack = 0;
	for (unsigned int ijet=0;ijet<JetVector.size();ijet++) {
		if (JetVector[ijet].pt_raw > minptcut && fabs(JetVector[ijet].eta_raw) < etacut && JetVector[ijet].pv) {
			if (JetVector[ijet].pt_raw > temppttrack) {
				temppttrack = JetVector[ijet].pt_raw;
				posLeadingTrackJetresult = ijet;
			}
		}
	}
	
	return posLeadingTrackJetresult;
}

