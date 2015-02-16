#ifndef JetAnalyzer_radii_strippedTree_h
#define JetAnalyzer_radii_strippedTree_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TFile.h>

#include "../src/MyJet.h"
#include "../src/MyGenJet.h"
#include "../src/MyTrackJet.h"

#include "HelperFunctions.h"

class JetAnalyzer_radii_strippedTree {
public:
	JetAnalyzer_radii_strippedTree(TString inputdir, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date, TString filename);
	virtual ~JetAnalyzer_radii_strippedTree();
    
    // Basic functions
	void Loop();
    void AfterLoopCalculations(TString file);
    
	// getter functions
    TString getOutputFile();
	TString getInputDir();
	TString getCurrentFile();
	
	//setter functions
	void setCurrentTFile();
	
	// static functions
	static void* OpenROOTFile(JetAnalyzer_radii_strippedTree* arg);
    
    // analysis specific helper functions
    int posLeadingGenJet(std::vector<MyGenJet> JetVector,double etacut,double minptcut);
	int posLeadingTrackJet(std::vector<MyTrackJet> JetVector,double etacut,double minptcut);
    
	
private:
    
    // basic variables
    TString inputdir_; 
    bool isData_; 
    const char* outputname_;
    TString gen_radius_;
    TString det_radius_;
    int totalEvents_;
    TString date_;
    TString filename_;
	
	TString currentfile_;
	TFile* currentTFile_;
	static TFile* currentStaticTFile_;
    
    TString LoopOutputFile_;
    
    HelperFunctions hf_;
    
    unsigned int Nruns;
    int nLumiBins;
    
    std::vector<int> runs;
    std::vector<double> LumiBins;
	
};

#endif
