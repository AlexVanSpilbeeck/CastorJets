#ifndef MainAnalyzer_h
#define MainAnalyzer_h

#include <TString.h>
#include <TRegexp.h>
#include <TObjArray.h>
#include <TCanvas.h>

#include "FileReader.h"
#include "HistoRetriever.h"
#include "JetAnalyzer.h"
#include "RadiusAnalyzer.h"
//#include "SystematicsMin.h"
//#include "SystematicsMax.h"
#include "JetAnalyzer_radii.h"
#include "JetAnalyzer_stripTheTree.h"
#include "JetAnalyzer_radii_strippedTree.h"


#include "TreeOutputCombiner.h"
#include "HadronAnalyzer.h"

class MainAnalyzer {
public:
	MainAnalyzer();
	virtual ~MainAnalyzer();
    
    void combineHistos(TString inputdir, TString regexpstr, double cmenergy);
    void makeHadronHistos(TString inputdir, TString regexpstr, double cmenergy, const char* outputname);
	
	void makeJetHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname);
      void makeJetAfterLoopHistos(TString inputfile, bool isData, const char* outputname);

	void makeRadiusHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname);	
        void makeSysMinHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname);
        void makeSysMaxHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname);
        void makeJetHistos_radii(TString inputdir, TString regexpstr, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date);
        void makeJetHistos_stripTheTree(TString inputdir, TString regexpstr, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date);
        void makeJetHistos_radii_strippedTree(TString inputdir, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date, TString filename);

	void plotSingleHistos(TString outputfile, TString selectname);
	
	void plotHistos(TString inputdir, TString regexpstr, TString selectname);
	void plotScaleHisto(TString inputdir,TString regexpstr,TString selectname);
	
	void compareHistogramContents(TString file1, TString file2);
    
	void saveAllCanvas(TString inputdir, TString name);
	void saveAllCanvasPDF(TString inputdir,TString name);
	
	void setPlotStyle();
	void setCMSStyle(); 
	void drawCMSLabels(double cmenergy, double lumi);
	
private:
	FileReader reader_;
	HistoRetriever histogetter_;

    TreeOutputCombiner treeoutputcombiner_;
	HadronAnalyzer hadronanalyzer_;
	
	std::vector<TCanvas*> canvasvector_;
    
};

#endif
