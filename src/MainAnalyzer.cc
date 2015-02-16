#include "MainAnalyzer.h"

#include <TROOT.h>
#include "TObjString.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TList.h"
#include "TKey.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TSystem.h"
#include <TStyle.h>
#include "TPaveText.h"

#include <iostream>

using namespace std;

MainAnalyzer::MainAnalyzer() {
	
	// initialize helper classes
    setCMSStyle();
	
}

MainAnalyzer::~MainAnalyzer() { 
	
	
}


void MainAnalyzer::combineHistos(TString inputdir, TString regexpstr, double cmenergy) {
	
	TObjArray *files = reader_.getFileList(inputdir,regexpstr);
	treeoutputcombiner_.Combine(inputdir,files,cmenergy);
	
}

void MainAnalyzer::makeHadronHistos(TString inputdir, TString regexpstr, double cmenergy, const char* outputname) {
	
	TObjArray *files = reader_.getFileList(inputdir,regexpstr);
	hadronanalyzer_.Loop(inputdir,files, cmenergy, outputname);
	
}

void MainAnalyzer::makeJetHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname) {
	
	TObjArray *files = reader_.getFileList(inputdir,regexpstr);
    JetAnalyzer jetanalyzer(inputdir,files, isData, outputname);
	jetanalyzer.Loop();
    TString input = jetanalyzer.getOutputFile();
    if (input != "") jetanalyzer.AfterLoopCalculations(input);
    if (input == "") std::cout << "empty input file for AfterLoopCalculations given" << std::endl;
	
}

// AVS.
void MainAnalyzer::makeRadiusHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname) {

        TObjArray *files = reader_.getFileList(inputdir,regexpstr);
    RadiusAnalyzer radiusanalyzer(inputdir,files, isData, outputname);
        radiusanalyzer.Loop();
    TString input = radiusanalyzer.getOutputFile();
//    if (input != "") radiusanalyzer.AfterLoopCalculations(input);
//    if (input == "") std::cout << "empty input file for AfterLoopCalculations given" << std::endl;

}
/*
void MainAnalyzer::makeSysMinHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname) {

        TObjArray *files = reader_.getFileList(inputdir,regexpstr);
        SystematicsMin sysminanalyzer(inputdir,files, isData, outputname);
        sysminanalyzer.Loop();
        TString input = sysminanalyzer.getOutputFile();

}
void MainAnalyzer::makeSysMaxHistos(TString inputdir, TString regexpstr, bool isData, const char* outputname) {

        TObjArray *files = reader_.getFileList(inputdir,regexpstr);
        SystematicsMax sysmaxanalyzer(inputdir,files, isData, outputname);
        sysmaxanalyzer.Loop();
        TString input = sysmaxanalyzer.getOutputFile();

}
*/
void MainAnalyzer::makeJetHistos_radii(TString inputdir, TString regexpstr, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date) {

        TObjArray *files = reader_.getFileList(inputdir,regexpstr);
	JetAnalyzer_radii jetanalyzer_radii(inputdir,files, isData, outputname, gen_radius, det_radius, totalEvents, date);
        jetanalyzer_radii.Loop();
    	TString input = jetanalyzer_radii.getOutputFile();
//    if (input != "") jetanalyzer_radii.AfterLoopCalculations(input);
    if (input == "") std::cout << "empty input file for AfterLoopCalculations given" << std::endl;
}
void MainAnalyzer::makeJetHistos_stripTheTree(TString inputdir, TString regexpstr, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date) {

        TObjArray *files = reader_.getFileList(inputdir,regexpstr);
        JetAnalyzer_stripTheTree jetanalyzer_stripTheTree(inputdir,files, isData, outputname, gen_radius, det_radius, totalEvents, date);
        jetanalyzer_stripTheTree.Loop();
        TString input = jetanalyzer_stripTheTree.getOutputFile();
//    if (input != "") jetanalyzer_stripTheTree.AfterLoopCalculations(input);
        if (input == "") std::cout << "empty input file for AfterLoopCalculations given" << std::endl;
}

void MainAnalyzer::makeJetHistos_radii_strippedTree(TString inputdir, bool isData, const char* outputname, TString gen_radius, TString det_radius, int totalEvents, TString date, TString filename) {

        JetAnalyzer_radii_strippedTree jetanalyzer_radii_strippedTree(inputdir, isData, outputname, gen_radius, det_radius, totalEvents, date, filename);
        jetanalyzer_radii_strippedTree.Loop();
        TString input = jetanalyzer_radii_strippedTree.getOutputFile();
//    if (input != "") jetanalyzer_radii_strippedTree.AfterLoopCalculations(input);
        if (input == "") std::cout << "empty input file for AfterLoopCalculations given" << std::endl;
}



// End AVS.

void MainAnalyzer::makeJetAfterLoopHistos(TString inputfile, bool isData, const char* outputname) {
	
    JetAnalyzer jetanalyzer("",NULL, isData, outputname);
	jetanalyzer.AfterLoopCalculations(inputfile);
	
}

void MainAnalyzer::plotSingleHistos(TString outputfile, TString selectname) {
	
	// reset canvas vector
	canvasvector_.clear();
	
	std::vector<TH1D*> histovector = histogetter_.getHistos(outputfile);
	std::vector<TH1D*> selectedhistos;
	selectedhistos.clear();
	
	for (unsigned int i=0;i<histovector.size();i++) {
		TString name = histovector[i]->GetName();
		if (name.Contains(selectname)) selectedhistos.push_back(histovector[i]);
	}
	
	TCanvas *c[selectedhistos.size()];
	
	// set plot styles
	gStyle->SetOptStat(111111);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetStatColor(kWhite);
	
	for (unsigned int i=0;i<selectedhistos.size();i++) {
		Char_t cname[50];
		sprintf(cname,"c_%s",selectedhistos[i]->GetName());
		
		c[i] = new TCanvas(cname,selectedhistos[i]->GetTitle());
		c[i]->Divide(1,1);
		c[i]->cd(1);
		selectedhistos[i]->Draw("elp");
		canvasvector_.push_back(c[i]);
	}
}


void MainAnalyzer::plotScaleHisto(TString inputdir,TString regexpstr,TString selectname) {
	
	//-- reset canvas vector
	canvasvector_.clear();
	
	TObjArray* filelist = reader_.getFileList(inputdir,regexpstr);	
	TIter next(filelist);
	TObjString* itfile = 0;
	
	std::vector<TString> filenamelist;
	std::vector<std::vector<TH1D*> > histolist;
	std::vector<std::vector<double> > entrylist;
	
	//-- loop over files
	while((itfile = (TObjString*)next())) {
		
		std::vector<TH1D*> histo = histogetter_.getHistos(inputdir+itfile->GetString());
		
		std::vector<TH1D*> selhisto;
		selhisto.clear();
		
		std::vector<double> nentry;
		nentry.clear();
		
		//-- loop to get selected histos
		for (unsigned int i = 0; i < histo.size(); i++) {
			
			TString name = histo[i]->GetName();
			
			if (name.CompareTo(selectname) != 0) continue;
			
			selhisto.push_back(histo[i]);
			nentry.push_back(histo[i]->GetEntries());
		}
		
		histolist.push_back(selhisto);
		entrylist.push_back(nentry);
		filenamelist.push_back(itfile->GetString());
	}
	
	
	//-- compute the weight
	
	int ifile_data = -1;
	std::vector<std::vector<double> > weightlist;
	
	for (unsigned int ifile = 0; ifile <histolist.size(); ifile++) {
		if(filenamelist[ifile].Contains("data")) ifile_data = ifile;
	}
	
	for (unsigned int ifile = 0; ifile <histolist.size(); ifile++) {
		
		std::vector<double> weight;
		weight.clear();
		
		for (unsigned int ihisto = 0; ihisto < histolist[0].size(); ihisto++) {
			weight.push_back(entrylist[ifile_data][ihisto] / entrylist[ifile][ihisto]);
		}
		
		weightlist.push_back(weight);
	}
	
	//-- print some information
	
	for (unsigned int ifile = 0; ifile <histolist.size(); ifile++) {
		
		cout<<endl<<"file: "<<filenamelist[ifile].Data()<<endl;
		
		for (unsigned int ihisto = 0; ihisto < histolist[0].size(); ihisto++) {
			
			cout<<"histo : "<<histolist[ifile][ihisto]->GetTitle()<<" entries: "<<entrylist[ifile][ihisto]<<" weight: "<<weightlist[ifile][ihisto]<<endl;
		}
	}
	
	//-- scale histograms
	
	for (unsigned int ifile = 0; ifile <histolist.size(); ifile++) {
		
		for (unsigned int ihisto = 0; ihisto < histolist[0].size(); ihisto++) {
			
			//histolist[ifile][ihisto]->Scale(weightlist[ifile][ihisto]);
			if (histolist[ifile][ihisto]->Integral() != 0) histolist[ifile][ihisto]->Scale(1./histolist[ifile][ihisto]->Integral()); // just scale them to their integral
		}
	}
	
	//-- plot everything
	
	TCanvas *c[histolist[0].size()];
	setPlotStyle();
	
	//-- loop over histos/Users/hans/tempEdits/CastorTree/Analyzer/CASTORTowerMulti_Cut1p46.pdf
	for (unsigned int ihisto = 0; ihisto < histolist[0].size(); ihisto++) { 
		
		double max = 0;
		double min = 10000000;
		
		//-- loop over files to calculate min and max in y
		for (unsigned int ifile = 0; ifile <histolist.size(); ifile++) { 
			double tempmax = histolist[ifile][ihisto]->GetBinContent(histolist[ifile][ihisto]->GetMaximumBin());
			double tempmin = histolist[ifile][ihisto]->GetBinContent(histolist[ifile][ihisto]->GetMinimumBin());
			if (tempmax >= max) max = tempmax;
			if (tempmin <= min) min = tempmin;
		}
		
		if (max-min > 1000 && min == 0) min = 0.001;
		if (max-min < 1000) min = 0;
		//min = 0.0001;
		//max = 2;
		
		// std::cout << "min = " << min << " max = " << max << std::endl;
		
		Char_t cname[50];
		sprintf(cname,"%s",histolist[0][ihisto]->GetName());
		
		c[ihisto] = new TCanvas(cname,histolist[0][ihisto]->GetTitle());
		c[ihisto]->Divide(1,2);
		
		TPad *p1 = (TPad*)c[ihisto]->cd(1);
		p1->SetPad(0.0,0.2,1,1);
		if (max-min > 1000) p1->SetLogy();
		//p1->SetLogy();
		
		if (ifile_data == -1) {
			// no data, just plot all MC
			for (unsigned int ifile = 0; ifile < histolist.size(); ifile++) { 
				histolist[ifile][ihisto]->GetYaxis()->SetRangeUser(min,max*1.2);
				histolist[ifile][ihisto]->SetLineColor(ifile+1);
				histolist[ifile][ihisto]->SetLineWidth(2);
				if (ifile==0) histolist[ifile][ihisto]->Draw("");
				if (ifile!=0) histolist[ifile][ihisto]->Draw("same");
			}
		} else {
			// there's data, plots it first
			histolist[ifile_data][ihisto]->GetYaxis()->SetRangeUser(min,max*1.2);
			histolist[ifile_data][ihisto]->SetLineColor(1);
			histolist[ifile_data][ihisto]->SetLineWidth(3);
			histolist[ifile_data][ihisto]->Draw("");
			// then MC
			int color = 1;
			for (unsigned int ifile = 0; ifile < histolist.size(); ifile++) { 
				if ((int)ifile != ifile_data) {
					histolist[ifile][ihisto]->GetYaxis()->SetRangeUser(min,max*1.2);
					histolist[ifile][ihisto]->SetLineColor(color+1);
					histolist[ifile][ihisto]->SetMarkerStyle(0);
					histolist[ifile][ihisto]->SetLineWidth(2);
					histolist[ifile][ihisto]->Draw("same");
					color++;
				}
				
			}
		}
		
		//-- get info for legend
		
		std::vector<TString> initial_title;
		
		//-- loop over files to retrieve initial title
		for (unsigned int ifile = 0; ifile < histolist.size(); ifile++)
			initial_title.push_back(histolist[ifile][ihisto]->GetTitle());
		
		//-- loop over files to set file title to the histo
		for (unsigned int ifile = 0; ifile < histolist.size(); ifile++) 
			histolist[ifile][ihisto]->SetTitle(filenamelist[ifile]);
		
		TPad *p2 = (TPad*)c[ihisto]->cd(2);
		p2->SetPad(0.0,0.0,1,0.2);	
		TLegend *legend = new TLegend(0.1,0.1,0.90,0.90);
		legend->SetMargin(0.1);
		legend->SetFillColor(kWhite);
		
		//-- loop over files to plot the legends
		// first data
		if (ifile_data != -1) legend->AddEntry(histolist[ifile_data][ihisto],histolist[ifile_data][ihisto]->GetTitle(),"lpf");
		// then MC
		for (unsigned int ifile = 0; ifile < histolist.size(); ifile++) {
			if ((int)ifile != ifile_data)legend->AddEntry(histolist[ifile][ihisto],histolist[ifile][ihisto]->GetTitle(),"lpf");
		}
		
		legend->Draw();
		
		//-- loop over files to set back initial title
		for (unsigned int ifile = 0; ifile < histolist.size(); ifile++)
			histolist[ifile][ihisto]->SetTitle(initial_title[ifile]);
		
		canvasvector_.push_back(c[ihisto]);
		if(ihisto == histolist[0].size() - 1) c[ihisto]->WaitPrimitive();
	}
	
}

void MainAnalyzer::plotHistos(TString inputdir,TString regexpstr, TString selectname) {
	
	// reset canvas vector
	canvasvector_.clear();
	
	TObjArray *files = reader_.getFileList(inputdir,regexpstr);
	
	std::vector<std::vector<TH1D*> > allhistovector;
	
	TIter       next(files); 
	TObjString* fn = 0;
	std::vector<TString> files_string;
	
	while((fn = (TObjString*)next())) {
		std::vector<TH1D*> histovector = histogetter_.getHistos(inputdir+fn->GetString());
		std::vector<TH1D*> selectedhistos;
		selectedhistos.clear();
		
		for (unsigned int i=0;i<histovector.size();i++) {
			TString name = histovector[i]->GetName();
			if (name.Contains(selectname)) selectedhistos.push_back(histovector[i]);
		}
		
		allhistovector.push_back(selectedhistos);
		files_string.push_back(fn->GetString());
	}
	
	
	TCanvas *c[allhistovector[0].size()];
	setPlotStyle(); // set plot style
	
	// plot everything
	for (unsigned int i=0;i<allhistovector[0].size();i++) { // loop over number of histos
		
		double max = 0;
		double min = 10000000;
		for (unsigned int j=0;j<allhistovector.size();j++) { // loop over number of files
			// calculate min and max ranges
			double tempmax = allhistovector[j][i]->GetBinContent(allhistovector[j][i]->GetMaximumBin());
			double tempmin = allhistovector[j][i]->GetBinContent(allhistovector[j][i]->GetMinimumBin());
			if (tempmax >= max) max = tempmax;
			if (tempmin <= min) min = tempmin;
		}
		
		if (max-min > 1000 && min == 0) min = 0.001;
		if (max-min < 1000) min = 0;
		
		//std::cout << "min = " << min << " max = " << max << std::endl;
		
		Char_t cname[50];
		sprintf(cname,"c_%s",allhistovector[0][i]->GetName());
		c[i] = new TCanvas(cname,allhistovector[0][i]->GetTitle());
		c[i]->Divide(1,2);
		TPad *p1 = (TPad*)c[i]->cd(1);
		p1->SetPad(0.0,0.2,1,1);
		if (max-min > 1000) p1->SetLogy();
		
		for (unsigned int j=0;j<allhistovector.size();j++) { // loop over number of files
			allhistovector[j][i]->GetYaxis()->SetRangeUser(min,max*1.2);
			allhistovector[j][i]->SetLineColor(j+1);
            allhistovector[j][i]->SetStats(0);
			if (j==0) allhistovector[j][i]->SetLineWidth(3);
			if (j==0) allhistovector[j][i]->Draw("elp");
			if (j!=0) allhistovector[j][i]->Draw("elpsame");
		}
		
		// get info for legend
		TString initial_title = allhistovector[0][i]->GetTitle();
		for (unsigned int j=0;j<allhistovector.size();j++) { // loop over number of files
			allhistovector[j][i]->SetTitle(files_string[j]);
		}
		
		TPad *p2 = (TPad*)c[i]->cd(2);
		p2->SetPad(0.0,0.0,1,0.2);
		TLegend *legend = new TLegend(0.1,0.1,0.90,0.90);
		legend->SetMargin(0.1);
		legend->SetFillColor(kWhite);
		
		for (unsigned int j=0;j<allhistovector.size();j++) { // loop over number of files
			legend->AddEntry(allhistovector[j][i],allhistovector[j][i]->GetTitle(),"lpf");
		}
		legend->Draw();
		
		// set back initial title to first histogram
		allhistovector[0][i]->SetTitle(initial_title);
		
		canvasvector_.push_back(c[i]);
		//if(i == allhistovector[0].size() - 1) c[i]->WaitPrimitive();
	}
}

void MainAnalyzer::compareHistogramContents(TString file1, TString file2) {
	
	HistoRetriever histogetter;
    std::vector<TH1D*> file1_histos = histogetter.getHistos(file1);
	std::vector<TH1D*> file2_histos = histogetter.getHistos(file2);
	
	std::cout << "Number of histograms in file1 = " << file1_histos.size() << " number of histograms in file2 = " << file2_histos.size() << std::endl;
	int Nbadhistos = 0;
	std::vector<TString> badhistos;
	if (file1_histos.size() >= file2_histos.size()) {
		std::cout << "file1 has more or equal number of histograms, let's loop on file2..." << std::endl;
		for (unsigned int i=0;i<file2_histos.size();i++) {
			for (unsigned int j=0;j<file1_histos.size();j++) {
				if (strcmp((*file2_histos[i]).GetName(),(*file1_histos[j]).GetName()) == 0) {
					std::cout << "	found histogram " << (*file2_histos[i]).GetName() << " in both files" << std::endl;
					bool contentsok = true;
					if ((*file2_histos[i]).GetNbinsX() != (*file1_histos[j]).GetNbinsX()) {
						std::cout << "		the histograms have a different number of bins" << std::endl;
						contentsok = false;
					}
					for (int ibin=0;ibin<(*file2_histos[i]).GetNbinsX();ibin++) {
						if ((*file2_histos[i]).GetBinContent(ibin+1) != (*file1_histos[j]).GetBinContent(ibin+1)) {
							contentsok = false;
							std::cout << "		found discrepancy in bin " << ibin+1 << std::endl;
						}
					}
					if (contentsok) std::cout << "	the histograms have exactly the same bin content" << std::endl;
					if (!contentsok) { Nbadhistos++; badhistos.push_back((*file2_histos[i]).GetName());}
				}
			}
		}
	}
	std::cout << " " << std::endl;
	std::cout << "Summary:" << std::endl;
	std::cout << "the comparison found " << Nbadhistos << " histograms that are not the same!" << std::endl;
	std::cout << "the following have discrepancies: " << std::endl;
	for (unsigned int i=0;i<badhistos.size();i++) {
		std::cout << badhistos[i] << ", ";
	}
	std::cout << "" << std::endl;
	
	if (file1_histos.size() < file2_histos.size()) std::cout << "file2 has more histograms, please interchange the arguments" << std::endl;
	
	
	
}

void MainAnalyzer::saveAllCanvasPDF(TString inputdir,TString name) {
	
	cout<<endl<<"begin to save canvas"<<endl;
	
	TString file_pdf = TString(inputdir) + TString("plot_") + TString(name) + TString(".pdf");
	
	canvasvector_[0]->Print(TString(TString(file_pdf)+TString("[")).Data());
	
	for (unsigned int i=0;i<canvasvector_.size();i++) {
		TCanvas *c = canvasvector_[i];
		canvasvector_[i]->Print(file_pdf.Data());
		TString cname;
		cname.Append(inputdir);
		cname.Append(c->GetName());
		cname.Append("_");
		cname.Append(name);
		cname.Append(".png");
		c->SaveAs(cname);
	}
	
	canvasvector_[0]->Print(TString(TString(file_pdf)+TString("]")).Data());
	cout<<"canvas saved !"<<endl;
}

void MainAnalyzer::saveAllCanvas(TString inputdir,TString name) {
	
	cout<<"begin to save canvas"<<endl;
	
	TString file_pdf = TString(inputdir) + TString("plot_") + TString(name) + TString(".pdf");
	
	canvasvector_[0]->Print(TString(TString(file_pdf)+TString("[")).Data());
	
	for (unsigned int i=0;i<canvasvector_.size();i++) {
		TCanvas *c = canvasvector_[i];
		canvasvector_[i]->Print(file_pdf.Data());
		TString cname;
		cname.Append(inputdir);
		cname.Append(c->GetName());
		cname.Append("_");
		cname.Append(name);
		cname.Append(".png");
		c->SaveAs(cname);
	}
	
	canvasvector_[0]->Print(TString(TString(file_pdf)+TString("]")).Data());
	cout<<"canvas saved !"<<endl;
}

void MainAnalyzer::setPlotStyle() {
	
	//-- set plot styles
	gStyle->SetOptStat(111111);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(kWhite);
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetStatColor(kWhite);
    
}


void MainAnalyzer::setCMSStyle(){
	
	std::cout << "CMS Style Loaded" << std::endl;
	
	TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
	
	// For the canvas:
	tdrStyle->SetCanvasBorderMode(0);
	tdrStyle->SetCanvasColor(kWhite);
	// tdrStyle->SetCanvasDefH(200); //Height of canvas
	//tdrStyle->SetCanvasDefW(200); //Width of canvas
	tdrStyle->SetCanvasDefX(0);   //POsition on screen
	tdrStyle->SetCanvasDefY(0);
	
	// For the Pad:
	tdrStyle->SetPadBorderMode(0);
	tdrStyle->SetPadBorderSize(0);
	tdrStyle->SetPadColor(kWhite);
	//   tdrStyle->SetPadGridX(false);
	//   tdrStyle->SetPadGridY(false);
	//   tdrStyle->SetGridColor(0);
	//   tdrStyle->SetGridStyle(3);
	//   tdrStyle->SetGridWidth(1);
	
	// For the frame:
	tdrStyle->SetFrameBorderMode(0);
	tdrStyle->SetFrameBorderSize(0);
	tdrStyle->SetFrameFillColor(0);
	tdrStyle->SetFrameFillStyle(0);
	tdrStyle->SetFrameLineColor(1);
	tdrStyle->SetFrameLineStyle(1);
	tdrStyle->SetFrameLineWidth(1);
	
	// For the histo:
	// tdrStyle->SetHistFillColor(1);
	// tdrStyle->SetHistFillStyle(0);
	tdrStyle->SetHistLineColor(1);
	tdrStyle->SetHistLineStyle(0);
	tdrStyle->SetHistLineWidth(1);
	// tdrStyle->SetLegoInnerR(Float_t rad =3D 0.5);
	// tdrStyle->SetNumberContours(Int_t number =3D 20);
	
	tdrStyle->SetEndErrorSize(2);
	//tdrStyle->SetErrorMarker(20);
	//tdrStyle->SetErrorX(0.);
	
	tdrStyle->SetMarkerStyle(20);
	
	//For the fit/function:
	tdrStyle->SetOptFit(0);
	tdrStyle->SetFitFormat("5.4g");
	tdrStyle->SetFuncColor(2);
	tdrStyle->SetFuncStyle(1);
	tdrStyle->SetFuncWidth(1);
	
	//For the date:
	tdrStyle->SetOptDate(0);
	// tdrStyle->SetDateX(Float_t x =3D 0.01);
	// tdrStyle->SetDateY(Float_t y =3D 0.01);
	
	// For the statistics box:
	//  tdrStyle->SetOptFile(0);
	tdrStyle->SetOptStat(0); // To display the mean and RMS:   = SetOptStat("mr");
	tdrStyle->SetStatColor(kWhite);
	tdrStyle->SetStatFont(42);
	tdrStyle->SetStatFontSize(0.025);
	tdrStyle->SetStatTextColor(kBlack);
	tdrStyle->SetStatFormat("6.4g");
	tdrStyle->SetStatBorderSize(0);
	tdrStyle->SetStatH(0.1);
	tdrStyle->SetStatW(0.15);
	// tdrStyle->SetStatStyle(Style_t style =3D 1001);
	// tdrStyle->SetStatX(Float_t x =3D 0);
	// tdrStyle->SetStatY(Float_t y =3D 0);
	
	// Margins:
	tdrStyle->SetPadTopMargin(0.05);
	tdrStyle->SetPadBottomMargin(0.13);
	tdrStyle->SetPadLeftMargin(0.13);
	tdrStyle->SetPadRightMargin(0.05);
	
	// For the Global title:
	
	tdrStyle->SetOptTitle(0);
	tdrStyle->SetTitleFont(42);
	tdrStyle->SetTitleColor(1);
	tdrStyle->SetTitleTextColor(kWhite);
	tdrStyle->SetTitleFillColor(kWhite);
	tdrStyle->SetTitleFontSize(0.05);
	// tdrStyle->SetTitleH(0); // Set the height of the title box
	// tdrStyle->SetTitleW(0); // Set the width of the title box
	// tdrStyle->SetTitleX(0); // Set the position of the title box
	// tdrStyle->SetTitleY(0.985); // Set the position of the title box
	// tdrStyle->SetTitleStyle(Style_t style =3D 1001);
	// tdrStyle->SetTitleBorderSize(2);
	
	// For the axis titles:
	
	tdrStyle->SetTitleColor(1, "XYZ");
	tdrStyle->SetTitleFont(42, "XYZ");
	tdrStyle->SetTitleSize(0.06, "XYZ");
	// tdrStyle->SetTitleXSize(Float_t size =3D 0.02); // Another way to = set the size?
	// tdrStyle->SetTitleYSize(Float_t size =3D 0.02);
	tdrStyle->SetTitleXOffset(0.9);
	tdrStyle->SetTitleYOffset(1.05);
	// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the = Offset
	
	// For the axis labels:
	
	tdrStyle->SetLabelColor(1, "XYZ");
	tdrStyle->SetLabelFont(42, "XYZ");
	tdrStyle->SetLabelOffset(0.007, "XYZ");
	tdrStyle->SetLabelSize(0.05, "XYZ");
	
	// For the axis:
	
	tdrStyle->SetAxisColor(1, "XYZ");
	tdrStyle->SetStripDecimals(kTRUE);
	tdrStyle->SetTickLength(0.03, "XYZ");
	tdrStyle->SetNdivisions(510, "XYZ");
	tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side = of the frame
	tdrStyle->SetPadTickY(1);
	
	// Change for log plots:
	//   tdrStyle->SetOptLogx(0);
	//   tdrStyle->SetOptLogy(0);
	//   tdrStyle->SetOptLogz(0);
	
	// Postscript options:
	//  tdrStyle->SetPaperSize(7.5,7.5);
	// tdrStyle->SetLineScalePS(Float_t scale =3D 3);
	// tdrStyle->SetLineStyleString(Int_t i, const char* text);
	// tdrStyle->SetHeaderPS(const char* header);
	// tdrStyle->SetTitlePS(const char* pstitle);
	
	// tdrStyle->SetBarOffset(Float_t baroff =3D 0.5);
	// tdrStyle->SetBarWidth(Float_t barwidth =3D 0.5);
	// tdrStyle->SetPaintTextFormat(const char* format =3D "g");
	// tdrStyle->SetPalette(Int_t ncolors =3D 0, Int_t* colors =3D 0);
	// tdrStyle->SetTimeOffset(Double_t toffset);
	// tdrStyle->SetHistMinimumZero(kTRUE);
	
	// for the legend
	tdrStyle->SetLegendBorderSize(0);
	//tdrStyle->SetLegendFillColor(0);
	
	tdrStyle->cd();
	// gSystem->Load("libRooStats");
	// using namespace RooFit;
	
}

void MainAnalyzer::drawCMSLabels(double cmenergy, double lumi) {
	
	TPaveText *plotlabel = new TPaveText(0.23,0.87,0.43,0.92,"NDC");
	plotlabel->SetTextColor(kBlack);
	plotlabel->SetFillColor(kWhite);
	plotlabel->SetBorderSize(0);
	plotlabel->SetTextAlign(12);
	plotlabel->SetTextSize(0.03);
	plotlabel->AddText("CMS Preliminary");
	plotlabel->Draw();
	
	TPaveText *plotlabel2 = new TPaveText(0.23,0.82,0.43,0.87,"NDC");
	plotlabel2->SetTextColor(kBlack);
	plotlabel2->SetFillColor(kWhite);
	plotlabel2->SetBorderSize(0);
	plotlabel2->SetTextAlign(12);
	plotlabel2->SetTextSize(0.03);
	if (cmenergy == 900) {
		plotlabel2->AddText("#sqrt{s} = 900 GeV");
	} else if (cmenergy == 2760) {
		plotlabel2->AddText("#sqrt{s} = 2.76 TeV");
	} else if (cmenergy == 7000) {
		plotlabel2->AddText("#sqrt{s} = 7 TeV");
	}
	plotlabel2->Draw();
	
	TPaveText *plotlabel3 = new TPaveText(0.23,0.75,0.43,0.80,"NDC");
	plotlabel3->SetTextColor(kBlack);
	plotlabel3->SetFillColor(kWhite);
	plotlabel3->SetBorderSize(0);
	plotlabel3->SetTextAlign(12);
	plotlabel3->SetTextSize(0.03);
	char temp[100];
	sprintf(temp, "%.4f", lumi);
	plotlabel3->AddText((string("#int#font[12]{L}dt = ") + temp + string(" pb^{ -1}")).c_str());
	plotlabel3->Draw();
	
}

