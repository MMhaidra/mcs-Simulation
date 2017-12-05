//******************************************************************************************************************************************************
// Author:	Daniel Parry *
// Version:1.0
// Date:	28-SEP-2015 *
// Description: This class sets the framework for making a histogram of scattering events
//******************************************************************************************************************************************************

#include "ScatteringHistogram.h"
//#include "ScatteringEvent.h"
#include "TROOT.h"
#include "TH1.h"


// Constructor for the Scattering Histogram which sets the histogram parameters, N the number of scattering events, the number of bins and the
// range in the x-axis. The muon properties of nominal momentum, variance and angular resolution of the source, and the length of the materials.
ScatteringHistogram::ScatteringHistogram(int NIn, int binsIn, double rangeIn, double MaterialLengthIn, double NominalMomentumIn, double MomentumVarianceIn, double AngularResolutionIn) {
		
		N = NIn;
		bins = binsIn;
		range = rangeIn;
		MaterialLength = MaterialLengthIn;
		NominalMomentum = NominalMomentumIn;
		MomentumVariance = MomentumVarianceIn;
		AngularResolution = AngularResolutionIn;
}


// Constructor for the Scattering Histogram which sets the histogram parameters, N the number of scattering events, the number of bins and the
// range in the x-axis. The muon properties of nominal momentum, variance and angular resolution of the source, and the length of the materials.
ScatteringHistogram::ScatteringHistogram(int NIn, int binsIn, double rangeIn, double MaterialLengthIn, double NominalMomentumIn, double MomentumVarianceIn) {

	N = NIn;
	bins = binsIn;
	range = rangeIn;
	MaterialLength = MaterialLengthIn;
	NominalMomentum = NominalMomentumIn;
	MomentumVariance = MomentumVarianceIn;
}

// Deconstructor for ScatteringHistogram.
ScatteringHistogram::~ScatteringHistogram() {}

// This method returns a histogram when the Radiation Length, X0 is specified. It also requires a String input to identifiy the material
// the histogram refers to and the line color.
TH1D *ScatteringHistogram::calculateHistogram(double X0,char* IDstring, int colour) {
	char titlestr[128];
	strcpy(titlestr, IDstring);
	strcat(titlestr, " vs Events");
	TH1D *Histogram = new TH1D(IDstring, titlestr, bins, 0, range);
	for (int i = 0; i < N; ++i) {
		ScatteringEvent *Event = new ScatteringEvent(X0, MaterialLength,NominalMomentum,MomentumVariance,AngularResolution);
		double ScatteringStrength = Event->getScatteringStrength();
		Histogram->Fill(ScatteringStrength);
	}
	Histogram->SetLineColor(colour);
	return Histogram;
}

TH1D *ScatteringHistogram::calculateHistogram(double X0, char* IDstring, int colour, double AngularResolutionIn) {
	char titlestr[128];
	strcpy(titlestr, IDstring);
	strcat(titlestr, " vs Events");
	TH1D *Histogram = new TH1D(IDstring, titlestr, bins, 0, range);
	for (int i = 0; i < N; ++i) {
		ScatteringEvent *Event = new ScatteringEvent(X0, MaterialLength, NominalMomentum, MomentumVariance, AngularResolutionIn);
		double ScatteringStrength = Event->getScatteringStrength();
		Histogram->Fill(ScatteringStrength);
	}
	Histogram->SetLineColor(colour);
	return Histogram;
}

double ScatteringHistogram::calculateProbability(double X01, double X02, char* IDString1, char* IDString2, int colour1, int colour2) {
	TH1D *hist1 = calculateHistogram(X01, IDString1, colour1);
	TH1D *hist2 = calculateHistogram(X02, IDString2, colour2);

	hist1->Sumw2();
	hist1->Divide(hist2);
	hist1->Draw();
	TF1 *f = new TF1("p", "pol0", 0.0,range);
	hist1->Fit(f, "QN");
	double probability = f->GetProb();
	hist1->Delete();
	hist2->Delete();
	return probability;

}

double ScatteringHistogram::calculateProbability(double X01, double X02, char* IDString1, char* IDString2, int colour1, int colour2, double AngularResolutionIn) {
	TH1D *hist1 = calculateHistogram(X01, IDString1, colour1, AngularResolutionIn);
	TH1D *hist2 = calculateHistogram(X02, IDString2, colour2, AngularResolutionIn);

	hist1->Sumw2();
	hist1->Divide(hist2);
	hist1->Draw();
	TF1 *f = new TF1("p", "pol0", 0.0, range);
	hist1->Fit(f, "QN");
	double probability = f->GetProb();
	hist1->Delete();
	hist2->Delete();
	return probability;

}
