#include "ScatteringHistogram.h"
#include "TROOT.h"
#include "TH1.h"
#include <string>

void ScatteringEventHistogram() {
	int N = 1000;
	int bins = 100;
	double range = 2;
	double MaterialLength = 0.1;
	double NominalMomentum = 3000;
	double MomentumVariance = 500;
	double AngularResolution = 0.00001;
	char *Material1 = "Iron";
	char *Material2 = "Lead";
	double totalprob;
	double probabilities[50];
	ScatteringHistogram *Default = new ScatteringHistogram(N, bins, range, MaterialLength, NominalMomentum, MomentumVariance);
	int i = 0;
	
	//if (AngularResolution>0){		
	for (int k = 0; k < 20;++k){
		for (int j = 0; j < 50; ++j){
			double probability = Default->calculateProbability(1.313, 0.5612, Material1, Material2, 1, 2, AngularResolution);
			totalprob = totalprob + probability;
		}

		cout << k << endl;
		AngularResolution = AngularResolution-0.000002;
		probabilities[k] = 1-totalprob/50;
		totalprob = 0;
	}
	
	for (int k = 0; k < 20; ++k) {
		cout <<  probabilities[k] * 100 << "%" << endl;
	}


}