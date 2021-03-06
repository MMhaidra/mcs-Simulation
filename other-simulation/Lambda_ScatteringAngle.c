// ================================================================================
// Author: Daniel Parry
// Current Version: 2.1
// Date: 30-Jul-2015
// 
// Description: The aim of this code is to plot a histogram of 
// 				(p/p_0)^2(ScatteringAngle_x^2 + ScatteringAngle_y^2)/2L,  where p is the momentum
// 				of the muon, ScatteringAngle_x/y is the change in angle for the muon, p0 is the
// 				nominal muon momentum and L is the path lenth that the muon traverses.
// ================================================================================
//
// ================================================================================
// ============================== CHANGE LOG ======================================
// ================================================================================
// Date: 21-Sep-2015
// Version: 2.2
// Description: Path length was calculated wrong, this has been resolved.
// ================================================================================
// Date: 20-Aug-2015
// Version: 2.1
// Description: - Changed the code so that the axis scales with material length
//				  and that the variance also uses material length.
// ================================================================================
// Date: 13-Aug-2015
// Version: 2.0
// Description: - Major overhaul to the code.
// 				- Previously all random values we taken from histograms, now random
//				  number generators are used.
//				- The code no longer displays theta as a plot 
// 				- Unit issues were fixed so values displayed are correct with 
//				  the units displayed on each plot and in comments
//				- Inclusion of muon flux is now also displayed
//				- Most variables and methods have been renamed to make more sense.
//				- Path length is now calculated correctly.
// ================================================================================

#include "TROOT.h"
#include "TH1.h"
#include <string>
#include "TFormula.h"
#include <iostream>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
using namespace std;

// ================================================================================
// This method uses equation 5 from the article Image Reconstruction and Material Z
// discrimination via cosmic ray muon radiagraphy by L.J Schultz et al. (2003) 
// to get the variance of the scattering angle, uses the constant 13.6^2 and 
// divides it by the nominal muon momentum  NomMom^2 then divides by the 
// radiation length L0 of the material in question. It is worth noting that this
// per unit depth so requires an input of material length of which it is multiplied
// by.

// ================================================================================
double Variance(double L0, double NomMom, double LMat) {
	TFormula *variance = new TFormula("variance", "((13.6*13.6)/([0]*[0]))*([1]/x)");
	variance->SetParameter(0,NomMom);  
	variance->SetParameter(1,LMat);
	double VarianceIn = variance->Eval(L0);
	cout << VarianceIn <<endl;
	return VarianceIn;
}

// ==================================================================================
// This method uses equation 6 from article Image Reconstruction and Material Z
// discrimination via cosmic ray muon radiagraphy by L.J Schultz et al. (2003) 
// to plot a histogram of the scattering strenght of a material for N events.
// It disregards the 1/N and summation. Where NominalMom the the nominal muon momentum,
// MomVariance is the nominal muon momentum variance, Variance is the variance in
// scattering angle and LPath is the path length that the muon traverses.
// ==================================================================================

TH1D *LambdaVertical(double NominalMom, double MomVariance, double Variance, double LMat, int N, const char* hname="Lambda"){
  	char titlestr[128];
  	strcpy( titlestr, hname );
  	strcat( titlestr, " vs Events");
	double XaxisRange = 3e-8*LMat; //changed from 3e-8
	//	TRandom3 *r = new TRandom3(0);
	TH1D *LambdaV = new TH1D(hname, titlestr, N/10, 0,XaxisRange);
	for (int i = 0; i<N; ++i){
		
		double Momentum = gRandom->Gaus(NominalMom,MomVariance); //Get a value for Muon Momentum
		
		double ScatteringAnglex = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_x
		double ScatteringAngley = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_y
		
		double Lambda = ((Momentum*Momentum)/(NominalMom*NominalMom))*(((ScatteringAnglex*ScatteringAnglex)+(ScatteringAngley*ScatteringAngley))/(2*LMat));
		LambdaV -> Fill(Lambda);
	}	
	return LambdaV;
}

// ==================================================================================
// This method uses a similar proceedure seen in VerticalFlux except it includes the
// muon flux when calculating scattering angle. To do this it takes a two random 
// values from a Gaussian with a mean at 0 and variance of 0.56 (which is a similar 
// distribution to cos^2), then it adds these values to the scattering angles in x and y
// then works out Path Length by sqrt((L/cos(thetax))^2 + (L/cos(thetay))^2).
// ==================================================================================

TH1D *LambdaFlux(double NominalMom, double MomVariance, double Variance, double LMat, int N, const char* hname="Lambda"){
  	char titlestr[128];
  	strcpy( titlestr, hname );
  	strcat( titlestr, " vs Events");
	//	TRandom3 *r = new TRandom3(0);
	double XaxisRange = 2./LMat; //changed from 2.
	
	TH1D *LambdaR = new TH1D(hname, titlestr, N/100, 0, XaxisRange);
	for (int i = 0; i<N; ++i){
		double Momentum = gRandom->Gaus(NominalMom,MomVariance);
		
		double ScatteringAnglex = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_x currently mean = 0 
		double ScatteringAngley = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_y
		
		double IncidenceAnglex = gRandom->Gaus(0,0.56); //changed from 0
		double IncidenceAngley = gRandom->Gaus(0,0.56); //changed from 0
		
		double ObsAnglex = ScatteringAnglex + IncidenceAnglex;
		double ObsAngley = ScatteringAngley + IncidenceAngley;
		
		double PathLengthx = LMat/(cos(ObsAnglex)*cos(ObsAnglex));
		double PathLengthy = LMat/(cos(ObsAngley)*cos(ObsAngley));
		
		double PathLength = sqrt(PathLengthx*PathLengthx + PathLengthy*PathLengthy);
		
		double Lambda = ((Momentum*Momentum)/(NominalMom*NominalMom))*(((ObsAnglex*ObsAnglex)+(ObsAngley*ObsAngley))/(2*PathLength));
		LambdaR -> Fill(Lambda);
	}	
	return LambdaR;
}

// ==================================================================================
// This method is identical to LambdaVertical except it includes the uncertainity
// in the momentum in the final calculation of lambda.
// ==================================================================================

TH1D *LambdaVertUncert(double NominalMom, double MomVariance, double MomUncert, double Variance, double LMat, int N, const char* hname="Lambda"){
  	char titlestr[128];
  	strcpy( titlestr, hname );
  	strcat( titlestr, " vs Events (with momentum uncertainity");
	double XaxisRange = 3e-8*LMat; //changed from 3e-8
	//	TRandom3 *r = new TRandom3(0);
	TH1D *LambdaVU = new TH1D(hname, titlestr, N/10, 0,XaxisRange);
	for (int i = 0; i<N; ++i){
		
		double Momentum = gRandom->Gaus(NominalMom,MomVariance); //Get a value for Muon Momentum
		
		double ScatteringAnglex = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_x
		double ScatteringAngley = gRandom->Gaus(0,Variance);  // Get ScatteringAngle_y

		double E_p = ((MomUncert*Momentum)/Momentum); //Calculate E_p 
		
		double Lambda1 = (1/((1+(E_p*E_p))));
		double Lambda2 = ((Momentum*Momentum)/(NominalMom*NominalMom))*(((ScatteringAnglex*ScatteringAnglex)+(ScatteringAngley*ScatteringAngley))/(2*LMat));
		
		double Lambda3 = (Lambda1*Lambda2);
		//double Lambda = (1/(N(1+(E_p*E_p))))(((Momentum*Momentum)/(NominalMom*NominalMom))*(((ScatteringAnglex*ScatteringAnglex)+(ScatteringAngley*ScatteringAngley))/(2*LMat)));
		

		//if(i==5) {
  			//cout << Lambda3 <<endl; 
		//} 
				
		LambdaVU -> Fill(Lambda3);
	}	
	return LambdaVU;
}

// ==================================================================================
// Main Method to plot Scattering Strength as a histogram.
// ==================================================================================

void Lambda_ScatteringAngle(){

//Creates the Canvas
TCanvas *c1 = new TCanvas("c1", "c1", 200,200, 1200,900);; //changed from TCanvas *c1 = new TCanvas("c1", "c1", 200,200, 800,600);
c1->Divide(1,3);

//Sets the Number of Events
double N = 100000;

//Set of material in Z plane 
double MaterialLength = 10000; // in cm

// Define Radiation Lengths of Materials
double L0_U  = 0.3166;  //Radiation Length of Uranium in cm
double L0_Pb = 0.5612;  //Radiation Length of Lead in cm
double L0_Fe = 1.757;  //Radiation Length of Iron in cm
double L0_H2O = 39.31; //Radiation Length of Water in cm

// Defines Mometum Properties
double NominalMuonMomentum = 3000; // Nominal Muon Momentum in MeV 
double MomentumVariance = 500; // Nominal Muon Momentum Variance in MeV	
double MomentumUncertainty = 0.1; //Uncertainty of Muon momentum as a percentage

// Evaluates the Variance of each Material
double Variance_Pb = Variance(L0_Pb, NominalMuonMomentum, MaterialLength); // Lead
double Variance_Fe = Variance(L0_Fe, NominalMuonMomentum, MaterialLength); // Iron
double Variance_U = Variance(L0_U, NominalMuonMomentum, MaterialLength);   // Uranium
double Variance_H2O = Variance(L0_H2O, NominalMuonMomentum, MaterialLength); // Water

//Create Histograms for the Scattering Strength of Lead where V denotes vertical, R denote Real (or the inclusion of Flux) and VU denotes vertical with uncertainties
TH1D *VScatteringStrength_Pb = LambdaVertical(NominalMuonMomentum, MomentumVariance,  Variance_Pb, MaterialLength, N, "Scattering Strength of Pb" );
VScatteringStrength_Pb->SetLineColor(2);
TH1D *RScatteringStrength_Pb = LambdaFlux(NominalMuonMomentum, MomentumVariance,  Variance_Pb, MaterialLength, N, "Scattering Strength of Pb (including Muon Flux)" );
RScatteringStrength_Pb->SetLineColor(2);
TH1D *VUScatteringStrength_Pb = LambdaVertUncert(NominalMuonMomentum, MomentumVariance, MomentumUncertainty,  Variance_Pb, MaterialLength, N, "Scattering Strength of Pb (with uncertainties)"  );
VUScatteringStrength_Pb->SetLineColor(2);

//Create Histograms for the Scattering Strength of Iron where V denotes vertical, R denote Real (or the inclusion of Flux) and VU denotes vertical with uncertainties
TH1D *VScatteringStrength_Fe = LambdaVertical(NominalMuonMomentum, MomentumVariance,  Variance_Fe, MaterialLength, N, "Scattering Strength of Fe" );
VScatteringStrength_Fe->SetLineColor(3);
TH1D *RScatteringStrength_Fe = LambdaFlux(NominalMuonMomentum, MomentumVariance,  Variance_Fe, MaterialLength, N, "Scattering Strength of Fe (including Muon Flux)" );
RScatteringStrength_Fe->SetLineColor(3);
TH1D *VUScatteringStrength_Fe = LambdaVertUncert(NominalMuonMomentum, MomentumVariance, MomentumUncertainty,  Variance_Fe, MaterialLength, N, "Scattering Strength of Fe (with uncertainties)" );
VUScatteringStrength_Fe->SetLineColor(3);

//Create Histograms for Scattering Strength of Uranium where V denotes vertical, R denote Real (or the inclusion of Flux) and VU denotes vertical with uncertainties
TH1D *VScatteringStrength_U = LambdaVertical(NominalMuonMomentum, MomentumVariance,  Variance_U, MaterialLength, N, "Scattering Strength of Uranium" );
VScatteringStrength_U->SetLineColor(4);
TH1D *RScatteringStrength_U = LambdaFlux(NominalMuonMomentum, MomentumVariance,  Variance_U, MaterialLength, N, "Scattering Strength of Uranium (including Muon Flux)" );
RScatteringStrength_U->SetLineColor(4);
TH1D *VUScatteringStrength_U = LambdaVertUncert(NominalMuonMomentum, MomentumVariance, MomentumUncertainty,  Variance_U, MaterialLength, N, "Scattering Strength of U (with uncertainties)" );
VUScatteringStrength_U->SetLineColor(4);

//Create Histograms for Scattering Strength of Water where V denotes vertical, R denote Real (or the inclusion of Flux) and VU denotes vertical with uncertainties
TH1D *VScatteringStrength_H2O = LambdaVertical(NominalMuonMomentum, MomentumVariance,  Variance_H2O, MaterialLength, N, "Scattering Strength of Water" );
VScatteringStrength_H2O->SetLineColor(5);
TH1D *RScatteringStrength_H2O = LambdaFlux(NominalMuonMomentum, MomentumVariance,  Variance_H2O, MaterialLength, N, "Scattering Strength of Water (including Muon Flux)" );
RScatteringStrength_H2O->SetLineColor(5);
TH1D *VUScatteringStrength_H2O = LambdaVertUncert(NominalMuonMomentum, MomentumVariance, MomentumUncertainty,  Variance_H2O, MaterialLength, N, "Scattering Strength of Water (with uncertainties)" );
VUScatteringStrength_H2O->SetLineColor(5);

//Create Legend
TLegend* legLambda =new TLegend(0.4,0.6,0.7,0.8); 
legLambda->AddEntry(VScatteringStrength_U, "U", "l");
legLambda->AddEntry(VScatteringStrength_Pb, "Pb", "l");
legLambda->AddEntry(VScatteringStrength_Fe, "Fe", "l");
legLambda->AddEntry(VScatteringStrength_H2O, "Water", "l");
legLambda->Draw();

//Draw the Histograms on to Canvas'

c1->cd(1);
VScatteringStrength_Fe->Draw();
VScatteringStrength_Fe->GetXaxis()->SetTitle("Scattering Strength in Rad^2/cm");
VScatteringStrength_Fe->GetYaxis()->SetTitle("Events");
VScatteringStrength_Fe->SetTitle("Scattering Strength vs Events ");
VScatteringStrength_Pb->Draw("same");
VScatteringStrength_U->Draw("same");
VScatteringStrength_H2O->Draw("same");
gPad->SetLogy();
legLambda->Draw();

c1->cd(2);
RScatteringStrength_Fe->Draw();
RScatteringStrength_Fe->GetXaxis()->SetTitle("Scattering Strength in Rad^2/cm");
RScatteringStrength_Fe->GetYaxis()->SetTitle("Events");
RScatteringStrength_Fe->SetTitle("Scattering Strength vs Events (including Muon Flux)");
RScatteringStrength_Pb->Draw("same");
RScatteringStrength_U ->Draw("same");
RScatteringStrength_H2O->Draw("same");
gPad->SetLogy();
legLambda->Draw();


c1->cd(3);
VUScatteringStrength_Fe->Draw();
VUScatteringStrength_Fe->GetXaxis()->SetTitle("Scattering Strength in Rad^2/cm");
VUScatteringStrength_Fe->GetYaxis()->SetTitle("Events");
VUScatteringStrength_Fe->SetTitle("Scattering Strength vs Events (Including uncertainty) ");
VUScatteringStrength_Pb->Draw("same");
VUScatteringStrength_U->Draw("same");
VUScatteringStrength_H2O->Draw("same");
gPad->SetLogy();
legLambda->Draw();
}

int main(){
   Lambda_ScatteringAngle();
   return 0;
}

