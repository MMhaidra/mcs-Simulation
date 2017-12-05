//Prototype6
//This version changes the number of detectors to 4 and includes a magnetic field which accelerates the muon
//in between detectors 1 and 2 in order the measure the momentum of the muon

#include "TROOT.h"
#include <iostream>
#include "TFormula.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TH1.h"
#include <math.h>  
#include "TRandom.h"
using namespace std;

//Setting up simulation

//the b field
//compare the displacement over a small and large distance, so if the muon scatters how high up will it go.
//make plots of position at the detectors and see how they change when they scatter. Use a large distance and a small distance
//detector infinitley small (large scale plastic scintillator, small scale silicon) for now
//2 detectors before and after Bfield so you can see the path
//have detectro close to material to maximise resoloution
//have magnetic field in small section of the 50cm gap (see paper)
//z is beam direction x and y detecotr plane (x and y +- 5 cm)
//for each detection plane plot x and y scatter
//and also plot r and z
//make detector based around zero
//plot scattering angle on material using detectors
//plot expectation line from angle of first position of muon then plot line of best fit through all the points, use TProfile

 
//Muon Controls
double nomMom = 1250; //Nominal muon momentum in Mev // * momentum by this to get physical value.
double muonMass = 1; //mass of the muon 
int numMuons = 0; //Number of muons in the simulation (currently) //10 lookthousand muons
int totalMuons = 1000; //Total number of muons in the simulation
TVector3 MuonMom[totalMuons];
TVector3 MuonPos[totalMuons];
TVector3 MuonAccel[totalMuons]; //This will be zero unless the muon is inside the Bfield

double setUpLength = 1.2; //The length of the simulation in m used to place the material, detectors and Bfield inside

//Detector controls
int numDetectors = 6;//Number of detectors will most likely be 4 and always 4 in the real thing
TVector3 DetectorPos[numDetectors];
double detectorLengthZ = 0.02; // this corresponds to z value in m
double detectorLengthX =0.02; //x
double detectorWidth = 0.1; // this corresponds with the y in m
double X0Det = 9.370; //Radiation length of silicon in cm

//Material controls
int NumMaterials = 1; //number of materials the muon will scatter through (excluding detector);
TVector3 Material1Pos(0,0,0);
double materialLength = 0.05; // this corresponds to x and z value in m
double materialWidth = 0.1; // this corresponds with the y in m was 0.07
double X0Mat = 0.3166; //Radiation length of said material (iron) in cm

//Bfield controls
TVector3 Bfield(0,0,0); //Bfield in tesla
TVector3 BfieldPos(0,0,0); 
double BfieldLength = 1; // length of the field (y) in m
double BfieldWidth = 0.2; // width of the field (x and z) in m

//Simulation controls
double varianceMat = 0; //variance of the material
double varianceDet = 0; //variance of the detector
bool* measured = new bool[totalMuons](); //has the particle been measured
bool* scattered = new bool[totalMuons](); //has the particle scattered already
double time = 0;
double timestep = 0.0001;
int step = 0;  //used to spawn muons and record position for plots
int detNumber = 0; //detector number which a particle has been measured in

double numScatter = 0; //records the number of muons that scatter off the material

//create array used to create graphs 
//records the position of the muons which hit every detector
TVector3 muPosDet[numDetectors][totalMuons];


//variables used for z-plot
TVector3 zPlot[totalMuons][100];
TVector3 lineMuonPlot[100];
TVector3 lineMuonPos(0,0,0);
TVector3 lineMuonMom(0,0,0);
int index=0; //used for zPlot array




int mcsSimulation(){	
	setupSim(setUpLength);	
	varianceMat = calculateVariance(X0Mat);
	varianceDet = calculateVariance(X0Det);	
	
	
	
	for (double i = 0; time<0.12; ++i){	//0.06	
		//used to spawn muons
		if(numMuons != (totalMuons)){ 
			generateMuon(numMuons);
		}	

		lineMuonPos = MuonPos[0];
		lineMuonMom = MuonMom[0];
		

		for (int j = 0; j<numMuons; ++j){

			//records muon position
			if(step % 5 == 0){		
				//zPlot[j][index] = MuonPos[j];
				//lineMuonPlot[index] = lineMuonPos; 
			}

			//Check to see if muon can be scattered 
			checkBField(MuonPos[j],j);
			checkMeasure(MuonPos[j],j);		
			checkScatter(MuonPos[j],j);
			updateEulerC(j);
			//updateEulerCLine();
		}
		if(step % 5 == 0){
			index = index + 1;
		}
		//increase time		
		time = time + timestep;
		//update position of muon		
		
		step = step + 1;

		if(numScatter == numMuons){
			//cout << time << endl;
		}	
	}
	
	//prints the scattered angle after material, used to calculate distance between u and h2o
	//double angle = (muPosDet[2][0].Angle(muPosDet[3][0]));
	//cout << angle << endl;		
	cout << numMuons << endl;
	cout << numScatter << endl;
	plotGraphs();
}

//generates the properties of the muons
//this includes position, velocity
//and accel
void generateMuon(int num){
	
	//placeholder pos generator, eventually need to do it generated with a spread (5deg)	
	TVector3 pos(0,0,0);	
	gRandom->SetSeed(0);
	double MomentumZ = 0.1;  // Set z mom = 1
	double MomentumY = gRandom->Gaus(0.0,0.01);  // Generate y and x momenta	
	double MomentumX = gRandom->Gaus(0.0,0.01);	
	//double MomentumX, MomentumY = 0; //no spread
	TVector3 mom(MomentumX, MomentumY, (MomentumZ*nomMom));	
	MuonMom[num] = mom;
	MuonPos[num] = pos;
	if(numMuons != (totalMuons)){
		numMuons = numMuons + 1;
	}
	//cout << "Muon " << numMuons << " created." << endl;
	
}
//places the detectors, material and bfield in their respective place
//based on the length of the entire setup. Uses 6 detectors
void setupSim(double setupLength){
	
	DetectorPos[0].SetZ(setupLength*0.02); //places detector1 2 percent of the way in 0.02
	DetectorPos[1].SetZ(setupLength*0.05); //places detector2 4 percent of the way in 0.04
	BfieldPos.SetZ(setupLength*0.25); //places the bfield 25 percent of the way in  0.25
	DetectorPos[2].SetZ(setupLength*0.47); //places detector3 47 percent of the way in 0.47
	Material1Pos.SetZ(setupLength*0.5); //places the material at halfway in
	DetectorPos[3].SetZ(setupLength*0.60); //places detector4 53 percent of the way in 0.53
	DetectorPos[4].SetZ(setupLength*0.95); //places detector5 96 percent of the way in 0.96
	DetectorPos[5].SetZ(setupLength*0.98); //places detector6 98 percent of the way in 0.98
	
	for (int i = 0; i<numDetectors; ++i){
		//cout << "detector " << i + 1 << " placed at" << DetectorPos[i].Z() << endl;

	}
	//cout << "sim set up" <<endl;
	
}


bool checkRangeDet(TVector3 MuPos){
	
	for (int i = 0; i<numDetectors; ++i){
		//TRYING TO CENTRE THE SIMULATION AROUND ZERO
		if((MuPos.X() >= (DetectorPos[i].X() - (detectorLengthX*0.5)) && MuPos.X() <= (DetectorPos[i].X() + (detectorLengthX*0.5))) &&
		(MuPos.Y() >= (DetectorPos[i].Y() - (detectorWidth*0.5)) && MuPos.Y() <= (DetectorPos[i].Y() + 	(detectorWidth*0.5))) &&
		(MuPos.Z() >= (DetectorPos[i].Z() - (detectorLengthZ*0.5)) && MuPos.Z() <= (DetectorPos[i].Z() + (detectorLengthZ*0.5)))) {
			
			detNumber = i;
			return true;
		}

	}

return false;
}

bool checkRangeBfield(TVector3 MuPos){
	
	if(MuPos.X() >= (BfieldPos.X() - BfieldWidth/2) && MuPos.X() <= (BfieldPos.X() + BfieldWidth/2) &&
	(MuPos.Y() >= (BfieldPos.Y() - BfieldLength/2) && MuPos.Y() <= (BfieldPos.Y() + BfieldLength/2)) &&
	(MuPos.Z() >= (BfieldPos.Z() - BfieldWidth/2) && MuPos.Z() <= (BfieldPos.Z() + BfieldWidth/2))) {
		
			return true;
	}



return false;
}

bool checkRangeScat(TVector3 MuonPos){

	if(MuonPos.X() >= (Material1Pos.X() - materialLength/2) && MuonPos.X() <= (Material1Pos.X() + materialLength) &&
	(MuonPos.Y() >= (Material1Pos.Y() - materialWidth/2) && MuonPos.Y() <= (Material1Pos.Y() + materialWidth)) &&
	(MuonPos.Z() >= (Material1Pos.Z() - materialLength/2) && MuonPos.Z() <= (Material1Pos.Z() + materialLength))) {

		return true;
	}

return false;
}
	

double calculateVariance(double L0){
	TFormula *variance = new TFormula("variance", "((13.6*13.6)/([0]*[0]))*([1]/x)");
	variance->SetParameter(0,nomMom);  
	variance->SetParameter(1,materialWidth);
	double VarianceIn = variance->Eval(L0);
	//cout << VarianceIn <<endl;
	return VarianceIn;
}
//Checks to see if the muon is in range of the Bfield
//if it is it changes its acceleration to that produced by lorentz force
//if it isnt it sets the acceleration = 0
void checkBField(TVector3 pos, int MuonNum){

	TVector3 noB(0,0,0);	
	if(checkRangeBfield(pos) == true){
		MuonAccel[MuonNum] = lorentz(MuonNum);
		//cout << "Muon " << MuonNum << "in B-field" << endl;
	}
	else{
		MuonAccel[MuonNum] = noB;
	}

}
//Checks to see if muon is in bounds of the detector
//Prints out the position of the muon if it is
//this is to check and see what value of x the muon
//has when it should be that of the detectors.
void checkMeasure(TVector3 pos, int MuonNum){
	//x is greater than as muon is coming from positive x direction, 
	
	if(checkRangeDet(pos) == true) {

		if(measured[MuonNum] == false){		
		printMuonPos(MuonNum);	
		measured[MuonNum] = true;
			
			if(scattered[MuonNum] == false){			
				simulateMCS(MuonNum, varianceDet);
				scattered[MuonNum] = true;
			}	
		}			
	}
	else{
		measured[MuonNum] = false;
	}
	
}

//prints the position of the muon to the console
//also records the position and which detector detected
// the muon in a 2x2 array
void printMuonPos(int MuonNum){
	//cout << "Muon " << MuonNum << " has been measured in detector " << detNumber << endl;	
	//cout << MuonPos[MuonNum].X() << " x" <<endl;
	//cout << MuonPos[MuonNum].Y() << " y" <<endl;
	//cout << MuonPos[MuonNum].Z() << " z" <<endl;
	muPosDet[detNumber][MuonNum] = MuonPos[MuonNum]; 
}

//Checks to see if muon is in bounds of the Material
//Prints out the position of the muon if it is
//this is to check and see what value of x the muon
//has when it should be that of the Materials.
void checkScatter(TVector3 pos, int MuonNum){
	//x is greater than as muon is coming from positive x direction, 
	
	if(checkRangeScat(pos) ==true) {

		if(scattered[MuonNum] == false){			
			simulateMCS(MuonNum, varianceMat);
			scattered[MuonNum] = true;
			numScatter = numScatter + 1;
		}		
		
	}
	
	else{
		scattered[MuonNum] = false;
	}
	
}
void simulateMCS(int MuonNum, double var){
	double ScatteringAnglex = gRandom->Gaus(0,var);  // Get ScatteringAngle_x
	double ScatteringAngley = gRandom->Gaus(0,var);  // Get ScatteringAngle_y
	double ScatteringAnglez = gRandom->Gaus(0,var);	 // Get ScatteringAngle_z

	MuonMom[MuonNum].RotateX(ScatteringAnglex); //rotate the velocity through the generated scattered angles
	MuonMom[MuonNum].RotateY(ScatteringAngley);
	MuonMom[MuonNum].RotateZ(ScatteringAnglez);
	
	//cout << "Particle "<< MuonNum << " scattered at time = " << time << endl;


}
//function which calculates the acceleration of the muon
//in the magnetic field
TVector3 lorentz(int num){
	TVector3 E(0,0,0);
	TVector3 vxB(0,0,0);
	vxB = MuonMom[num].Cross(Bfield);
	E = E + vxB;
	//E *= 1/MuonMass;
	return E;
}


//uses the Euler Cromer method to update the MuonPos of the particle
void updateEulerC(int MuonNum){
	TVector3 accel = MuonAccel[MuonNum];
	TVector3 tempVel = MuonMom[MuonNum];
	accel *= timestep;
	MuonMom[MuonNum] = tempVel + accel; //v = u + at
	MuonPos[MuonNum] = (MuonPos[MuonNum] + (tempVel *= timestep)); //old MuonPos + ut
}

//uses the Euler Cromer method to update lineMuonPos
void updateEulerCLine(){
	TVector3 accel = MuonAccel[0];
	TVector3 tempVel = lineMuonMom;
	accel *= timestep;
	lineMuonMom = tempVel + accel; //v = u + at
	lineMuonPos = (lineMuonPos + (tempVel *= timestep)); //old MuonPos + ut
}



void plotGraphs(){

//Plot x y spread on each detector	
	TCanvas *c1 = new TCanvas("c1", "c1", 200,200, 2000,2000);
	c1->Divide(2,3);
	TCanvas *c2 = new TCanvas("c2", "c2", 200,200, 2000,2000);
	c2->Divide(2,3);
	TH2* h = new TH2D(
		/* name */ "h",
		/* title */ "Muon position in 2D detector 1",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);	
	for (int i = 0; i<totalMuons; ++i){
		h->Fill(muPosDet[0][i].X(), muPosDet[0][i].Y());
	}
	
	TH2* h1 = new TH2D(
		/* name */ "h1",
		/* title */ "Muon position in 2D detector 2",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);	
	for (int i = 0; i<totalMuons; ++i){
		h1->Fill(muPosDet[1][i].X(), muPosDet[1][i].Y());
	}
	
	TH2* h2 = new TH2D(
		/* name */ "h2",
		/* title */ "Muon position in 2D detector 3",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);	
	for (int i = 0; i<totalMuons; ++i){
		h2->Fill(muPosDet[2][i].X(), muPosDet[2][i].Y());
	}
	
	TH2* h3 = new TH2D(
		/* name */ "h3",
		/* title */ "Muon position in 2D detector 4",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);		
	for (int i = 0; i<totalMuons; ++i){
		h3->Fill(muPosDet[3][i].X(), muPosDet[3][i].Y());
	}
	
	TH2* h4 = new TH2D(
		/* name */ "h4",
		/* title */ "Muon position in 2D detector 5",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);		
	for (int i = 0; i<totalMuons; ++i){
		h4->Fill(muPosDet[4][i].X(), muPosDet[4][i].Y());
	}
	
	TH2* h5 = new TH2D(
		/* name */ "h5",
		/* title */ "Muon position in 2D detector 6",
		/* x-axis */ 100,-(detectorLengthX/2), detectorLengthX/2,
		/* y-axis */ 100, -(detectorWidth/2), detectorWidth/2);	
	for (int i = 0; i<totalMuons; ++i){
		h5->Fill(muPosDet[5][i].X(), muPosDet[5][i].Y());
	}
	
	
/*
	c1->cd(1);
	h->Draw();
	h->GetXaxis()->SetTitle("x/m");
	h->GetYaxis()->SetTitle("y/m");	
	h->SetMarkerStyle(2);
	h->SetMarkerSize(2);

	c1->cd(2);
	h1->Draw();
	h1->SetMarkerStyle(2);
	h1->SetMarkerSize(2);
	h1->GetXaxis()->SetTitle("x/m");
	h1->GetYaxis()->SetTitle("y/m");	

	c1->cd(3);
	h2->Draw();
	h2->SetMarkerStyle(2);
	h2->SetMarkerSize(2);
	h2->GetXaxis()->SetTitle("x/m");
	h2->GetYaxis()->SetTitle("y/m");	

	c2->cd(1);
	h3->Draw();
	h3->SetMarkerStyle(2);
	h3->SetMarkerSize(2);
	h3->GetXaxis()->SetTitle("x/m");
	h3->GetYaxis()->SetTitle("y/m");	

	c2->cd(2);
	h4->Draw();
	h4->SetMarkerStyle(2);
	h4->SetMarkerSize(2);
	h4->GetXaxis()->SetTitle("x/m");
	h4->GetYaxis()->SetTitle("y/m");	

	c2->cd(3);
	h5->Draw();
	h5->SetMarkerStyle(2);
	h5->SetMarkerSize(2);
	h5->GetXaxis()->SetTitle("x/m");
	h5->GetYaxis()->SetTitle("y/m");	
*/
//Calculate angle between position vectors before and after b-field and material scattering

	TCanvas *c3 = new TCanvas("c3", "c3", 200,200, 2000,2000);
	c3->Divide(1,2);

	TH1D* g = new TH1D("g", "Scattered angle after b-field", totalMuons/10, 0,1e-5);	

	for (int i = 0; i<totalMuons; ++i){
		double angle = (muPosDet[1][i].Angle(muPosDet[2][i]));	
		//cout << angle << endl;		
		g->Fill(angle);
		
	}
	
	//cout << " " << endl;
	TH1D* U = new TH1D("U", "Scattered angle after material collision", totalMuons/10, 0,5e-5);	//changed from 10e-6
	double angle2 = 0;
	for (int i = 0; i<totalMuons; ++i){
		double angle = (muPosDet[2][i].Angle(muPosDet[3][i]));			
		
		if(angle > angle2){
			angle2 = angle;
		}		
		
		U->Fill(angle);
	}
	cout << "angle " << angle2 << endl;

	c3->cd(1);
	g->Draw();
	g->GetXaxis()->SetTitle("Scattered angle in rads");
	g->GetYaxis()->SetTitle("Events");
	gPad->SetLogy();

	c3->cd(2);
	U->GetXaxis()->SetTitle("Scattered angle in rads");
	U->GetYaxis()->SetTitle("Events");	
	U->Draw();
	gPad->SetLogy();

//Plot r against z graph
/*
 	TCanvas *c4 = new TCanvas("c4", "c4", 200,200, 2000,2000);
	TH2* z = new TH2D(
		/* name */ //"z",
		/* title */ //"r against z for muon 1",
		/* x-axis */ //100 ,0, 1.25,
		/* y-axis */ //100, -0.0002, 0.0002);	
/*
	for (int i = 0; i<index; ++i){
		double r = sqrt ((zPlot[0][i].X()*zPlot[0][i].X()) + (zPlot[0][i].Y()*zPlot[0][i].Y()));		
		z->Fill(zPlot[0][i].Z(), r);
		//cout << r << endl;
	}

	
	TH2* z1 = new TH2D(

*/		/* name */ //"z1",
		/* title */ //"r against z for muon 1",
		/* x-axis */ //100 ,0, 1.25,
		/* y-axis */ //100, -0.0002, 0.0002);	
/*	for (int i = 0; i<index; ++i){
		double r = sqrt ((lineMuonPlot[i].X()*lineMuonPlot[i].X()) + (lineMuonPlot[i].Y()*lineMuonPlot[i].Y()));		
		z1->Fill(lineMuonPlot[i].Z(), r);
		//cout << r << endl;
	}
/*	
	c4->cd(1);
	z->SetMarkerStyle(2);
	z->SetMarkerSize(2);
	z->GetXaxis()->SetTitle("z/m");
	z->GetYaxis()->SetTitle("r/m");
	z->Draw("C");
	z1->Draw("same");
	z1->SetMarkerColor(2);
	z1->SetMarkerSize(2);
	z1->SetMarkerStyle(2);
	
	//use other file to plot path of others muons and change line colour	
*/
	
//Plot 1DH of R at detector 4
	TCanvas *c5 = new TCanvas("c5", "c5", 200,200, 2000,2000);
	TH1D* f = new TH1D("f", "R at detector 5", totalMuons/2, 0,0.015); //0.00023
	double r2 = 0;
	for (int i = 0; i<totalMuons; ++i){
		double r = sqrt ((muPosDet[4][i].X()*muPosDet[4][i].X())+(muPosDet[4][i].Y()*muPosDet[4][i].Y()));		
		
		if(r > r2){
			r2 = r;
		}		
		
		f->Fill(r);
	}
	cout << "r " << r2 << endl;	

	c5->cd(1);
	f->Draw();
	f->GetXaxis()->SetTitle("R/m");
	f->GetYaxis()->SetTitle("Events");
	
}



