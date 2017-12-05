#include "ScatteringEvent.h"
// Constructor for a ScatteringEvent
ScatteringEvent::ScatteringEvent(double X0In, double MaterialLengthIn, double NominalMomentumIn, double MomentumVarianceIn, double AngularResolutionIn){
	X0 = X0In;
	MaterialLength = MaterialLengthIn;
	AngularResolution = AngularResolutionIn;
	Momentum = calculateMomentum(NominalMomentumIn, MomentumVarianceIn);
	ScatteringVariance = calculateScatteringVariance(NominalMomentumIn);
	ScatteringAngleX = calculateScatteringAngle();
	ScatteringAngleY = calculateScatteringAngle();
	IncidenceAngleX = calculateIncidenceAngle();	
	IncidenceAngleY = calculateIncidenceAngle();	
	TotalScatterInX = calculateTotalScatter(ScatteringAngleX, IncidenceAngleX);
	TotalScatterInY = calculateTotalScatter(ScatteringAngleY, IncidenceAngleY);
	PathLengthInX = calculatePathLength(TotalScatterInX, MaterialLength);
	PathLengthInY = calculatePathLength(TotalScatterInY, MaterialLength);
	PathLengthTotal = calculateTotalPathLength(PathLengthInX, PathLengthInY);
	ScatteringStrength = calculateScatteringStrength(NominalMomentumIn);	
}
// Destructor for a ScatteringEvent
ScatteringEvent::~ScatteringEvent(){};

// Calculates Momentum
double ScatteringEvent::calculateMomentum(double NominalMomentumIn, double MomentumVarianceIn){
	TRandom3 *r = new TRandom3(0);
	double calculatedMomentum = r->Gaus(NominalMomentumIn,MomentumVarianceIn);
	return calculatedMomentum;
}

// Calculates the variances in material scattering
double ScatteringEvent::calculateScatteringVariance(double NominalMomentumIn){
	double Variance = ((13.6*13.6)/(NominalMomentumIn*NominalMomentumIn))*(MaterialLength/X0);
	return Variance;	
} 

// Calculates the Scattering in a plane
double ScatteringEvent::calculateScatteringAngle(){
	TRandom3 *r = new TRandom3(0);
	double ScatteringAngle = r->Gaus(0,ScatteringVariance);
	return ScatteringAngle;
}

// Calculates Incidence Angle in a plane
double ScatteringEvent::calculateIncidenceAngle(){
	TRandom3 *r = new TRandom3(0);
	double IncidenceAngle = r->Gaus(0,AngularResolution);
	return IncidenceAngle;
}

// Calculates the total Scattering in a plane (Incidence + Scattering)
double ScatteringEvent::calculateTotalScatter(double x1,double x2){
	double xtot = x1 + x2;
	return xtot;
}

// calculates the path length in a given plane
double ScatteringEvent::calculatePathLength(double TotalScatterIn, double MaterialLength){
	double PathLengthXY = MaterialLength/(cos(TotalScatterIn*TotalScatterIn)*cos(TotalScatterIn*TotalScatterIn));
	return PathLengthXY;
}

// calculates the total path length traversed
double ScatteringEvent::calculateTotalPathLength(double x, double y){
	double PathLength = sqrt(x*x + y*y);
	return PathLength;
}

// calculates the scattering strength
double ScatteringEvent::calculateScatteringStrength(double NominalMomentumIn){
	double a; 
	a = (Momentum*Momentum)/(NominalMomentumIn*NominalMomentumIn);
	double b; 
	b = TotalScatterInX*TotalScatterInX + TotalScatterInY*TotalScatterInY;
	double calcScatteringStrength = a*(b/(2*PathLengthTotal));
	return calcScatteringStrength;
}
// returns the scattering strength
double ScatteringEvent::getScatteringStrength(){
	return ScatteringStrength;
}

