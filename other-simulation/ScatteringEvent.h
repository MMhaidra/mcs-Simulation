#ifndef ScatteringEvent_h
#define ScatteringEvent_h

#include <ctime>
#include "TRandom3.h"

class ScatteringEvent{
    public:
        ScatteringEvent(double X0In, double MaterialLengthIn, double NominalMomentumIn, double MomentumVarianceIn, double AngularResolutionVarianceIn);
        ~ScatteringEvent();
        double getScatteringStrength();
     
        
    private:
        TRandom3 *r;        
        double X0;
		double AngularResolution;
        double MaterialLength;        
        double Momentum;        
        double IncidenceAngleX;
        double IncidenceAngleY;
        double ScatteringAngleX;
        double ScatteringAngleY;
        double ScatteringStrength;
        double ScatteringVariance;
        double TotalScatterInX;
        double TotalScatterInY;
        double PathLengthInX;
        double PathLengthInY;
        double PathLengthTotal;
        
        double calculateMomentum(double NominalMomentumIn, double MomentumVarianceIn);
        double calculateScatteringVariance(double NominalMomentumIn);
        double calculateScatteringAngle();
        double calculateIncidenceAngle();
        double calculateTotalScatter(double x1,double x2);
        double calculatePathLength(double TotalScatterIn, double MaterialLength);
        double calculateTotalPathLength(double x, double y);
        double calculateScatteringStrength(double NominalMomentumIn);
};
#endif
