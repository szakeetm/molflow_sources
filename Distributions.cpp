#include <malloc.h>
#include <cstring>
#include <math.h>
#include "Distributions.h"
#include "Simulation.h"
#include "Random.h"
//#include "GLApp\GLTypes.h"

Distribution2D::Distribution2D(int N){
	if (!(N>0) && (N<10000000)) N=1; //don't create 0 size distributions
	valuesX=(double*)malloc(N*sizeof(double));
	valuesY=(double*)malloc(N*sizeof(double));
	size=N;
	average=0.0;
}

Distribution2D::Distribution2D(const Distribution2D &copy_src){ //copy constructor to avoid shallow copy
	valuesX=(double*)malloc(copy_src.size*sizeof(double));
	valuesY=(double*)malloc(copy_src.size*sizeof(double));
	memcpy(valuesX,copy_src.valuesX,copy_src.size*sizeof(double));
	memcpy(valuesY,copy_src.valuesY,copy_src.size*sizeof(double));
	size=copy_src.size;
	average=copy_src.average;
	average1=copy_src.average1;
}

Distribution2D::~Distribution2D(){
	SAFE_FREE(valuesX);
	SAFE_FREE(valuesY);
}

Distribution2D& Distribution2D::operator= (const Distribution2D &copy_src) {
	if (this != &copy_src) // protect against invalid self-assignment
	{
		//*this=Distribution2D_tiny(copy_src.size);
		valuesX=(double*)malloc(copy_src.size*sizeof(double));
		valuesY=(double*)malloc(copy_src.size*sizeof(double));
		memcpy(valuesX,copy_src.valuesX,copy_src.size*sizeof(double));
		memcpy(valuesY,copy_src.valuesY,copy_src.size*sizeof(double));
		size=copy_src.size;
	}
	// by convention, always return *this
	return *this;
}

double Distribution2D::InterpolateY(double x) {
	int inferior_index,superior_index;
	double slope, overshoot;

	for (superior_index=0;valuesX[superior_index]<x && superior_index<size;superior_index++);
	if (superior_index==size) superior_index--; //not found, x too large
	if (superior_index==0)    superior_index++; //not found, x too small
	inferior_index=superior_index-1;

	double diffX=valuesX[superior_index]-valuesX[inferior_index];
	double diffY=valuesY[superior_index]-valuesY[inferior_index];
	slope=diffY/diffX;
	overshoot=x-valuesX[inferior_index];

	return valuesY[inferior_index]+slope*overshoot;
}

double Distribution2D::InterpolateX(double y) {
	int inferior_index,superior_index;
	double slope, overshoot;

	for (superior_index=0;valuesY[superior_index]<y && superior_index<size;superior_index++);
	if (superior_index==size) return valuesX[size-1]; //not found, y too large
	if (superior_index==0)    return valuesX[0]; //not found, y too small
	inferior_index=superior_index-1;

	double diffX=valuesX[superior_index]-valuesX[inferior_index];
	double diffY=valuesY[superior_index]-valuesY[inferior_index];
	slope=diffY/diffX;
	if (slope==0.0) return valuesX[inferior_index];
	overshoot=y-valuesY[inferior_index];

	return valuesX[inferior_index]+overshoot/slope;
}

int Distribution2D::findXindex(double x) {
	int superior_index;
	for (superior_index=0;valuesX[superior_index]<x && superior_index<size;superior_index++);
	return superior_index;
}

Distribution2D Generate_CDF(double gasTempKelvins,double gasMassGramsPerMol,int size){
	
	double Kb=1.38E-23;
	double R=8.3144621;
	double a=sqrt(Kb*gasTempKelvins/(gasMassGramsPerMol*1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

	//Generate cumulative distribution function
	Distribution2D CDF(size);
	double mostProbableSpeed=sqrt(2*R*gasTempKelvins/(gasMassGramsPerMol/1000.0));
	double binSize=4.0*mostProbableSpeed/(double)size; //distribution generated between 0 and 3*V_prob
	double coeff1=1.0/sqrt(2.0)/a;
	double coeff2=sqrt(2.0/PI)/a;
	double coeff3=1.0/(2.0*pow(a,2));

	for (int i=0;i<size;i++) {
		double x=(double)i*binSize;
		CDF.valuesX[i]=x;
		CDF.valuesY[i]=erf(x*coeff1)-coeff2*x*exp(-pow(x,2)*coeff3);
	}

	/*
	//CDF created, let's generate its inverse
	Distribution2D inverseCDF(size);
	binSize=1.0/(double)size; //Divide probability to bins
	for (int i=0;i<size;i++) {
		double p=(double)i*binSize;
		inverseCDF.valuesX[i]=p;
		inverseCDF.valuesY[i]=CDF.InterpolateX(p);
	}*/
	return CDF;
}

double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}