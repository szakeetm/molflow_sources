#ifndef _DISTRIBUTIONS_
#define _DISTRIBUTIONS_
/*#define NUMBER_OF_DISTRO_VALUES 100
#define NUMBER_OF_INTEGR_VALUES 250
#define VERY_SMALL 1.0E-30
#define UPPER_LIMIT 100.0
#define LOWER_LIMIT 1.0E-10

#define INTEGRAL_MODE_N_PHOTONS 1
#define INTEGRAL_MODE_SR_POWER  2

#define SYNGEN_MODE_FLUXWISE  0
#define SYNGEN_MODE_POWERWISE 1*/

//#include "File.h" //FileReader for LoadCSV
#include <vector>
#include "Simulation.h" //for BOOL and FALSE macros

double my_erf(double x);
double InterpolateY(double x,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds=FALSE,BOOL logarithmic=FALSE);
double InterpolateX(double y,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds=FALSE);
double FastLookupY(double x,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds=FALSE);

/*class Averages {
public:
	double average;
	double average1;
	double average_;
};*/

/*class Distribution2D {
public:
	Distribution2D(int size);
	Distribution2D::Distribution2D(const Distribution2D &copy_src); //copy constructor
	Distribution2D& operator= (const Distribution2D & other); //assignment op
	~Distribution2D();
	double InterpolateY(double x); //interpolates the Y value corresponding to X (allows extrapolation)
	double InterpolateX(double y); //(no extrapolation, first/last X values are the output limits)
	double *valuesX,*valuesY;
	int findXindex(double x);
	int size;
	double average,average1;
	double Interval_Mean(double x1,double x2);
};

class Indexes {
public:
	double inferior_flux,inferior_power;
	double superior_flux,superior_power;
	int i1;
};*/

/*
double g0ki(double x, double order, int kind);
double SYNRAD_FAST(double x);
double Gi(double x,int order);
double H(double x, int order);
double calc_polarization_percentage(double energy,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
double find_psi(double x,double gamma_square,double f_times_g1h2,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
double find_chi(double psi,double gamma_square,double f_times_g1h2,bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
double SYNGEN1(double x_min,double x_max,int mode);

Distribution2D Generate_K_Distribution(double order);
Distribution2D Generate_G1_H2_Distribution();
Distribution2D Generate_LN_Distribution(); //precalculated ln(x) values for the most used energies
Distribution2D Generate_Polarization_Distribution(bool calculate_parallel_polarization, bool calculate_orthogonal_polarization);
Distribution2D Generate_Integral(double x1,double x2,int mode);*/

#endif