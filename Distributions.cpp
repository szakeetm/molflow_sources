//#include <malloc.h>
//#include <cstring>
#include "Distributions.h"
#include <algorithm>
//#include "Random.h"
//#include "GLApp\GLTypes.h"

/*Distribution2D::Distribution2D(int N){
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
*/
std::vector<std::pair<double,double>> Generate_CDF(double gasTempKelvins,double gasMassGramsPerMol,size_t size){
	std::vector<std::pair<double,double>> cdf;cdf.reserve(size);
	double Kb=1.38E-23;
	double R=8.3144621;
	double a=sqrt(Kb*gasTempKelvins/(gasMassGramsPerMol*1.67E-27)); //distribution a parameter. Converting molar mass to atomic mass

	//Generate cumulative distribution function
	double mostProbableSpeed=sqrt(2*R*gasTempKelvins/(gasMassGramsPerMol/1000.0));
	double binSize=4.0*mostProbableSpeed/(double)size; //distribution generated between 0 and 4*V_prob
	/*double coeff1=1.0/sqrt(2.0)/a;
	double coeff2=sqrt(2.0/PI)/a;
	double coeff3=1.0/(2.0*pow(a,2));

	for (size_t i=0;i<size;i++) {
		double x=(double)i*binSize;
		cdf.push_back(std::make_pair(x,erf(x*coeff1)-coeff2*x*exp(-pow(x,2)*coeff3)));
	}*/
	for (size_t i = 0; i<size; i++) {
		double x = (double)i*binSize;
		double x_square_per_2_a_square = pow(x, 2) / (2 * pow(a, 2));
		cdf.push_back(std::make_pair(x, 1 - exp(-x_square_per_2_a_square)*(x_square_per_2_a_square + 1)));
	}

	/* //UPDATE: not generating inverse since it was introducing sampling problems at the large tail for high speeds
	//CDF created, let's generate its inverse
	std::vector<std::pair<double,double>> inverseCDF;inverseCDF.reserve(size);
	binSize=1.0/(double)size; //Divide probability to bins
	for (size_t i=0;i<size;i++) {
		double p=(double)i*binSize;
		//inverseCDF.push_back(std::make_pair(p,InterpolateX(p,cdf,TRUE)));
		inverseCDF.push_back(std::make_pair(p, InterpolateX(p, cdf, FALSE)));
	}
	return inverseCDF;
	*/
	return cdf;
}
double my_erf(double x)
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

BOOL compare_second(const std::pair<double,double>& lhs, const std::pair<double,double>& rhs) {
	return (lhs.second<rhs.second);
}


double InterpolateY(double x, const std::vector<std::pair<double, double>>& table, BOOL limitToBounds, BOOL logarithmic){
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size() == 1) return table[0].second; //constant value

	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	std::vector<std::pair<double, double> >::const_iterator lower, upper;
	bool outOfLimits = false;

	if (x >= table.back().first) {
		if (limitToBounds) return table.back().second;
		else {
			outOfLimits = true;
			lower = upper = table.end() - 1;
			lower--;
		}
	}
	else if (x < table[0].first) {
		if (limitToBounds) return table[0].second;
		else {
			outOfLimits = true;
			lower = upper = table.begin();
			upper++;
		}
	}

	// INFINITY is defined in math.h in the glibc implementation
	if (!outOfLimits) {
		lower = upper = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -MY_INFINITY));
		// Corner case
		if (upper == table.begin()) return upper->second;
		lower--;
	}
	if (logarithmic) return exp(log(lower->second) + (log(upper->second) - log(lower->second))
		*(log(x) - log(lower->first)) / (log(upper->first) - log(lower->first)));
	else return lower->second + (upper->second - lower->second)*(x - lower->first) / (upper->first - lower->first);

}

double InterpolateX(double y,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds){
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size()==1) return table[0].second; //constant value

	// Assumes that "table" is sorted by .second
    // Check if y is out of bound
    std::vector<std::pair<double, double> >::const_iterator lower, upper;
	BOOL outOfLimits = FALSE;
	
	if (y >= table.back().second) {
		if (limitToBounds) return table.back().first;
		else {
			outOfLimits = TRUE;
			lower = upper = table.end()-1;
			lower--;
		}
	} else if (y < table[0].second) {
		if (limitToBounds) return table[0].first;
		else {
			outOfLimits = TRUE;
			lower = upper = table.begin();
			upper++;
		}
	}

	// INFINITY is defined in math.h in the glibc implementation
	if (!outOfLimits) {
		lower = upper = std::lower_bound(table.begin(), table.end(), std::make_pair(MY_INFINITY, y),compare_second);
		// Corner case
		if (upper == table.begin()) return upper->first;
		lower--;
	}
	return lower->first + (upper->first - lower->first)*(y - lower->second)/(upper->second - lower->second);	
}

double FastLookupY(double x,const std::vector<std::pair<double,double>>& table,BOOL limitToBounds){
	//Function inspired by http://stackoverflow.com/questions/11396860/better-way-than-if-else-if-else-for-linear-interpolation
	_ASSERTE(table.size());
	if (table.size()==1) return table[0].second; //constant value

	// Assumes that table .first is SORTED AND EQUIDISTANT
    // Check if x is out of bound
    std::vector<std::pair<double, double> >::const_iterator lower, upper;
	BOOL outOfLimits = FALSE;

	if (x >= table.back().first) {
		if (limitToBounds) return table.back().second;
		else {
			outOfLimits = TRUE;
			lower = upper = table.end()-1;
			lower--;
		}
	} else if (x < table[0].first) {
		if (limitToBounds) return table[0].second;
		else {
			outOfLimits = TRUE;
			lower = upper = table.begin();
			upper++;
		}
	} 

	if (!outOfLimits) {
		double distanceX = table[1].first-table[0].first;
		size_t lowerIndex = (int)((x-table[0].first)/distanceX);
		lower = upper = table.begin()+(lowerIndex+1);
		// Corner case
		if (upper == table.begin()) return upper->second;
		lower--;
	}
	double result= lower->second + (upper->second - lower->second)*(x - lower->first)/(upper->first - lower->first);
	return result;
}