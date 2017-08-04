#ifndef _PARAMETERH_
#define _PARAMETERH_

#include <string>
#include <vector>

class Parameter {
	

public:
	std::string name;
	std::vector<std::pair<double,double>> values;
	
	Parameter();
	void AddValue(const std::pair<double,double> &value); //looks up correct insert position and adds the value
	void AddValue(const double &moment,const double &value); //looks up correct insert position and adds the value

	void RemoveValue(size_t index);
	void SetValues(std::vector<std::pair<double,double>> values,bool sort=true);

	double GetValueAt(double time);
};

#endif