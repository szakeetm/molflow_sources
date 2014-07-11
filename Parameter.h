#ifndef _PARAMETERH_
#define _PARAMETERH_

#include <string>
#include <vector>

class Parameter {
	Parameter();

public:
	std::string name;
	std::vector<std::pair<double,double>> values;

	void AddValue(const std::pair<double,double> &value); //looks up correct insert position and adds the value
	void RemoveValue(size_t index);
	void SetValues(std::vector<std::pair<double,double>> values);

	double GetValueAt(double time);
};

#endif