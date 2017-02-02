#include "Parameter.h"
#include <algorithm>
#include "GLApp/MathTools.h"

Parameter::Parameter() { //initialize with empty values
	name="";
	values=std::vector<std::pair<double,double>>();
}

void Parameter::AddValue(const std::pair<double,double> &value) {
	//Assuming existing values are stored in order
	size_t pos=0;
	for (;pos<this->values.size() && value.first>this->values[pos].first;pos++); //find pos to insert
	values.insert(this->values.begin()+pos,value); //inserts or appends, depending on pos
}

void Parameter::AddValue(const double &moment,const double &value){
	AddValue(std::make_pair(moment,value));
}

void Parameter::RemoveValue(size_t pos) {
	values.erase(values.begin()+pos);
}

struct myclass {
  bool operator() (const std::pair<double,double> &i,const std::pair<double,double> &j) {
	  return (i.first<j.first);
  }
} sorter;

void Parameter::SetValues(std::vector<std::pair<double,double>> insertValues,BOOL sort) { //sort then set
	if (sort) std::sort(insertValues.begin(),insertValues.end(),sorter); //sort pairs by time first
	this->values=insertValues;
}

double Parameter::GetValueAt(double time) {
	return InterpolateY(time,values,TRUE,FALSE); //linear interpolation, limited to bounds
}

