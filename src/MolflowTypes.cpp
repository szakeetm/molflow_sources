
#include "MolflowTypes.h"

/**
* \brief Overloaded += operator summing up all parts of a ProfileSlice: this += rhs
* \param rhs reference to a ProfileSlice right hand side
* \return reference of the lhs ProfileSlice
*/
ProfileSlice& ProfileSlice::operator+=(const ProfileSlice& rhs)
{
	this->countEquiv += rhs.countEquiv;
	this->sum_v_ort += rhs.sum_v_ort;
	this->sum_1_per_ort_velocity += rhs.sum_1_per_ort_velocity;
	return *this;
}

ProfileSlice operator+(const ProfileSlice& lhs,const ProfileSlice& rhs)
{
	ProfileSlice result(lhs);
	result += rhs;
	return result;
}

/**
* \brief Overloaded += operator summing up all parts of a TextureCell: this += rhs
* \param rhs reference to a TextureCell right hand side
* \return reference of the lhs TextureCell
*/
TextureCell& TextureCell::operator+=(const TextureCell& rhs)
{
	this->countEquiv += rhs.countEquiv;
	this->sum_v_ort_per_area += rhs.sum_v_ort_per_area;
	this->sum_1_per_ort_velocity += rhs.sum_1_per_ort_velocity;
	return *this;
}

/**
* \brief Overloaded + operator summing up all parts of a TextureCell: lhs + rhs
* \param rhs reference to a TextureCell right hand side
* \return reference of the lhs TextureCell
*/
TextureCell operator+(const TextureCell& lhs,const TextureCell& rhs)
{
	TextureCell result(lhs);
	result += rhs;
	return result;
}