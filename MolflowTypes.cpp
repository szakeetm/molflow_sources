#include "MolflowTypes.h"

ProfileSlice& ProfileSlice::operator+=(const ProfileSlice& rhs)
{
	this->countEquiv += rhs.countEquiv;
	this->sum_v_ort += rhs.sum_v_ort;
	this->sum_1_per_ort_velocity += rhs.sum_1_per_ort_velocity;
	return *this;
}

TextureCell& TextureCell::operator+=(const TextureCell& rhs)
{
	this->countEquiv += rhs.countEquiv;
	this->sum_v_ort_per_area += rhs.sum_v_ort_per_area;
	this->sum_1_per_ort_velocity += rhs.sum_1_per_ort_velocity;
	return *this;
}
