/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY / Pascal BAEHR
Copyright:   E.S.R.F / CERN
Website:     https://cern.ch/molflow

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/
#pragma once
#include "GLApp/GLTypes.h"
#include "Parameter.h"
#include <stddef.h> //size_t for gcc
#include <string>
#include <vector>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <Distributions.h>
#include "GeometryTypes.h" // SelectionGroup

//#include "Buffer_shared.h"

// Desorption type
#define DES_NONE    0   // No desorption
#define DES_UNIFORM 1   // Uniform
#define DES_COSINE  2   // cos(theta)
#define DES_COSINE_N 3 // cos(theta)^N
#define DES_ANGLEMAP 4 //imported file

// (Old) Reflection types
#define REFLECTION_DIFFUSE 0   // Diffuse (cosine law)
#define REFLECTION_SPECULAR  1   // Mirror
#define REFLECTION_UNIFORM 2   // Uniform (for testing)

// Profile type
#define PROFILE_NONE       0  // No recording
#define PROFILE_U  1  // Pressure and density profile (U direction)
#define PROFILE_V  2  // Pressure and density profile (V direction)
#define PROFILE_ANGULAR    3  // Angular profile
#define PROFILE_VELOCITY   4 //Velocity distribution
#define PROFILE_ORT_VELOCITY 5 //Orthogonal velocity component
#define PROFILE_TAN_VELOCITY 6 //Tangential velocity (experimental)

// Hit type
#define HIT_DES   1
#define HIT_ABS   2
#define HIT_REF   3
#define HIT_TRANS 4
#define HIT_TELEPORTSOURCE 5
#define HIT_TELEPORTDEST 6
#define HIT_MOVING 7
#define HIT_LAST 10

#define MC_MODE 0         // Monte Carlo simulation mode
#define AC_MODE 1         // Angular coefficient simulation mode

#define MBARLS_TO_PAM3S 0.1 //Multiply by this to convert a value expressed in mbar.l/s to Pa.m3/s
#define MBAR_TO_PA 100 //Multiply by this to convert a value expressed in mbar to Pa

typedef float ACFLOAT;
typedef std::pair<std::string,double> UserMoment;
typedef std::pair<double,double> Moment; //mid-time,window_size
typedef std::pair<double,double> ID_p;
typedef std::pair<double,double> CDF_p;

struct UserInput {
    std::vector<SelectionGroup> selections;
    std::vector<std::tuple<bool, bool>> facetViewSettings;
    std::vector<UserMoment> userMoments;
    std::vector<Parameter> parameters;
};

class IntegratedDesorption : public Distribution2D {

};

struct OutgassingMap {
    OutgassingMap() = default;
    size_t   outgassingMapWidth; //rounded up outgassing file map width
    size_t   outgassingMapHeight; //rounded up outgassing file map height
    double outgassingMapWidth_precise; //actual outgassing file map width
    double outgassingMapHeight_precise; //actual outgassing file map height
    double outgassingFileRatioU; //desorption file's sample/unit ratio in U direction
    double outgassingFileRatioV; //desorption file's sample/unit ratio in V direction
    std::vector<double>   outgassingMap_cdf; // Cumulative outgassing map when desorption is based on imported file
    std::vector<double>   outgassingMap; // Cumulative outgassing map when desorption is based on imported file

    // Analytic properties
    double totalDose;
    double totalFlux;
};

// Density/Hit field stuff
#define HITMAX 1E38
class ProfileSlice {
public:
	double countEquiv=0.0;
	double sum_v_ort=0.0;
	double sum_1_per_ort_velocity=0.0;
	ProfileSlice& operator+=(const ProfileSlice& rhs);
	ProfileSlice& operator+(const ProfileSlice& rhs);
    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                countEquiv,
                sum_v_ort,
                sum_1_per_ort_velocity
        );
    }
};

class TextureCell {
public:
	double countEquiv=0.0;
	double sum_v_ort_per_area=0.0;
	double sum_1_per_ort_velocity=0.0;
	TextureCell& operator+=(const TextureCell& rhs);
	TextureCell& operator+(const TextureCell& rhs);
    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(
                countEquiv,
                sum_v_ort_per_area,
                sum_1_per_ort_velocity
        );
    }
};

//Texture limit types
typedef struct {
	double steady_state;
	double moments_only;
} TEXTURE_MOMENT_TYPE;

typedef struct {
	TEXTURE_MOMENT_TYPE min;
	TEXTURE_MOMENT_TYPE max;
} TEXTURE_MIN_MAX;

typedef struct {
	TEXTURE_MIN_MAX manual;
	TEXTURE_MIN_MAX autoscale;
} TEXTURE_SCALE_TYPE;

class AnglemapParams {
public:
    AnglemapParams(){
        record = false;
        phiWidth = 0;
        thetaLimit = 0.0;
        thetaLimit = 0;
		thetaLowerRes = 0;
        thetaHigherRes = 0;
    }
	bool   record; // Record incident angle 2-dim distribution
	//bool hasRecorded;
	size_t phiWidth; //resolution between -PI and +PI
	double thetaLimit; //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2
	size_t thetaLowerRes; //resolution between 0 and angleMapThetaLimit
	size_t thetaHigherRes; //resolution between angleMapThetaLimit and PI/2
	
	template<class Archive>
	void serialize(Archive & archive)
	{
		archive(
			   record, // Record incident angle 2-dim distribution
		 phiWidth, //resolution between -PI and +PI
		 thetaLimit, //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2
		 thetaLowerRes, //resolution between 0 and angleMapThetaLimit
		 thetaHigherRes //resolution between angleMapThetaLimit and PI/2
		);
	}

	size_t GetMapSize() const{
		return phiWidth * (thetaLowerRes + thetaHigherRes);
	}
	size_t GetRecordedMapSize() const{
		/*if (!hasRecorded) return 0;
		else */return GetMapSize();
	}
	size_t GetDataSize() const {
		return sizeof(size_t)*GetMapSize();
	}
	size_t GetRecordedDataSize() const {
		return sizeof(size_t)*GetRecordedMapSize();
	}
};

class Reflection {
public:
	double diffusePart;
	double specularPart;
	double cosineExponent; //Cos^N part: 1-diffuse-specular
	
	template<class Archive>
	void serialize(Archive & archive)
	{
		archive(diffusePart, specularPart, cosineExponent);
	}
};

//Just for AC matrix calculation in Molflow, old mesh structure:
typedef struct {

	float   area;     // Area of element
	float   uCenter;  // Center coordinates
	float   vCenter;  // Center coordinates
	int     elemId;   // Element index (MESH array)
	bool    full;     // Element is full

} SHELEM_OLD;

