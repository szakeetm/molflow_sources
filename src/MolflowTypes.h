
#pragma once
#include "GeometryTypes.h"
//#include "Parameter.h"
//#include "Interface.h"
#include <stddef.h> //size_t for gcc
#include <string>
#include <vector>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <Distributions.h>
//#include "GeometryTypes.h" // SelectionGroup

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
#define HIT_SCATTER 8
#define HIT_VOLUME_DECAY 9
#define HIT_LAST 10

#define MC_MODE 0         // Monte Carlo simulation mode
#define AC_MODE 1         // Angular coefficient simulation mode

#define MBARLS_TO_PAM3S 0.1 //Multiply by this to convert a value expressed in mbar.l/s to Pa.m3/s
#define PAM3S_TO_MBARLS 10.0 //Multiply by this to convert a value expressed in mbar.l/s to Pa.m3/s
#define MBAR_TO_PA 100 //Multiply by this to convert a value expressed in mbar to Pa

typedef float ACFLOAT;
struct UserMoment {
	std::string content;
	double timeWindow;
};
struct Moment {
	double time;
	double window;
	/*
	inline double GetElement(const bool first) const {
		return first ? time : window;
	}
	*/
};
struct Interval {
	double startTime;
	double endTime;
};
struct DesorptionEntry {
	double time;
	double desValue;
	DesorptionEntry(double _time, double _desVal) :time(_time), desValue(_desVal) {};
	inline double GetElement(const bool first) const {
		return first ? time : desValue;
	}
};
struct IntegratedDesorptionEntry {
	double time;
	double cumulativeDesValue;
	IntegratedDesorptionEntry(double _time, double _cumDesVal) :time(_time), cumulativeDesValue(_cumDesVal) {};
	inline double GetElement(const bool first) const {
		return first ? time : cumulativeDesValue;
	}
};

struct IntegratedVelocityEntry {
	double velocity; // m/s
	double cumulativeProbability; //monotonous increase from 0 to 1
	inline double GetElement(const bool first) const {
		return first ? velocity : cumulativeProbability;
	}
};


struct OutgassingMap {
	OutgassingMap() = default;
	size_t   outgassingMapWidth = 0; //rounded up outgassing file map width
	size_t   outgassingMapHeight = 0; //rounded up outgassing file map height
	double outgassingMapWidth_precise = 0.0; //actual outgassing file map width
	double outgassingMapHeight_precise = 0.0; //actual outgassing file map height
	double outgassingFileRatioU = 0.0; //desorption file's sample/unit ratio in U direction
	double outgassingFileRatioV = 0.0; //desorption file's sample/unit ratio in V direction
	std::vector<double>   outgassingMap_cdf; // Cumulative outgassing map when desorption is based on imported file
	std::vector<double>   outgassingMap; // Cumulative outgassing map when desorption is based on imported file

	// Analytic properties
	double totalDose = 0.0;
	double totalFlux = 0.0;
};

// Density/Hit field stuff

class ProfileSlice {
public:
	double countEquiv = 0.0;
	double sum_v_ort = 0.0;
	double sum_1_per_ort_velocity = 0.0;
	ProfileSlice& operator+=(const ProfileSlice& rhs);

	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(
			countEquiv,
			sum_v_ort,
			sum_1_per_ort_velocity
		);
	}
};
ProfileSlice operator+(const ProfileSlice& lhs, const ProfileSlice& rhs);

class TextureCell {
public:
	double countEquiv = 0.0;
	double sum_v_ort_per_area = 0.0;
	double sum_1_per_ort_velocity = 0.0;
	TextureCell& operator+=(const TextureCell& rhs);
	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(
			countEquiv,
			sum_v_ort_per_area,
			sum_1_per_ort_velocity
		);
	}
};
TextureCell operator+(const TextureCell& lhs, const TextureCell& rhs); //non-member function so std::plus can use it

//Texture limit types
struct TEXTURE_MOMENT_TYPE {
	double steady_state;
	double moments_only;
	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(
			CEREAL_NVP(steady_state),
			CEREAL_NVP(moments_only)
		);
	}
};

struct TEXTURE_MIN_MAX {
	TEXTURE_MOMENT_TYPE min;
	TEXTURE_MOMENT_TYPE max;
	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(
			CEREAL_NVP(min),
			CEREAL_NVP(max)
		);
	}
};

struct TEXTURE_SCALE_TYPE {
	TEXTURE_MIN_MAX manual;
	TEXTURE_MIN_MAX autoscale;
};


enum AutoScaleMode : int {
	AutoscaleMomentsOnly,
	AutoscaleMomentsAndConstFlow,
	AutoscaleConstFlow
};

class AnglemapParams {
public:

	bool   record = false; // Record incident angle 2-dim distribution
	size_t phiWidth = 0; //resolution between -PI and +PI
	double thetaLimit = 1.570796326; //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2. Slightly lower than PI/2 def value
	size_t thetaLowerRes = 0; //resolution between 0 and angleMapThetaLimit
	size_t thetaHigherRes = 0; //resolution between angleMapThetaLimit and PI/2

	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(
			record, // Record incident angle 2-dim distribution
			phiWidth, //resolution between -PI and +PI
			thetaLimit, //angle map can have a different resolution under and over the limit. Must be between 0 and PI/2
			thetaLowerRes, //resolution between 0 and angleMapThetaLimit
			thetaHigherRes //resolution between angleMapThetaLimit and PI/2
		);
	}

	size_t GetMapSize() const {
		return phiWidth * (thetaLowerRes + thetaHigherRes);
	}
	size_t GetRecordedMapSize() const {
		/*if (!hasRecorded) return 0;
		else */return GetMapSize();
	}
	size_t GetDataSize() const {
		return sizeof(size_t) * GetMapSize();
	}
	size_t GetRecordedDataSize() const {
		return sizeof(size_t) * GetRecordedMapSize();
	}
};

class ReflectionParam {
public:
	double diffusePart = 1.0;
	double specularPart = 0.0;
	double cosineExponent = 0.0; //Cos^N part: 1-diffuse-specular

	template<class Archive>
	void serialize(Archive& archive)
	{
		archive(diffusePart, specularPart, cosineExponent);
	}
};

struct MolflowInterfaceSettings : InterfaceSettings {
	std::vector<UserMoment> userMoments;
	TEXTURE_SCALE_TYPE textureLimits[3]; //auto- and manual scaling limits for textures
	AutoScaleMode texAutoscaleMode = AutoscaleMomentsAndConstFlow;
};