//
// Created by pbahr on 29/10/2019.
//

#pragma once

#ifndef MOLFLOW_PROJ_OPTIXPOLYGON_H
#define MOLFLOW_PROJ_OPTIXPOLYGON_H

#include <cuda_runtime.h>
#include <limits>
#include <cstdint> // uint8_t

typedef unsigned long long int uint64_cu;

namespace flowgpu {

    struct Texel64;
    struct Texel32;

    // Texel typedef
#ifdef HIT64
        using Texel = Texel64;
#else
        using Texel = Texel32;
#endif

    enum RayType : uint8_t
    {
        RAY_TYPE_MOLECULE  = 0,
        RAY_TYPE_COUNT
    };

    enum FacetType : uint8_t
    {
        FACET_TYPE_SOLID  = 0,
#ifdef WITH_TRANS
        FACET_TYPE_TRANS  = 1,
#else
        FACET_TYPE_TRANS  = 0,
#endif // WITH_TRANS

        FACET_TYPE_COUNT
    };

    enum TEXTURE_FLAGS {
        noTexture = 0u,
        countAbs  = (1u << 0u),
        countRefl           = (1u << 1u),
        countTrans         = (1u << 2u),
        countDirection         = (1u << 3u)
        ,countDes = (1u << 4u) // Molflow only
    };

    enum PROFILE_FLAGS {
        noProfile = 0,
        profileU  = (1u << 0u),
        profileV           = (1u << 1u),
        // TODO: Not yet implemented
        profileAngular         = (1u << 2u),
        profileVelocity         = (1u << 3u)
        ,profileOrtVelocity = (1u << 4u)
        ,profileTanVelocity = (1u << 5u)
    };

    struct TempFacetProperties { //Formerly SHFACET
    public:
        //For sync between interface and subprocess
        double sticking;       // Sticking (0=>reflection  , 1=>absorption)   - can be overridden by time-dependent parameter
        double opacity;        // opacity  (0=>transparent , 1=>opaque)
        double area;           // Facet area (m^2)

        int    profileType;    // Profile type
        int    superIdx;       // Super structure index (Indexed from 0) -1: facet belongs to all structures (typically counter facets)
        size_t    superDest;      // Super structure destination index (Indexed from 1, 0=>current)
        int	 teleportDest;   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none, -1=>teleport to where it came from)

        bool   countAbs;       // Count absorption (MC texture)
        bool   countRefl;      // Count reflection (MC texture)
        bool   countTrans;     // Count transparent (MC texture)
        bool   countDirection;

        // Flags
        bool   is2sided;     // 2 sided
        bool   isProfile;    // Profile facet
        bool   isTextured;   // texture
        bool   isVolatile;   // Volatile facet (absorbtion facet which does not affect particule trajectory)

        // Geometry
        size_t nbIndex;   // Number of index/vertex
        //double sign;      // Facet vertex rotation (see Facet::DetectOrientation())

        // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
        float3   O;  // Origin
        float3   U;  // U vector
        float3   V;  // V vector
        float3   nU; // Normalized U
        float3   nV; // Normalized V

        // Normal vector
        float3    N;    // normalized
        float3    Nuv;  // normal to (u,v) not normlized

        // Hit/Abs/Des/Density recording on 2D texture map
        size_t    texWidth;    // Rounded texture resolution (U)
        size_t    texHeight;   // Rounded texture resolution (V)
        float texWidthD;   // Actual texture resolution (U)
        float texHeightD;  // Actual texture resolution (V)

        // ----
        // Molflow-specific facet parameters
        float temperature;    // Facet temperature (Kelvin) for velocity and outgassing calculation and soujourn time
        float outgassing;           // (in unit *m^3/s)  used to calculate the true facetOutgassing value

        /*int sticking_paramId;    // -1 if use constant value, 0 or more if referencing time-dependent parameter
        int opacity_paramId;     // -1 if use constant value, 0 or more if referencing time-dependent parameter
        int outgassing_paramId;  // -1 if use constant value, 0 or more if referencing time-dependent parameter

        int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
        int IDid;  //If time-dependent desorption, which is its ID*/

        int    desorbType;     // Desorption type
        float desorbTypeN;    // Exponent in Cos^N desorption type

        bool   countDes;       // Count desoprtion (MC texture)


        /*double maxSpeed;       // Max expected particle velocity (for velocity histogram)
        double accomodationFactor; // Thermal accomodation factor [0..1]
        bool   enableSojournTime;
        double sojournFreq, sojournE;*/

        // Facet hit counters
        // FacetHitBuffer tmpCounter; - removed as now it's time-dependent and part of the hits buffer

        // Moving facets
        //bool isMoving;

        //Outgassing map
        /*bool   useOutgassingFile;   //has desorption file for cell elements
        double outgassingFileRatio; //desorption file's sample/unit ratio
        size_t   outgassingMapWidth; //rounded up outgassing file map width
        size_t   outgassingMapHeight; //rounded up outgassing file map height*/

        float totalOutgassing; //total outgassing for the given facet

        //AnglemapParams anglemapParams;//Incident angle map
        // ----

    };

    struct Texel64 {
        Texel64() : countEquiv(0), sum_v_ort_per_area(0.0), sum_1_per_ort_velocity(0.0){}

        uint64_cu countEquiv;
        double sum_v_ort_per_area;
        double sum_1_per_ort_velocity;
    };

    struct Texel32 {
        Texel32() : countEquiv(0), sum_v_ort_per_area(0.0f), sum_1_per_ort_velocity(0.0f){}

        //TODO: 32bit counter for temporary
        uint32_t countEquiv;
        float sum_v_ort_per_area;
        float sum_1_per_ort_velocity;
    };

    struct FacetTexture {
        FacetTexture() : texelOffset(0), texWidth(0), texHeight(0), texWidth_precise(0.0f), texHeight_precise(0.0f), bbMin(), bbMax(){}

        unsigned int texelOffset;
        // Hit/Abs/Des/Density recording on 2D texture map
        unsigned int    texWidth;    // Rounded texture resolution (U)
        unsigned int    texHeight;   // Rounded texture resolution (V)
        float texWidth_precise;   // Actual texture resolution (U)
        float texHeight_precise;  // Actual texture resolution (V)

        float3 bbMin;
        float3 bbMax;

        //TODO: Maybe fit texCellInc with texels, needs smart offsets then with time-dependence added
        //float*   textureCellIncrements; // texWidth*texHeight
        //Texel*   texels; // texWidth*texHeight
    };

    /*struct FacetProfile {
        FacetProfile() : profileOffset(0), profileType(0){}

        unsigned int profileOffset;
        int profileType;
    };*/

    // TODO: Maybe save some variables (e.g. texture WxH) not with the Polygon to save memory when Facets are actually not textured
    struct TextureProperties{
        TextureProperties() : textureOffset(0), textureSize(0),textureFlags(TEXTURE_FLAGS::noTexture){}
        TextureProperties& operator=(const TextureProperties& o) = default;/*(const TextureProperties& o){
            this->textureOffset = o.textureOffset;
            this->textureSize = o.textureSize;
            this->textureFlags = o.textureFlags;

            return *this;
        }*/

        unsigned int textureOffset; // map it to the original polygon index
        unsigned int textureSize; // check ==0 to find if texture exists, could also get its own SBT
        unsigned int textureFlags; // enum TextureCounters

    };

    // TODO: Maybe save some variables (e.g. texture WxH) not with the Polygon to save memory when Facets are actually not textured
    struct ProfileProperties{
        ProfileProperties() : profileOffset(0), profileType(PROFILE_FLAGS::noProfile){}
        ProfileProperties& operator=(const ProfileProperties& o)=default;/*{
            this->profileOffset = o.profileOffset;
            //this->textureSize = o.textureSize;
            this->profileType = o.profileType;

            return *this;
        }*/

        unsigned int profileOffset; // map it to the original polygon index
        //unsigned int textureSize; // check ==0 to find if texture exists, could also get its own SBT
        unsigned int profileType; // enum TextureCounters

    };

    struct DesorpProperties{
        DesorpProperties() : desorbType(0u), cosineExponent(-1.0f){}
        DesorpProperties& operator=(const DesorpProperties& o)= default;

        uint8_t desorbType; // type of desorption
        float cosineExponent; // for Cosine^N
    };

    struct SimProperties{
        SimProperties() : stickingFactor(-1.0f), temperature(-1.0f), is2sided(false), opacity(-1.0f){}
        SimProperties& operator=(const SimProperties& o)=default;/*{
            this->stickingFactor = o.stickingFactor;
            this->temperature = o.temperature;
            this->is2sided = o.is2sided;
            this->opacity = o.opacity;
            return *this;
        }*/

        float stickingFactor;
        float temperature;

        // only for transparent / 2sided
        bool is2sided;
        float opacity;
    };

    class Polygon {
    public:
        Polygon()
        : nbVertices(0), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(),
        parentIndex(std::numeric_limits<unsigned int>::max()), texProps(), facProps(),desProps(){
        }
        Polygon(unsigned int nbOfVertices)
        : nbVertices(nbOfVertices), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(),
        parentIndex(std::numeric_limits<unsigned int>::max()), texProps(), facProps(), desProps(){
        }
        Polygon(Polygon&& o){
            *this = std::move(o);
        }
        Polygon(const Polygon& o){
            *this = o;
        }

        ~Polygon()=default;

        Polygon& operator=(Polygon&& o){
            if (this != &o)
            {
                this->nbVertices = o.nbVertices;
                this->indexOffset = o.indexOffset;
                this->O = o.O;
                this->U = o.U;
                this->V = o.V;
                this->Nuv = o.Nuv;

                this->Ox64 = o.Ox64;
                this->Ux64 = o.Ux64;
                this->Vx64 = o.Vx64;
                this->Nuvx64 = o.Nuvx64;

                this->nU = o.nU;
                this->nV = o.nV;
                this->N = o.N;
                this->nUx64 = o.nUx64;
                this->nVx64 = o.nVx64;
                this->Nx64 = o.Nx64;

                this->parentIndex = o.parentIndex;
                this->texProps = o.texProps;
                this->profProps = o.profProps;
                this->facProps = o.facProps;
                this->desProps = o.desProps;

                o.nbVertices = 0;
                o.indexOffset = 0;
                o.O = float3();
                o.U = float3();
                o.V = float3();
                o.Nuv = float3();

                o.Ox64 = double3();
                o.Ux64 = double3();
                o.Vx64 = double3();
                o.Nuvx64 = double3();

                o.nU = float3();
                o.nV = float3();
                o.N = float3();
                o.nUx64 = double3();
                o.nVx64 = double3();
                o.Nx64 = double3();

                o.parentIndex = std::numeric_limits<unsigned int>::max();
                o.texProps = TextureProperties();
                o.profProps = ProfileProperties();
                o.facProps = SimProperties();
                o.desProps = DesorpProperties();

            }
            return *this;
        }
        Polygon& operator=(const Polygon& o)=default;/*{

            this->nbVertices = o.nbVertices;
            this->indexOffset = o.indexOffset;
            this->O = o.O;
            this->U = o.U;
            this->V = o.V;
            this->Nuv = o.Nuv;
            this->N = o.N;
            this->nU = o.nU;
            this->nV = o.nV;

            this->parentIndex = o.parentIndex;
            this->texProps = o.texProps;
            this->profProps = o.profProps;
            this->facProps = o.facProps;
            this->desProps = o.desProps;

            return *this;
        }*/

        void copyParametersFrom(const Polygon& o){

            this->O = o.O;
            this->U = o.U;
            this->V = o.V;
            this->Nuv = o.Nuv;
            this->N = o.N;
            this->nU = o.nU;
            this->nV = o.nV;

            this->Ox64 = o.Ox64;
            this->Ux64 = o.Ux64;
            this->Vx64 = o.Vx64;
            this->Nuvx64 = o.Nuvx64;

            this->parentIndex = o.parentIndex;
            this->texProps = o.texProps;
            this->profProps = o.profProps;
            this->facProps = o.facProps;
            this->desProps = o.desProps;

        }

        // attributes that don't describe the geometry
        TextureProperties texProps;
        ProfileProperties profProps;
        SimProperties facProps;
        DesorpProperties desProps;

        unsigned int parentIndex; // map it to the original polygon index

        // variables for access to  global memory (indices, vertices)
        unsigned int nbVertices;
        unsigned int indexOffset;

        float3 center; // center for offset

        // variables for ray-plane (3d space) intersection
        float3 O;
        float3 U;
        float3 V;
        float3 Nuv;

        double3 Ox64;
        double3 Ux64;
        double3 Vx64;
        double3 Nuvx64;

        // normalized facet vectors
        float3 nU;
        float3 nV;
        float3 N;

        double3 nUx64;
        double3 nVx64;
        double3 Nx64;
    };
}

#endif //MOLFLOW_PROJ_OPTIXPOLYGON_H
