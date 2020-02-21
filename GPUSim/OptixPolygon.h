//
// Created by pbahr on 29/10/2019.
//

#pragma once

#ifndef MOLFLOW_PROJ_OPTIXPOLYGON_H
#define MOLFLOW_PROJ_OPTIXPOLYGON_H

#include <cereal/cereal.hpp>

namespace flowgeom {

    enum TextureCounters {
        noTexture = (1 << 0),
        countAbs  = (1 << 1),
        countRefl           = (1 << 2),
        countTrans         = (1 << 3),
        countDirection         = (1 << 4)
                ,countDes = (1 << 5) // Molflow only
    };

    struct FacetProperties { //Formerly SHFACET
    public:
        //For sync between interface and subprocess
        double sticking;       // Sticking (0=>reflection  , 1=>absorption)   - can be overridden by time-dependent parameter
        double opacity;        // opacity  (0=>transparent , 1=>opaque)
        double area;           // Facet area (m^2)

        int    profileType;    // Profile type
        int    superIdx;       // Super structure index (Indexed from 0) -1: facet belongs to all structures (typically counter facets)
        size_t    superDest;      // Super structure destination index (Indexed from 1, 0=>current)
        int	 teleportDest;   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none, -1=>teleport to where it came from)

        bool   countAbs;       // Count absoprtion (MC texture)
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
        double texWidthD;   // Actual texture resolution (U)
        double texHeightD;  // Actual texture resolution (V)

        // ----
        // Molflow-specific facet parameters
        double temperature;    // Facet temperature (Kelvin) for velocity and outgassing calculation and soujourn time
        double outgassing;           // (in unit *m^3/s)  used to calculate the true facetOutgassing value

        /*int sticking_paramId;    // -1 if use constant value, 0 or more if referencing time-dependent parameter
        int opacity_paramId;     // -1 if use constant value, 0 or more if referencing time-dependent parameter
        int outgassing_paramId;  // -1 if use constant value, 0 or more if referencing time-dependent parameter

        int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
        int IDid;  //If time-dependent desorption, which is its ID*/

        int    desorbType;     // Desorption type
        double desorbTypeN;    // Exponent in Cos^N desorption type

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

        double totalOutgassing; //total outgassing for the given facet

        //AnglemapParams anglemapParams;//Incident angle map
        // ----

        template<class Archive>
        void serialize(Archive & archive)
        {
            archive(
                    CEREAL_NVP(sticking),       // Sticking (0=>reflection  , 1=>absorption)   - can be overridden by time-dependent parameter
                    CEREAL_NVP(opacity),        // opacity  (0=>transparent , 1=>opaque)
                    CEREAL_NVP(area),          // Facet area (m^2)

                    CEREAL_NVP(profileType),    // Profile type
                    CEREAL_NVP(superIdx),       // Super structure index (Indexed from 0)
                    CEREAL_NVP(superDest),      // Super structure destination index (Indexed from 1, 0=>current)
                    CEREAL_NVP(teleportDest),   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none, -1=>teleport to where it came from)

                    CEREAL_NVP(countAbs),       // Count absoprtion (MC texture)
                    CEREAL_NVP(countRefl),      // Count reflection (MC texture)
                    CEREAL_NVP(countTrans),     // Count transparent (MC texture)
                    CEREAL_NVP(countDirection),

                    // Flags
                    CEREAL_NVP(is2sided),     // 2 sided
                    CEREAL_NVP(isProfile),    // Profile facet
                    CEREAL_NVP(isTextured),   // texture
                    CEREAL_NVP(isVolatile),   // Volatile facet (absorbtion facet which does not affect particule trajectory)

                    // Geometry
                    CEREAL_NVP(nbIndex),   // Number of index/vertex
                    //CEREAL_NVP(sign),      // Facet vertex rotation (see Facet::DetectOrientation())

                    // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
                    CEREAL_NVP(O),  // Origin
                    CEREAL_NVP(U),  // U vector
                    CEREAL_NVP(V),  // V vector
                    CEREAL_NVP(nU), // Normalized U
                    CEREAL_NVP(nV), // Normalized V

                    // Normal vector
                    CEREAL_NVP(N),    // normalized
                    CEREAL_NVP(Nuv),  // normal to (u,v) not normlized

                    // Hit/Abs/Des/Density recording on 2D texture map
                    CEREAL_NVP(texWidth),    // Rounded texture resolution (U)
                    CEREAL_NVP(texHeight),   // Rounded texture resolution (V)
                    CEREAL_NVP(texWidthD),   // Actual texture resolution (U)
                    CEREAL_NVP(texHeightD),  // Actual texture resolution (V)


                    // Molflow-specific facet parameters
                    CEREAL_NVP(temperature),    // Facet temperature (Kelvin)                  - can be overridden by time-dependent parameter
                    CEREAL_NVP(outgassing),           // (in unit *m^3/s)                      - can be overridden by time-dependent parameter

                    CEREAL_NVP(desorbType),     // Desorption type
                    CEREAL_NVP(desorbTypeN),    // Exponent in Cos^N desorption type

                    CEREAL_NVP(countDes),       // Count desoprtion (MC texture)

                    CEREAL_NVP(totalOutgassing) //total outgassing for the given facet
            );
        }
    };

    class Polygon {
    public:
        Polygon()
        : stickingFactor(-1.0), nbVertices(0), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(),
        parentIndex(std::numeric_limits<uint32_t>::max()), textureOffset(0), textureSize(0), textureFlags(TextureCounters::noTexture){
        }
        Polygon(int32_t nbOfVertices)
        : stickingFactor(-1.0), nbVertices(nbOfVertices), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(),
        parentIndex(std::numeric_limits<uint32_t>::max()), textureOffset(0), textureSize(0), textureFlags(TextureCounters::noTexture){
        }
        Polygon(Polygon&& o){
            *this = std::move(o);
        }
        Polygon(const Polygon& o){
            *this = o;
        }

        ~Polygon(){

        }

        Polygon& operator=(Polygon&& o){
            if (this != &o)
            {
                this->stickingFactor = o.stickingFactor;

                this->nbVertices = o.nbVertices;
                this->indexOffset = o.indexOffset;
                this->O = o.O;
                this->U = o.U;
                this->V = o.V;
                this->Nuv = o.Nuv;
                this->nU = o.nU;
                this->nV = o.nV;
                this->N = o.N;
                this->parentIndex = o.parentIndex;
                this->textureOffset = o.textureOffset;
                this->textureSize = o.textureSize;
                this->textureFlags = o.textureFlags;

                o.nbVertices = 0;
                o.indexOffset = 0;
                o.O = float3();
                o.U = float3();
                o.V = float3();
                o.Nuv = float3();
                o.nU = float3();
                o.nV = float3();
                o.N = float3();
                o.parentIndex = std::numeric_limits<uint32_t>::max();
                o.textureOffset = 0;
                o.textureSize = 0;
                o.textureFlags = TextureCounters::noTexture;

            }
            return *this;
        }
        Polygon& operator=(const Polygon& o){
            this->stickingFactor = o.stickingFactor;

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
            this->textureOffset = o.textureOffset;
            this->textureSize = o.textureSize;
            this->textureFlags = o.textureFlags;

            return *this;
        }

        void copyParametersFrom(const Polygon& o){
            this->stickingFactor = o.stickingFactor;

            this->O = o.O;
            this->U = o.U;
            this->V = o.V;
            this->Nuv = o.Nuv;
            this->N = o.N;
            this->nU = o.nU;
            this->nV = o.nV;

            this->parentIndex = o.parentIndex;
            this->textureOffset = o.textureOffset;
            this->textureSize = o.textureSize;
            this->textureFlags = o.textureFlags;

        }

        uint32_t parentIndex; // map it to the original polygon index
        uint32_t textureOffset; // map it to the original polygon index
        uint32_t textureSize; // check ==0 to find if texture exists, could also get its own SBT
        uint32_t textureFlags; // enum TextureCounters

        // attributes that don't describe the geometry
        float stickingFactor;

        // variables for access to  global memory (indices, vertices)
        uint32_t nbVertices;
        uint32_t indexOffset;

        // variables for ray-plane (3d space) intersection
        float3 O;
        float3 U;
        float3 V;
        float3 Nuv;

        // normalized facet vectors
        float3 nU;
        float3 nV;
        float3 N;

        template<class Archive>
        void serialize(Archive & archive)
        {
            archive(
                    stickingFactor,
                    nbVertices,
                    indexOffset,
                    O,
                    U,
                    V,
                    Nuv,
                    nU,
                    nV,
                    N
            );
        }
    };
}

#endif //MOLFLOW_PROJ_OPTIXPOLYGON_H
