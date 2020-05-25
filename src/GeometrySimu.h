//
// Created by Pascal Baehr on 28.04.20.
//

#ifndef MOLFLOW_PROJ_GEOMETRYSIMU_H
#define MOLFLOW_PROJ_GEOMETRYSIMU_H

#include <vector>
#include "MolflowTypes.h"
#include "Buffer_shared.h"

class Anglemap {
public:
    std::vector<size_t>   pdf;		  // Incident angle distribution, phi and theta, not normalized. Used either for recording or for 2nd order interpolation
    std::vector<double>   phi_CDFs;    // A table containing phi distributions for each theta, starting from 0 for every line (1 line = 1 theta value). For speed we keep it in one memory block, 1 pointer
    std::vector<size_t>   phi_CDFsums; // since CDF runs only to the middle of the last segment, for each theta a line sum is stored here. Also a pdf for theta
    std::vector<double>   theta_CDF;	  // Theta CDF, not normalized. nth value is the CDF at the end of region n (beginning of first section is always 0)
    size_t   theta_CDFsum; // since theta CDF only runs till the middle of the last segment, the map sum is here
};

// Local facet structure
struct SubprocessFacet{
    FacetProperties sh;

    std::vector<size_t>      indices;          // Indices (Reference to geometry vertex)
    std::vector<Vector2d> vertices2;        // Vertices (2D plane space, UV coordinates)
    std::vector<std::vector<TextureCell>>     texture;            // Texture hit recording (taking area, temperature, mass into account), 1+nbMoments
    std::vector<double>   textureCellIncrements;              // Texure increment
    std::vector<bool>     largeEnough;      // cells that are NOT too small for autoscaling
    double   fullSizeInc;       // Texture increment of a full texture element
    std::vector<std::vector<DirectionCell>>     direction;       // Direction field recording (average), 1+nbMoments
    //bool     *fullElem;         // Direction field recording (only on full element)
    std::vector<std::vector<ProfileSlice>> profile;         // Distribution and hit recording
    std::vector<double>   outgassingMap; // Cumulative outgassing map when desorption is based on imported file
    double outgassingMapWidthD; //actual outgassing file map width
    double outgassingMapHeightD; //actual outgassing file map height
    Anglemap angleMap;

    // Temporary var (used in Intersect for collision)
    double colDist;
    double colU;
    double colV;
    double rw;
    double rh;
    double iw;
    double ih;

    // Temporary var (used in FillHit for hit recording)
    bool   isHit;
    bool   isReady;         // Volatile state
    size_t    textureSize;   // Texture size (in bytes)
    size_t    profileSize;   // profile size (in bytes)
    size_t    directionSize; // direction field size (in bytes)
    size_t    angleMapSize;  // incidentangle map size (in bytes)

    /*int CDFid; //Which probability distribution it belongs to (one CDF per temperature)
    int IDid;  //If time-dependent desorption, which is its ID*/
    size_t globalId; //Global index (to identify when superstructures are present)

    // Facet hit counters
    std::vector<FacetHitBuffer> tmpCounter; //1+nbMoment
    std::vector<FacetHistogramBuffer> tmpHistograms; //1+nbMoment

    void  ResetCounter();
    void	ResizeCounter(size_t nbMoments);
    bool InitializeOnLoad(const size_t &id, const size_t &nbMoments, size_t &histogramTotalSize);

    void InitializeHistogram(const size_t &nbMoments, size_t &histogramTotalSize);

    bool InitializeDirectionTexture(const size_t &nbMoments);

    bool InitializeProfile(const size_t &nbMoments);

    bool InitializeTexture(const size_t &nbMoments);

    bool InitializeAngleMap();

    void InitializeOutgassingMap();

    bool InitializeLinkAndVolatile(const size_t & id);

    //void RegisterTransparentPass(SubprocessFacet *facet); //Allows one shared Intersect routine between MolFlow and Synrad
};

// Local simulation structure

class AABBNODE;

class SuperStructure {
public:
    SuperStructure();
    ~SuperStructure();
    std::vector<SubprocessFacet>  facets;   // Facet handles
    AABBNODE* aabbTree; // Structure AABB tree
};

#endif //MOLFLOW_PROJ_GEOMETRYSIMU_H
