//
// Created by Pascal Baehr on 28.04.20.
//

#include <sstream>
#include <Helper/MathTools.h>
#include <cmath>
#include <set>
#include <Simulation/CDFGeneration.h>
#include <Simulation/IDGeneration.h>
#include <Helper/Chronometer.h>
#include <Helper/ConsoleLogger.h>
#include <Polygon.h>
#include "GeometrySimu.h"
#include "IntersectAABB_shared.h" // include needed for recursive delete of AABBNODE

SuperStructure::~SuperStructure()
{
    aabbTree.reset();
}

bool ParameterSurface::IsHardHit(const Ray &r) {
    const double td_opacity = dist->InterpolateY(r.time, false);
    if(td_opacity >= 1.0)
        return true;
    else
        return (r.rng->rnd() < td_opacity);
};

bool SubprocessFacet::InitializeOnLoad(const size_t &id, const size_t &nbMoments) {
    globalId = id;
    //ResizeCounter(nbMoments); //Initialize counter
    if (!InitializeLinkAndVolatile(id)) return false;
    InitializeOutgassingMap();
    InitializeAngleMap();
    InitializeTexture(nbMoments);
    InitializeProfile(nbMoments);
    InitializeDirectionTexture(nbMoments);
    InitializeHistogram(nbMoments);

    return true;
}

size_t SubprocessFacet::InitializeHistogram(const size_t &nbMoments) const
{
    //FacetHistogramBuffer hist;
    //hist.Resize(sh.facetHistogramParams);
    //tmpHistograms = std::vector<FacetHistogramBuffer>(1 + nbMoments, hist);
    size_t histSize = (1 + nbMoments) *
                                   (sh.facetHistogramParams.GetBouncesDataSize()
                                    + sh.facetHistogramParams.GetDistanceDataSize()
                                    + sh.facetHistogramParams.GetTimeDataSize());

    return histSize;
}

size_t SubprocessFacet::InitializeDirectionTexture(const size_t &nbMoments)
{
    size_t directionSize = 0;

        //Direction
    if (sh.countDirection) {
        directionSize = sh.texWidth*sh.texHeight * sizeof(DirectionCell);
        try {
            //direction = std::vector<std::vector<DirectionCell>>(1 + nbMoments, std::vector<DirectionCell>(sh.texWidth*sh.texHeight));
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load direction textures");
        }
    }
    else
        directionSize = 0;
    return directionSize;
}

size_t SubprocessFacet::InitializeProfile(const size_t &nbMoments)
{    size_t profileSize = 0;

    //Profiles
    if (sh.isProfile) {
        profileSize = PROFILE_SIZE * sizeof(ProfileSlice);
        try {
            //profile = std::vector<std::vector<ProfileSlice>>(1 + nbMoments, std::vector<ProfileSlice>(PROFILE_SIZE));
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load profiles");
            return false;
        }
    }
    else
        profileSize = 0;
    return profileSize;
}

std::vector<double> SubprocessFacet::InitTextureMesh()
{
    GLAppPolygon P1, P2;
    double sx, sy;
    double iw = 1.0 / (double)sh.texWidth_precise;
    double ih = 1.0 / (double)sh.texHeight_precise;
    double rw = sh.U.Norme() * iw;
    double rh = sh.V.Norme() * ih;
    double fullCellArea = iw*ih;

    std::vector<Vector2d>(4).swap(P1.pts);
    //P1.sign = 1;
    P2.pts = vertices2;
    //P2.sign = -sign;

    std::vector<double> interCellArea;
    interCellArea.resize(sh.texWidth * sh.texHeight, -1.0); //will shrink at the end

    for (size_t j = 0;j < sh.texHeight;j++) {
        sy = (double)j;
        for (size_t i = 0;i < sh.texWidth;i++) {
            sx = (double)i;

            bool allInside = false;
            double u0 = sx * iw;
            double v0 = sy * ih;
            double u1 = (sx + 1.0) * iw;
            double v1 = (sy + 1.0) * ih;
            //mesh[i + j*wp.texWidth].elemId = -1;

            if (sh.nbIndex <= 4) {

                // Optimization for quad and triangle
                allInside = IsInPoly(Vector2d(u0,v0), vertices2)
                            && IsInPoly(Vector2d(u0, v1), vertices2)
                            && IsInPoly(Vector2d(u1, v0), vertices2)
                            && IsInPoly(Vector2d(u1, v1), vertices2);

            }

            if (!allInside) {
                double area{};

                // Intersect element with the facet (facet boundaries)
                P1.pts[0] = {u0, v0};
                P1.pts[1] = {u1, v0};
                P1.pts[2] = {u1, v1};
                P1.pts[3] = {u0, v1};

                std::vector<bool>visible(P2.pts.size());
                std::fill(visible.begin(), visible.end(), true);
                auto [A,center,vList] = GetInterArea(P1, P2, visible);
                if (!IsZero(A)) {

                    if (A > (fullCellArea + 1e-10)) {

                        // Polyon intersection error !
                        // Switch back to brute force
                        auto [bfArea, center_loc] = GetInterAreaBF(P2, Vector2d(u0, v0), Vector2d(u1, v1));
                        bool fullElem = IsZero(fullCellArea - bfArea);
                        if (!fullElem) {
                            interCellArea[i + j*sh.texWidth] = (bfArea*(rw*rh) / (iw*ih));
                        }
                        else {
                            interCellArea[i + j*sh.texWidth] = -1.0;
                        }

                        //cellprop.full = IsZero(fullCellArea - A);

                    }
                    else {

                        bool fullElem = IsZero(fullCellArea - A);
                        if (!fullElem) {
                            // !! P1 and P2 are in u,v coordinates !!
                            interCellArea[i + j*sh.texWidth] = (A*(rw*rh) / (iw*ih));
                        }
                        else {
                            interCellArea[i + j*sh.texWidth] = -1.0;
                        }

                    }

                }
                else interCellArea[i + j*sh.texWidth] = -2.0; //zero element

            }
            else {  //All indide and triangle or quad
                interCellArea[i + j*sh.texWidth] = -1.0;
            }
        }
    }

    return interCellArea;
}

size_t SubprocessFacet::InitializeTexture(const size_t &nbMoments)
{
    size_t textureSize = 0;
    //Textures
    if (sh.isTextured) {
        size_t nbE = sh.texWidth*sh.texHeight;
        largeEnough.resize(nbE);
        textureSize = nbE * sizeof(TextureCell);
        /*try {
            texture = std::vector<std::vector<TextureCell>>(1 + nbMoments, std::vector<TextureCell>(nbE));
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load textures");
            return false;
        }
        fullSizeInc = (sh.texWidth_precise * sh.texHeight_precise) / (sh.U.Norme() * sh.V.Norme());
        }*/
        // Texture increment of a full texture element
        double fullSizeInc = (sh.texWidth_precise * sh.texHeight_precise) / (sh.U.Norme() * sh.V.Norme());
        for (size_t j = 0; j < nbE; j++) { //second pass, filter out very small cells
            largeEnough[j] = textureCellIncrements[j] < (5.0*fullSizeInc);
        }

        //double iw = 1.0 / (double)sh.texWidth_precise;
        //double ih = 1.0 / (double)sh.texHeight_precise;
        //double rw = sh.U.Norme() * iw;
        //double rh = sh.V.Norme() * ih;
    }
    else
        textureSize = 0;
    return textureSize;
}


size_t SubprocessFacet::InitializeAngleMap()
{
    //Incident angle map
    size_t angleMapSize = 0;
    if (sh.desorbType == DES_ANGLEMAP) { //Use mode
        //if (angleMapCache.empty()) throw Error(("Facet " + std::to_string(globalId + 1) + ": should generate by angle map but has none recorded.").c_str());

        //Construct CDFs
        try {
            angleMap.phi_CDFsums.resize(sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (phi CDF line sums)");
            return false;
        }
        try {
            angleMap.theta_CDF.resize(sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (line sums, CDF)");
            return false;
        }
        try {
            angleMap.phi_CDFs.resize(sh.anglemapParams.phiWidth * (sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes));
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (CDF)");
            return false;
        }

        //First pass: determine sums
        angleMap.theta_CDFsum = 0;
        memset(angleMap.phi_CDFsums.data(), 0, sizeof(size_t) * (sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes));
        for (size_t thetaIndex = 0; thetaIndex < (sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes); thetaIndex++) {
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                angleMap.phi_CDFsums[thetaIndex] += angleMap.pdf[thetaIndex*sh.anglemapParams.phiWidth + phiIndex];
            }
            angleMap.theta_CDFsum += angleMap.phi_CDFsums[thetaIndex];
        }
        if (!angleMap.theta_CDFsum) {
            auto err = fmt::format("Facet {} has all-zero recorded angle map.", globalId + 1);
            throw std::runtime_error(err.c_str());
            return false;
        }

        //Second pass: write CDFs
        double thetaNormalizingFactor = 1.0 / (double)angleMap.theta_CDFsum;
        for (size_t thetaIndex = 0; thetaIndex < (sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes); thetaIndex++) {
            if (angleMap.theta_CDFsum == 0) { //no hits in this line, generate CDF of uniform distr.
                angleMap.theta_CDF[thetaIndex] = (0.5 + (double)thetaIndex) / (double)(sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes);
            }
            else {
                if (thetaIndex == 0) {
                    //First CDF value, covers half of first segment
                    angleMap.theta_CDF[thetaIndex] = 0.5 * (double)angleMap.phi_CDFsums[0] * thetaNormalizingFactor;
                }
                else {
                    //value covering second half of last segment and first of current segment
                    angleMap.theta_CDF[thetaIndex] = angleMap.theta_CDF[thetaIndex - 1] + (double)(angleMap.phi_CDFsums[thetaIndex - 1] + angleMap.phi_CDFsums[thetaIndex])*0.5*thetaNormalizingFactor;
                }
            }
            double phiNormalizingFactor = 1.0 / (double)angleMap.phi_CDFsums[thetaIndex];
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                size_t index = sh.anglemapParams.phiWidth * thetaIndex + phiIndex;
                if (angleMap.phi_CDFsums[thetaIndex] == 0) { //no hits in this line, create CDF of uniform distr.
                    angleMap.phi_CDFs[index] = (0.5 + (double)phiIndex) / (double)sh.anglemapParams.phiWidth;
                }
                else {
                    if (phiIndex == 0) {
                        //First CDF value, covers half of first segment
                        angleMap.phi_CDFs[index] = 0.5 * (double)angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex] * phiNormalizingFactor;
                    }
                    else {
                        //value covering second half of last segment and first of current segment
                        angleMap.phi_CDFs[index] = angleMap.phi_CDFs[sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + (double)(angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex])*0.5*phiNormalizingFactor;
                    }
                }
            }
        }
    }
    else {
        //Record mode, create pdf vector
        angleMap.pdf.resize(sh.anglemapParams.GetMapSize());
    }

    if(sh.anglemapParams.record)
        angleMapSize += sh.anglemapParams.GetDataSize();
    return angleMapSize;
}

void SubprocessFacet::InitializeOutgassingMap()
{
    if (sh.useOutgassingFile) {
        //Precalc actual outgassing map width and height for faster generation:
        ogMap.outgassingMapWidth_precise = sh.U.Norme() * ogMap.outgassingFileRatioU;
        ogMap.outgassingMapHeight_precise = sh.V.Norme() * ogMap.outgassingFileRatioV;
        size_t nbE = ogMap.outgassingMapWidth*ogMap.outgassingMapHeight;
        // TODO: Check with molflow_threaded e10c2a6f and 66b89ac7 if right
        // making a copy shouldn't be necessary as i will never get changed before use
        //outgassingMapWindow = facetRef->outgassingMapWindow; //init by copying pdf
        ogMap.outgassingMap_cdf = ogMap.outgassingMap;
        for (size_t i = 1; i < nbE; i++) {
            ogMap.outgassingMap_cdf[i] = ogMap.outgassingMap_cdf[i - 1] + ogMap.outgassingMap_cdf[i]; //Convert p.d.f to cumulative distr.
        }
    }
}

bool SubprocessFacet::InitializeLinkAndVolatile(const size_t & id)
{
    if (sh.superDest || sh.isVolatile) {
        // Link or volatile facet, overides facet settings
        // Must be full opaque and 0 sticking
        // (see SimulationMC.c::PerformBounce)
        //sh.isOpaque = true;
        sh.opacity = 1.0;
        sh.opacity_paramId = -1;
        sh.sticking = 0.0;
        sh.sticking_paramId = -1;
    }
    return true;
}

/*void SubprocessFacet::ResetCounter() {
    std::fill(tmpCounter.begin(), tmpCounter.end(), FacetHitBuffer());
}

void SubprocessFacet::ResizeCounter(size_t nbMoments) {
    tmpCounter = std::vector<FacetHitBuffer>(nbMoments + 1); //Includes 0-init
    tmpCounter.shrink_to_fit();
    tmpHistograms = std::vector<FacetHistogramBuffer>(nbMoments + 1);
    tmpHistograms.shrink_to_fit();
}*/

/**
* \brief Constructor for cereal initialization
*/
SubprocessFacet::SubprocessFacet() : Facet() {
    isReady = false;
    globalId = 0;
}

/**
* \brief Constructor with initialisation based on the number of indices/facets
* \param nbIndex number of indices/facets
*/
SubprocessFacet::SubprocessFacet(size_t nbIndex) : Facet(nbIndex) {
    isReady = false;
    globalId = 0;
    indices.resize(nbIndex);                    // Ref to Geometry Vector3d
    vertices2.resize(nbIndex);
}

SubprocessFacet::SubprocessFacet(const SubprocessFacet& cpy)  : Facet(cpy) {
    *this = cpy;
}

SubprocessFacet::SubprocessFacet(SubprocessFacet&& cpy) noexcept : Facet(cpy){
    *this = std::move(cpy);
}

SubprocessFacet& SubprocessFacet::operator=(const SubprocessFacet& cpy){
    this->angleMap = cpy.angleMap;
    this->ogMap = cpy.ogMap;
    this->largeEnough = cpy.largeEnough;
    this->textureCellIncrements = cpy.textureCellIncrements;
    this->sh = cpy.sh;

    isReady = cpy.isReady;
    globalId = cpy.globalId;
    indices = cpy.indices;                    // Ref to Geometry Vector3d
    vertices2 = cpy.vertices2;
    if(cpy.surf) surf = cpy.surf;
    else surf = nullptr;

    return *this;
}

SubprocessFacet& SubprocessFacet::operator=(SubprocessFacet&& cpy) noexcept {
    this->angleMap = cpy.angleMap;
    this->ogMap = cpy.ogMap;
    this->largeEnough = std::move(cpy.largeEnough);
    this->textureCellIncrements = std::move(cpy.textureCellIncrements);
    this->sh = cpy.sh;

    isReady = cpy.isReady;
    globalId = cpy.globalId;
    indices = std::move(cpy.indices);                    // Ref to Geometry Vector3d
    vertices2 = std::move(cpy.vertices2);
    surf = cpy.surf;
    cpy.surf = nullptr;

    return *this;
}

/**
* \brief Calculates the hits size for a single facet which is necessary for hits dataport
* \param nbMoments amount of moments
* \return calculated size of the facet hits
*/
size_t SubprocessFacet::GetHitsSize(size_t nbMoments) const { //for hits dataport
    return   (1 + nbMoments)*(
            sizeof(FacetHitBuffer) +
            +(sh.isTextured ? (sh.texWidth*sh.texHeight * sizeof(TextureCell)) : 0)
            + (sh.isProfile ? (PROFILE_SIZE * sizeof(ProfileSlice)) : 0)
            + (sh.countDirection ? (sh.texWidth*sh.texHeight * sizeof(DirectionCell)) : 0)
            + sh.facetHistogramParams.GetDataSize()
    ) + (sh.anglemapParams.record ? (sh.anglemapParams.GetRecordedDataSize()) : 0);

}

size_t SubprocessFacet::GetMemSize() const {
    size_t sum = 0;
    sum += sizeof (SubprocessFacet);
    sum += sizeof (size_t) * indices.capacity();
    sum += sizeof (Vector2d) * vertices2.capacity();
    sum += sizeof (double) * textureCellIncrements.capacity();
    sum += sizeof (bool) * largeEnough.capacity();
    sum += sizeof (double) * ogMap.outgassingMap.capacity();
    sum += angleMap.GetMemSize();
    return sum;
}

/**
* \brief Initialises geometry properties that haven't been loaded from file
* \return error code: 0=no error, 1=error
*/
int SimulationModel::InitialiseFacets() {
    if (!m.try_lock()) {
        return 1;
    }

    for (const auto& f : facets) {
        auto& facet = *f;
        // Main facet params
        // Current facet
        //SubprocessFacet *f = model->facets[i];
        CalculateFacetParams(&facet);

        // Set some texture parameters
        // bool Facet::SetTexture(double width, double height, bool useMesh)
        if (facet.sh.texWidth_precise * facet.sh.texHeight_precise > 0.0000001) {
            const double ceilCutoff = 0.9999999;
            facet.sh.texWidth = (int) std::ceil(facet.sh.texWidth_precise *
                                                ceilCutoff); //0.9999999: cut the last few digits (convert rounding error 1.00000001 to 1, not 2)
            facet.sh.texHeight = (int) std::ceil(facet.sh.texHeight_precise * ceilCutoff);
        } else {
            facet.sh.texWidth = 0;
            facet.sh.texHeight = 0;
            facet.sh.texWidth_precise = 0.0;
            facet.sh.texHeight_precise = 0.0;
        }
    }

    m.unlock();
    return 0;
}
/*!
 * @brief Calculates various facet parameters without sanity checking @see Geometry::CalculateFacetParams(Facet* f)
 * @param f individual subprocess facet
 */
void SimulationModel::CalculateFacetParams(SubprocessFacet* f) {
    // Calculate facet normal
    Vector3d p0 = vertices3[f->indices[0]];
    Vector3d v1;
    Vector3d v2;
    bool consecutive = true;
    size_t ind = 2;

    // TODO: Handle possible collinear consequtive vectors
    size_t i0 = f->indices[0];
    size_t i1 = f->indices[1];
    while (ind < f->sh.nbIndex && consecutive) {
        size_t i2 = f->indices[ind++];

        v1 = vertices3[i1] - vertices3[i0]; // v1 = P0P1
        v2 = vertices3[i2] - vertices3[i1]; // v2 = P1P2
        f->sh.N = CrossProduct(v1, v2);              // Cross product
        consecutive = (f->sh.N.Norme() < 1e-11);
    }
    f->sh.N = f->sh.N.Normalized();                  // Normalize

    // Calculate Axis Aligned Bounding Box
    f->sh.bb.min = Vector3d(1e100, 1e100, 1e100);
    f->sh.bb.max = Vector3d(-1e100, -1e100, -1e100);

    for (const auto& i : f->indices) {
        const Vector3d& p = vertices3[i];
        f->sh.bb.min.x = std::min(f->sh.bb.min.x,p.x);
        f->sh.bb.min.y = std::min(f->sh.bb.min.y, p.y);
        f->sh.bb.min.z = std::min(f->sh.bb.min.z, p.z);
        f->sh.bb.max.x = std::max(f->sh.bb.max.x, p.x);
        f->sh.bb.max.y = std::max(f->sh.bb.max.y, p.y);
        f->sh.bb.max.z = std::max(f->sh.bb.max.z, p.z);
    }

    // Facet center (AxisAlignedBoundingBox center)
    f->sh.center = 0.5 * (f->sh.bb.max + f->sh.bb.min);

    // Plane equation
    //double A = f->sh.N.x;
    //double B = f->sh.N.y;
    //double C = f->sh.N.z;
    //double D = -Dot(f->sh.N, p0);

    Vector3d p1 = vertices3[f->indices[1]];

    Vector3d U, V;

    U = (p1 - p0).Normalized(); //First side

    // Construct a normal vector V:
    V = CrossProduct(f->sh.N, U); // |U|=1 and |N|=1 => |V|=1

    // u,v vertices (we start with p0 at 0,0)
    f->vertices2[0].u = 0.0;
    f->vertices2[0].v = 0.0;
    Vector2d BBmin; BBmin.u = 0.0; BBmin.v = 0.0;
    Vector2d BBmax; BBmax.u = 0.0; BBmax.v = 0.0;

    for (size_t j = 1; j < f->sh.nbIndex; j++) {
        Vector3d p = vertices3[f->indices[j]];
        Vector3d v = p - p0;
        f->vertices2[j].u = Dot(U, v);  // Project p on U along the V direction
        f->vertices2[j].v = Dot(V, v);  // Project p on V along the U direction

        // Bounds
        BBmax.u  = std::max(BBmax.u , f->vertices2[j].u);
        BBmax.v = std::max(BBmax.v, f->vertices2[j].v);
        BBmin.u = std::min(BBmin.u, f->vertices2[j].u);
        BBmin.v = std::min(BBmin.v, f->vertices2[j].v);
    }

    // Calculate facet area (Meister/Gauss formula)
    double area = 0.0;
    for (size_t j = 0; j < f->sh.nbIndex; j++) {
        size_t j_next = Next(j,f->sh.nbIndex);
        area += f->vertices2[j].u*f->vertices2[j_next].v - f->vertices2[j_next].u*f->vertices2[j].v; //Equal to Z-component of vectorial product
    }
    if (area > 0.0) {

    }
    else if (area < 0.0) {
        //This is a case where a concave facet doesn't obey the right-hand rule:
        //it happens when the first rotation (usually around the second index) is the opposite as the general outline rotation

        //Do a flip
        f->sh.N = -1.0 * f->sh.N;
        V = -1.0 * V;
        BBmin.v = BBmax.v = 0.0;
        for (auto& v : f->vertices2) {
            v.v = -1.0 * v.v;
            BBmax.v = std::max(BBmax.v, v.v);
            BBmin.v = std::min(BBmin.v, v.v);
        }
    }

    f->sh.area = std::abs(0.5 * area);

    // Compute the 2D basis (O,U,V)
    double uD = (BBmax.u - BBmin.u);
    double vD = (BBmax.v - BBmin.v);

    // Origin
    f->sh.O = p0 + BBmin.u * U + BBmin.v * V;

    // Rescale U and V vector
    f->sh.nU = U;
    f->sh.U = U * uD;

    f->sh.nV = V;
    f->sh.V = V * vD;

    f->sh.Nuv = CrossProduct(f->sh.U,f->sh.V); //Not normalized normal vector

    // Rescale u,v coordinates
    for (auto& p : f->vertices2) {
        p.u = (p.u - BBmin.u) / uD;
        p.v = (p.v - BBmin.v) / vD;
    }

#if defined(MOLFLOW)
    f->sh.maxSpeed = 4.0 * std::sqrt(2.0*8.31*f->sh.temperature / 0.001 / wp.gasMass);
#endif
}

int SimulationModel::BuildAccelStructure(GlobalSimuState *globState, int accel_type, BVHAccel::SplitMethod split,
                                         int bvh_width) {
    Chronometer timer;
    timer.Start();

    if (!m.try_lock()) {
        return 1;
    }

#if defined(USE_OLD_BVH)
    std::vector<std::vector<SubprocessFacet*>> facetPointers;
    facetPointers.resize(this->sh.nbSuper);
    for(auto& sFac : this->facets){
        // TODO: Build structures
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : facetPointers) {
                fp_vec.push_back(sFac.get());
            }
        }
        else {
            facetPointers[sFac->sh.superIdx].push_back(sFac.get()); //Assign to structure
        }
    }

    // Build all AABBTrees
    size_t maxDepth=0;
    for (size_t s = 0; s < this->sh.nbSuper; ++s) {
        auto& structure = this->structures[s];
        if(structure.aabbTree)
            structure.aabbTree.reset();
        AABBNODE* tree = BuildAABBTree(facetPointers[s], 0, maxDepth);
        structure.aabbTree = std::make_shared<AABBNODE>(*tree);
        //delete tree; // pointer unnecessary because of make_shared
    }

#else
    std::vector<std::vector<std::shared_ptr<Primitive>>> primPointers;
    primPointers.resize(this->sh.nbSuper);
    for(auto& sFac : this->facets){
        if (sFac->sh.superIdx == -1) { //Facet in all structures
            for (auto& fp_vec : primPointers) {
                fp_vec.push_back(sFac);
            }
        }
        else {
            primPointers[sFac->sh.superIdx].push_back(sFac); //Assign to structure
        }
    }

    for(auto& sFac : this->facets){
        if (sFac->sh.opacity_paramId == -1){ //constant sticking
            sFac->sh.opacity = std::clamp(sFac->sh.opacity, 0.0, 1.0);
            sFac->surf = this->GetSurface(sFac->sh.opacity);
        }
        else {
            auto* par = &this->tdParams.parameters[sFac->sh.opacity_paramId];
            sFac->surf = this->GetParameterSurface(sFac->sh.opacity_paramId, par);
        }
    }

    this->accel.clear();
    if(BVHAccel::SplitMethod::ProbSplit == split && globState && globState->initialized && globState->globalHits.globalHits.nbDesorbed > 0){
        if(globState->facetStates.size() != this->facets.size())
            return 1;
        std::vector<double> probabilities;
        probabilities.reserve(globState->facetStates.size());
        for(auto& state : globState->facetStates) {
            probabilities.emplace_back(state.momentResults[0].hits.nbHitEquiv / globState->globalHits.globalHits.nbHitEquiv);
        }
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == 1)
                this->accel.emplace_back(std::make_shared<KdTreeAccel>(primPointers[s], probabilities));
            else
                this->accel.emplace_back(std::make_shared<BVHAccel>(primPointers[s], bvh_width, BVHAccel::SplitMethod::ProbSplit, probabilities));
        }
    }
    else {
        for (size_t s = 0; s < this->sh.nbSuper; ++s) {
            if(accel_type == 1)
                this->accel.emplace_back(std::make_shared<KdTreeAccel>(primPointers[s]));
            else
                this->accel.emplace_back(std::make_shared<BVHAccel>(primPointers[s], bvh_width, split));
        }
    }
#endif // old_bvb

    timer.Stop();
    m.unlock();

    return 0;
}

/**
* \brief Do calculations necessary before launching simulation
* determine latest moment
* Generate integrated desorption functions
* match parameters
* Generate speed distribution functions
* Angle map
*/
int SimulationModel::PrepareToRun() {
    if (!m.try_lock()) {
        return 1;
    }
    initialized = false;

    std::string errLog;

    //determine latest moment
    wp.latestMoment = 1E-10;
    if(!tdParams.moments.empty())
        wp.latestMoment = (tdParams.moments.end()-1)->second;
        //wp.latestMoment = (tdParams.moments.end()-1)->first + (tdParams.moments.end()-1)->second / 2.0;

    std::set<size_t> desorptionParameterIDs;

    //Check and calculate various facet properties for time dependent simulations (CDF, ID )
    for (size_t i = 0; i < sh.nbFacet; i++) {
        SubprocessFacet& facet = *facets[i];
        // TODO: Find a solution to integrate catalog parameters
        if(facet.sh.outgassing_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Outgassing parameter \"%d\" isn't defined.", i + 1, facet.sh.outgassing_paramId);
            errLog.append(tmp);
        }
        if(facet.sh.opacity_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Opacity parameter \"%d\" isn't defined.", i + 1, facet.sh.opacity_paramId);
            errLog.append(tmp);
        }
        if(facet.sh.sticking_paramId >= (int) tdParams.parameters.size()){
            char tmp[256];
            sprintf(tmp, "Facet #%zd: Sticking parameter \"%d\" isn't defined.", i + 1, facet.sh.sticking_paramId);
            errLog.append(tmp);
        }

        if (facet.sh.outgassing_paramId >= 0) { //if time-dependent desorption
            int id = IDGeneration::GetIDId(desorptionParameterIDs, facet.sh.outgassing_paramId);
            if (id >= 0)
                facet.sh.IDid = id; //we've already generated an ID for this temperature
            else {
                auto[id_new, id_vec] = IDGeneration::GenerateNewID(desorptionParameterIDs, facet.sh.outgassing_paramId, this);
                facet.sh.IDid = id_new;
                tdParams.IDs.emplace_back(std::move(id_vec));
            }
        }

        // Generate speed distribution functions
        std::list<double> temperatureList;
        int id = CDFGeneration::GetCDFId(temperatureList, facet.sh.temperature);
        if (id >= 0)
            facet.sh.CDFid = id; //we've already generated a CDF for this temperature
        else {
            auto[cdf_id, cdf_vec] = CDFGeneration::GenerateNewCDF(temperatureList, facet.sh.temperature, wp.gasMass);
            facet.sh.CDFid = cdf_id;
            tdParams.CDFs.emplace_back(cdf_vec);
        }
        //Angle map
        if (facet.sh.desorbType == DES_ANGLEMAP) {
            if (!facet.sh.anglemapParams.hasRecorded) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Uses angle map desorption but doesn't have a recorded angle map.", i + 1);
                errLog.append(tmp);
            }
            if (facet.sh.anglemapParams.record) {
                char tmp[256];
                sprintf(tmp, "Facet #%zd: Can't RECORD and USE angle map desorption at the same time.", i + 1);
                errLog.append(tmp);
            }
        }
    }

    if(!errLog.empty()){
        m.unlock();
        return 1;
    }

    CalcTotalOutgassing();

    initialized = true;
    m.unlock();

    return 0;
}

/**
* \brief Compute the outgassing of all source facet depending on the mode (file, regular, time-dependent) and set it to the global settings
*/
void SimulationModel::CalcTotalOutgassing() {
    // Compute the outgassing of all source facet
    double totalDesorbedMolecules = 0.0;
    double finalOutgassingRate_Pa_m3_sec = 0.0;
    double finalOutgassingRate = 0.0;

    const double latestMoment = wp.latestMoment;

    for (size_t i = 0; i < sh.nbFacet; i++) {
        SubprocessFacet& facet = *facets[i];
        if (facet.sh.desorbType != DES_NONE) { //there is a kind of desorption
            if (facet.sh.useOutgassingFile) { //outgassing file
                auto& ogMap = facet.ogMap;
                for (size_t l = 0; l < (ogMap.outgassingMapWidth * ogMap.outgassingMapHeight); l++) {
                    totalDesorbedMolecules += latestMoment * ogMap.outgassingMap[l] / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate += ogMap.outgassingMap[l] / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate_Pa_m3_sec += ogMap.outgassingMap[l];
                }
            } else { //regular outgassing
                if (facet.sh.outgassing_paramId == -1) { //constant outgassing
                    totalDesorbedMolecules += latestMoment * facet.sh.outgassing / (1.38E-23 * facet.sh.temperature);
                    finalOutgassingRate +=
                            facet.sh.outgassing / (1.38E-23 * facet.sh.temperature);  //Outgassing molecules/sec
                    finalOutgassingRate_Pa_m3_sec += facet.sh.outgassing;
                } else { //time-dependent outgassing
                    totalDesorbedMolecules += tdParams.IDs[facet.sh.IDid].back().second / (1.38E-23 * facet.sh.temperature);
                    size_t lastIndex = tdParams.parameters[facet.sh.outgassing_paramId].GetSize() - 1;
                    double finalRate_mbar_l_s = tdParams.parameters[facet.sh.outgassing_paramId].GetY(lastIndex);
                    finalOutgassingRate +=
                            finalRate_mbar_l_s * 0.100 / (1.38E-23 * facet.sh.temperature); //0.1: mbar*l/s->Pa*m3/s
                    finalOutgassingRate_Pa_m3_sec += finalRate_mbar_l_s * 0.100;
                }
            }
        }
    }

    wp.totalDesorbedMolecules = totalDesorbedMolecules;
    wp.finalOutgassingRate_Pa_m3_sec = finalOutgassingRate_Pa_m3_sec;
    wp.finalOutgassingRate = finalOutgassingRate;
}

SimulationModel::~SimulationModel() = default;

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator=(const GlobalSimuState & src) {
    //Copy all but mutex
    facetStates = src.facetStates;
    globalHistograms = src.globalHistograms;
    globalHits = src.globalHits;
    initialized = src.initialized;
    return *this;
}

/**
* \brief Assign operator
* \param src reference to source object
* \return address of this
*/
GlobalSimuState& GlobalSimuState::operator+=(const GlobalSimuState & src) {
    //Copy all but mutex
    facetStates += src.facetStates;
    globalHistograms += src.globalHistograms;
    globalHits += src.globalHits;
    return *this;
}

/**
* \brief Clears simulation state
*/
void GlobalSimuState::clear() {
    tMutex.lock();
    globalHits = GlobalHitBuffer();
    globalHistograms.clear();
    facetStates.clear();
    initialized = false;
    tMutex.unlock();
}

/**
* \brief Constructs the 'dpHit' structure to hold all results, zero-init
* \param w Worker handle
*/
void GlobalSimuState::Resize(const SimulationModel &model) { //Constructs the 'dpHit' structure to hold all results, zero-init
    //LockMutex(mutex);
    tMutex.lock();
    size_t nbF = model.sh.nbFacet;
    size_t nbMoments = model.tdParams.moments.size();
    facetStates.assign(nbF, FacetState());

    if(!model.facets.empty()) {
        for (size_t i = 0; i < nbF; i++) {
            auto sFac = model.facets[i];
            if (sFac->globalId != i) {
                std::cerr << fmt::format("Facet ID mismatch! : {} / {}\n", sFac->globalId , i);
                tMutex.unlock();
                exit(0);
            }

            FacetMomentSnapshot facetMomentTemplate{};
            facetMomentTemplate.histogram.Resize(sFac->sh.facetHistogramParams);
            facetMomentTemplate.direction.assign((sFac->sh.countDirection ? sFac->sh.texWidth*sFac->sh.texHeight : 0), DirectionCell());
            facetMomentTemplate.profile.assign((sFac->sh.isProfile ? PROFILE_SIZE : 0), ProfileSlice());
            facetMomentTemplate.texture.assign((sFac->sh.isTextured ? sFac->sh.texWidth*sFac->sh.texHeight : 0), TextureCell());

            //No init for hits
            facetStates[i].momentResults.assign(1 + nbMoments, facetMomentTemplate);
            if (sFac->sh.anglemapParams.record)
                facetStates[i].recordedAngleMapPdf.assign(sFac->sh.anglemapParams.GetMapSize(), 0);
        }
    }
    /*for (size_t i = 0; i < nbF; i++) {
        const SubprocessFacet& sFac = facets[i];
        FacetMomentSnapshot facetMomentTemplate;
        facetMomentTemplate.histogram.Resize(sFac.sh.facetHistogramParams);
        facetMomentTemplate.direction = std::vector<DirectionCell>(sFac.sh.countDirection ? sFac.sh.texWidth*sFac.sh.texHeight : 0);
        facetMomentTemplate.profile = std::vector<ProfileSlice>(sFac.sh.isProfile ? PROFILE_SIZE : 0);
        facetMomentTemplate.texture = std::vector<TextureCell>(sFac.sh.isTextured ? sFac.sh.texWidth*sFac.sh.texHeight : 0);
        //No init for hits
        facetStates[i].momentResults = std::vector<FacetMomentSnapshot>(1 + nbMoments, facetMomentTemplate);
        if (sFac.sh.anglemapParams.record) facetStates[i].recordedAngleMapPdf = std::vector<size_t>(sFac.sh.anglemapParams.GetMapSize());
    }*/
    //Global histogram

    FacetHistogramBuffer globalHistTemplate{};
    globalHistTemplate.Resize(model.wp.globalHistogramParams);
    globalHistograms.assign(1 + nbMoments, globalHistTemplate);
    initialized = true;
    //ReleaseMutex(mutex);
    tMutex.unlock();
}

/**
* \brief zero-init for all structures
*/
void GlobalSimuState::Reset() {
    //LockMutex(mutex);
    tMutex.lock();
    for (auto& h : globalHistograms) {
        ZEROVECTOR(h.distanceHistogram);
        ZEROVECTOR(h.nbHitsHistogram);
        ZEROVECTOR(h.timeHistogram);
    }
    memset(&globalHits, 0, sizeof(globalHits)); //Plain old data
    for (auto& state : facetStates) {
        ZEROVECTOR(state.recordedAngleMapPdf);
        for (auto& m : state.momentResults) {
            ZEROVECTOR(m.histogram.distanceHistogram);
            ZEROVECTOR(m.histogram.nbHitsHistogram);
            ZEROVECTOR(m.histogram.timeHistogram);
            std::fill(m.direction.begin(), m.direction.end(), DirectionCell());
            std::fill(m.texture.begin(), m.texture.end(), TextureCell());
            std::fill(m.profile.begin(), m.profile.end(), ProfileSlice());
            memset(&(m.hits), 0, sizeof(m.hits));
        }
    }
    tMutex.unlock();
    //ReleaseMutex(mutex);
}

/**
 * @brief Compare function for two simulation states
 * @param lhsGlobHit first simulation state
 * @param rhsGlobHit second simulation state
 * @param globThreshold threshold for relative difference on global counters
 * @param locThreshold threshold for relative difference on facet local counters
 * @return a tuple containing the number of global, local (facet), and fine (facet profile/texture) errors
 */
std::tuple<int, int, int>
GlobalSimuState::Compare(const GlobalSimuState &lhsGlobHit, const GlobalSimuState &rhsGlobHit, double globThreshold,
                         double locThreshold) {

    const double velocityThresholdFactor = 40.0;
    //std::ofstream cmpFile("cmpFile.txt");
    size_t globalErrNb = 0;
    size_t facetErrNb = 0;
    size_t fineErrNb = 0;

    std::string cmpFile;
    std::string cmpFileFine; // extra stream to silence after important outputs

    // Sanity check
    {
        if(lhsGlobHit.globalHits.globalHits.nbDesorbed == 0 && rhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Neither state has recorded desorptions\n");
            ++globalErrNb;
        }
        else if (lhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] First state has no recorded desorptions\n");
            ++globalErrNb;
        }
        else if (rhsGlobHit.globalHits.globalHits.nbDesorbed == 0){
            cmpFile += fmt::format("[Global][desorp] Second state has no recorded desorptions\n");
            ++globalErrNb;
        }

        if(globalErrNb){
            fmt::print("{}\n", cmpFile);
            return std::make_tuple(static_cast<int>(globalErrNb), -1, -1);
        }
    }

    {
        double absRatio = lhsGlobHit.globalHits.globalHits.nbAbsEquiv / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed);
        double absRatio_rhs = rhsGlobHit.globalHits.globalHits.nbAbsEquiv / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
        if (!IsEqual(absRatio, absRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][absRatio] has large difference: {} ({})\n",
                                   std::abs(absRatio - absRatio_rhs), std::abs(absRatio - absRatio_rhs) / std::max(absRatio, absRatio_rhs));
            ++globalErrNb;
        }
    }

    {
        double hitRatio = static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed);
        double hitRatio_rhs = static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit) / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
        if (!IsEqual(hitRatio, hitRatio_rhs, globThreshold)) {
            cmpFile += fmt::format("[Global][hitRatio] has large difference: {} ({}) --> "
                                   "{} / {} vs {} / {}\n",
                                   std::abs(hitRatio - hitRatio_rhs), std::abs(hitRatio - hitRatio_rhs) / std::max(hitRatio, hitRatio_rhs),
                                   lhsGlobHit.globalHits.globalHits.nbMCHit, lhsGlobHit.globalHits.globalHits.nbDesorbed,
                                   rhsGlobHit.globalHits.globalHits.nbMCHit, rhsGlobHit.globalHits.globalHits.nbDesorbed);

            ++globalErrNb;
        }
    }

    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_v_ort, rhsGlobHit.globalHits.globalHits.sum_v_ort, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_v_ort] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_v_ort - rhsGlobHit.globalHits.globalHits.sum_v_ort));
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_1_per_velocity, rhsGlobHit.globalHits.globalHits.sum_1_per_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_velocity] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_1_per_velocity - rhsGlobHit.globalHits.globalHits.sum_1_per_velocity));
        ++globalErrNb;
    }
    if (!IsEqual(lhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity, rhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity, globThreshold)) {
        cmpFile += fmt::format("[Global][sum_1_per_ort_velocity] has large difference: {}\n",
                               std::abs(lhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity - rhsGlobHit.globalHits.globalHits.sum_1_per_ort_velocity));
        ++globalErrNb;
    }


    // Histogram
    {
        auto& hist_lhs = lhsGlobHit.globalHistograms;
        auto& hist_rhs = rhsGlobHit.globalHistograms;
        for(size_t tHist = 0; tHist < hist_lhs.size(); tHist++) {
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].nbHitsHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].nbHitsHistogram[hIndex], hist_rhs[tHist].nbHitsHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Bounces][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].nbHitsHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].distanceHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].distanceHistogram[hIndex], hist_rhs[tHist].distanceHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Dist][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].distanceHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].distanceHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));
                    ++globalErrNb;
                }
            }
            for (size_t hIndex = 0; hIndex < hist_lhs[tHist].timeHistogram.size(); ++hIndex) {
                if(std::sqrt(std::max(1.0,std::min(hist_lhs[tHist].timeHistogram[hIndex], hist_rhs[tHist].timeHistogram[hIndex]))) < 80) {
                    // Sample size not large enough
                    continue;
                }
                if (!IsEqual(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit), hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit), locThreshold)) {
                    cmpFile += fmt::format("[Global][Hist][Time][Ind={}] has large difference: {}\n",
                                           hIndex, std::abs(hist_lhs[tHist].timeHistogram[hIndex] / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbMCHit) - hist_rhs[tHist].timeHistogram[hIndex] / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbMCHit)));

                    ++globalErrNb;
                }
            }
        }
    }
    // facets

    auto locThreshold_bak = locThreshold;
    for(int facetId = 0; facetId < lhsGlobHit.facetStates.size(); ++facetId)
    {//cmp
        if(lhsGlobHit.facetStates[facetId].momentResults.size() != rhsGlobHit.facetStates[facetId].momentResults.size()){
            cmpFile += fmt::format("[Facet][{}] Different amount of moments for each state: {} vs {}\n", facetId, lhsGlobHit.facetStates[facetId].momentResults.size(), rhsGlobHit.facetStates[facetId].momentResults.size());
            ++facetErrNb;
            continue;
        }

        for(int m = 0; m < lhsGlobHit.facetStates[facetId].momentResults.size(); ++m) {
            auto &facetCounter_lhs = lhsGlobHit.facetStates[facetId].momentResults[m];
            auto &facetCounter_rhs = rhsGlobHit.facetStates[facetId].momentResults[m];

            // If one facet doesn't have any hits recorded, comparison is pointless, so just skip to next facet
            if ((facetCounter_lhs.hits.nbMCHit == 0 && facetCounter_rhs.hits.nbMCHit == 0)
            || (facetCounter_lhs.hits.nbMCHit < 40000 || facetCounter_rhs.hits.nbMCHit < 40000)) {
                // Skip facet comparison if not enough hits have been recorded for both states
                continue;
            } else if (facetCounter_lhs.hits.nbMCHit == 0 && facetCounter_rhs.hits.nbMCHit > 0) {
                cmpFile += fmt::format("[Facet][{}][hits][{}] First state has no recorded hits for this facet\n", facetId, m);
                ++facetErrNb;
                continue;
            } else if (facetCounter_lhs.hits.nbMCHit > 0 && facetCounter_rhs.hits.nbMCHit == 0) {
                cmpFile += fmt::format("[Facet][{}][hits][{}] Second state has no recorded hits for this facet\n", facetId, m);
                ++facetErrNb;
                continue;
            }

            // Adjust threshold to reasonable limit for moments due to lower amount of hits
            if(m >= 0){
                auto nMC = std::min(lhsGlobHit.facetStates[facetId].momentResults[m].hits.nbMCHit, rhsGlobHit.facetStates[facetId].momentResults[m].hits.nbMCHit);
                locThreshold = std::max(locThreshold_bak, 1.0 / std::sqrt(nMC));
            }
            else{
                locThreshold = locThreshold_bak;
            }

            double scale = 1.0 / static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed); // getmolpertp
            double scale_rhs = 1.0 / static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed);
            double fullScale = 1.0;

            // density correction value for scale
            if (facetCounter_lhs.hits.nbMCHit > 0 || facetCounter_lhs.hits.nbDesorbed > 0) {
                if (facetCounter_lhs.hits.nbAbsEquiv > 0.0 ||
                    facetCounter_lhs.hits.nbDesorbed > 0) { //otherwise, save calculation time
                    fullScale = 1.0 - (facetCounter_lhs.hits.nbAbsEquiv + (double) facetCounter_lhs.hits.nbDesorbed) /
                                      (facetCounter_lhs.hits.nbHitEquiv + (double) facetCounter_lhs.hits.nbDesorbed) /
                                      2.0;
                }
            }

            double fullScale_rhs = 1.0;
            if (facetCounter_rhs.hits.nbMCHit > 0 || facetCounter_rhs.hits.nbDesorbed > 0) {
                if (facetCounter_rhs.hits.nbAbsEquiv > 0.0 ||
                    facetCounter_rhs.hits.nbDesorbed > 0) { //otherwise, save calculation time
                    fullScale_rhs = 1.0 -
                                    (facetCounter_rhs.hits.nbAbsEquiv + (double) facetCounter_rhs.hits.nbDesorbed) /
                                    (facetCounter_rhs.hits.nbHitEquiv + (double) facetCounter_rhs.hits.nbDesorbed) /
                                    2.0;
                }
            }

            fullScale *= scale;
            fullScale_rhs *= scale_rhs;

            scale = 1.0 / lhsGlobHit.globalHits.globalHits.nbHitEquiv;
            scale_rhs = 1.0 / rhsGlobHit.globalHits.globalHits.nbHitEquiv;
            fullScale = 1.0 /
                        (lhsGlobHit.globalHits.globalHits.nbHitEquiv + lhsGlobHit.globalHits.globalHits.nbAbsEquiv +
                         static_cast<double>(lhsGlobHit.globalHits.globalHits.nbDesorbed));
            fullScale_rhs = 1.0 /
                            (rhsGlobHit.globalHits.globalHits.nbHitEquiv + rhsGlobHit.globalHits.globalHits.nbAbsEquiv +
                             static_cast<double>(rhsGlobHit.globalHits.globalHits.nbDesorbed));
            double sumHitDes = facetCounter_lhs.hits.nbHitEquiv + static_cast<double>(facetCounter_lhs.hits.nbDesorbed);
            double sumHitDes_rhs =
                    facetCounter_rhs.hits.nbHitEquiv + static_cast<double>(facetCounter_rhs.hits.nbDesorbed);

            if (!(std::sqrt(
                    std::max(1.0, std::min(facetCounter_lhs.hits.nbHitEquiv, facetCounter_rhs.hits.nbHitEquiv))) <
                  80)) {
                double hitRatio = facetCounter_lhs.hits.nbHitEquiv * scale;
                double hitRatio_rhs = facetCounter_rhs.hits.nbHitEquiv * scale_rhs;
                if (!IsEqual(hitRatio, hitRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][hitRatio][{}] has large difference: "
                                           "{} (normalized: {}) --> "
                                           "{} / {} vs {} / {}\n",
                                           facetId, m,
                                           std::abs(hitRatio - hitRatio_rhs), std::abs(hitRatio - hitRatio_rhs) / std::max(hitRatio, hitRatio_rhs),
                                           facetCounter_lhs.hits.nbHitEquiv, lhsGlobHit.globalHits.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.nbHitEquiv, lhsGlobHit.globalHits.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs,
                             locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][sum_v_ort][{}] has large difference: "
                                           "{} (normalized: {}) --> "
                                           "{} / {} vs {} / {}\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_v_ort * scale - facetCounter_rhs.hits.sum_v_ort * scale_rhs), std::abs(facetCounter_lhs.hits.sum_v_ort * scale - facetCounter_rhs.hits.sum_v_ort * scale_rhs) / std::max(facetCounter_lhs.hits.sum_v_ort * scale, facetCounter_rhs.hits.sum_v_ort * scale_rhs),
                                           facetCounter_lhs.hits.sum_v_ort, lhsGlobHit.globalHits.globalHits.nbHitEquiv,
                                           facetCounter_rhs.hits.sum_v_ort, lhsGlobHit.globalHits.globalHits.nbHitEquiv);
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_velocity][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs),
                                           std::abs(facetCounter_lhs.hits.sum_1_per_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs) / std::max(facetCounter_lhs.hits.sum_1_per_velocity * fullScale, facetCounter_rhs.hits.sum_1_per_velocity * fullScale_rhs));
                    ++facetErrNb;
                }
                if (!IsEqual(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale,
                             facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs,
                             locThreshold * velocityThresholdFactor)) {
                    cmpFile += fmt::format("[Facet][{}][sum_1_per_ort_velocity][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs),
                                           std::abs(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale - facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs) / std::max(facetCounter_lhs.hits.sum_1_per_ort_velocity * fullScale, facetCounter_rhs.hits.sum_1_per_ort_velocity * fullScale_rhs));
                    ++facetErrNb;
                }
            }

            if (!(std::sqrt(
                    std::max(1.0, std::min(facetCounter_lhs.hits.nbAbsEquiv, facetCounter_rhs.hits.nbAbsEquiv))) <
                  80)) {
                double absRatio = facetCounter_lhs.hits.nbAbsEquiv / static_cast<double>(facetCounter_lhs.hits.nbMCHit);
                double absRatio_rhs = facetCounter_rhs.hits.nbAbsEquiv / static_cast<double>(facetCounter_rhs.hits.nbMCHit);
                if (!IsEqual(absRatio, absRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][absRatio][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(absRatio - absRatio_rhs), std::abs(absRatio - absRatio_rhs) / std::max(absRatio, absRatio_rhs));
                    ++facetErrNb;
                }
            }

            if (!(std::sqrt(std::max((size_t) 1,
                                     std::min(facetCounter_lhs.hits.nbDesorbed, facetCounter_rhs.hits.nbDesorbed))) <
                  80)) {
                double desRatio = (double) facetCounter_lhs.hits.nbDesorbed / static_cast<double>(facetCounter_lhs.hits.nbMCHit);
                double desRatio_rhs = (double) facetCounter_rhs.hits.nbDesorbed / static_cast<double>(facetCounter_rhs.hits.nbMCHit);
                if (!IsEqual(desRatio, desRatio_rhs, locThreshold)) {
                    cmpFile += fmt::format("[Facet][{}][desRatio][{}] has large difference: "
                                           "{} (normalized: {})\n",
                                           facetId, m,
                                           std::abs(desRatio - desRatio_rhs), std::abs(desRatio - desRatio_rhs) / std::max(desRatio, desRatio_rhs));
                    ++facetErrNb;
                }
            }

            //profile
            {
                auto &prof_lhs = facetCounter_lhs.profile;
                auto &prof_rhs = facetCounter_rhs.profile;

                for (int id = 0; id < prof_lhs.size(); ++id) {
                    if (std::sqrt(std::max(1.0, std::min(prof_lhs[id].countEquiv, prof_rhs[id].countEquiv))) < 10) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual(prof_lhs[id].countEquiv / sumHitDes, prof_rhs[id].countEquiv / sumHitDes_rhs,
                                 locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][countEquiv][{}] has large difference: "
                                               "{} : {} - {}\n",
                                               facetId, id, m,
                                               std::abs(prof_lhs[id].countEquiv / sumHitDes - prof_rhs[id].countEquiv / sumHitDes_rhs) / (prof_lhs[id].countEquiv / sumHitDes), 
                                               std::abs(prof_lhs[id].countEquiv / sumHitDes), (prof_rhs[id].countEquiv / sumHitDes_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(prof_lhs[id].sum_1_per_ort_velocity * scale,
                                 prof_rhs[id].sum_1_per_ort_velocity * scale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][sum_1_per_ort_velocity][{}] has large rel difference: "
                                               "{} : {} - {}\n",
                                               facetId, id, m,
                                               std::abs(prof_lhs[id].sum_1_per_ort_velocity * scale - prof_rhs[id].sum_1_per_ort_velocity * scale_rhs) / (prof_lhs[id].sum_1_per_ort_velocity * scale),
                                               std::abs(prof_lhs[id].sum_1_per_ort_velocity * scale), (prof_rhs[id].sum_1_per_ort_velocity * scale_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(prof_lhs[id].sum_v_ort * scale, prof_rhs[id].sum_v_ort * scale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Profile][Ind={}][sum_v_ort][{}] has large rel difference: "
                                               "{} : {} - {}\n",
                                               facetId, id, m,
                                               std::abs(prof_lhs[id].sum_v_ort * scale - prof_rhs[id].sum_v_ort * scale_rhs) / (prof_lhs[id].sum_v_ort * scale),
                                               std::abs(prof_lhs[id].sum_v_ort * scale), (prof_rhs[id].sum_v_ort * scale_rhs));

                        ++fineErrNb;
                    }
                }
            }

            //texture
            {
                auto &tex_lhs = facetCounter_lhs.texture;
                auto &tex_rhs = facetCounter_rhs.texture;
                int ix = 0;
                for (int iy = 0; iy < tex_lhs.size(); iy++) {
                    //for (int ix = 0; ix < texWidth_file; ix++) {
                    if (std::sqrt(std::max(1.0, std::min(tex_lhs[iy].countEquiv, tex_rhs[iy].countEquiv))) < 80) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual(tex_lhs[iy].countEquiv / sumHitDes, tex_rhs[iy].countEquiv / sumHitDes_rhs,
                                 locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][countEquiv][{}] has large rel difference: "
                                               "{} : {} - {}\n",
                                               facetId, ix, iy, m,
                                               std::abs( tex_lhs[iy].countEquiv / sumHitDes - tex_rhs[iy].countEquiv / sumHitDes_rhs) / (tex_lhs[iy].countEquiv / sumHitDes),
                                               std::abs(tex_lhs[iy].countEquiv / sumHitDes), (tex_rhs[iy].countEquiv / sumHitDes_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(tex_lhs[iy].sum_1_per_ort_velocity * fullScale,
                                 tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][sum_1_per_ort_velocity][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs(tex_lhs[iy].sum_1_per_ort_velocity * fullScale - tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs) / (tex_lhs[iy].sum_1_per_ort_velocity * fullScale),
                                                   std::abs(tex_lhs[iy].sum_1_per_ort_velocity * fullScale), (tex_rhs[iy].sum_1_per_ort_velocity * fullScale_rhs));

                        ++fineErrNb;
                    }
                    if (!IsEqual(tex_lhs[iy].sum_v_ort_per_area * scale, tex_rhs[iy].sum_v_ort_per_area * scale_rhs,
                                 locThreshold * velocityThresholdFactor)) {
                        cmpFileFine += fmt::format("[Facet][{}][Texture][{},{}][sum_v_ort_per_area][{}] has large rel difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs(tex_lhs[iy].sum_v_ort_per_area * scale - tex_rhs[iy].sum_v_ort_per_area * scale_rhs) / (tex_lhs[iy].sum_v_ort_per_area * scale),
                                                   std::abs(tex_lhs[iy].sum_v_ort_per_area * scale), (tex_rhs[iy].sum_v_ort_per_area * scale_rhs));

                        ++fineErrNb;
                    }
                    //}
                } // end for comp texture
            }

            //Directions
            {
                auto &dir_lhs = facetCounter_lhs.direction;
                auto &dir_rhs = facetCounter_rhs.direction;
                int ix = 0;
                for (int iy = 0; iy < dir_lhs.size(); iy++) {
                    //for (int ix = 0; ix < dirWidth_file; ix++) {
                    if (std::sqrt(std::max(1.0, (double) std::min(dir_lhs[iy].count, dir_rhs[iy].count))) < 80) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual((double)dir_lhs[iy].count, (double)dir_rhs[iy].count, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][count][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].count - (int) dir_rhs[iy].count),
                                                   (int) dir_lhs[iy].count, (int) dir_rhs[iy].count);

                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.x, dir_rhs[iy].dir.x, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.x][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.x - (int) dir_rhs[iy].dir.x),
                                                   (int) dir_lhs[iy].dir.x, (int) dir_rhs[iy].dir.x);
                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.y, dir_rhs[iy].dir.y, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.y][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.y - (int) dir_rhs[iy].dir.y),
                                                   (int) dir_lhs[iy].dir.y, (int) dir_rhs[iy].dir.y);
                        ++fineErrNb;
                    }
                    if (!IsEqual(dir_lhs[iy].dir.z, dir_rhs[iy].dir.z, locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][dirs][{},{}][dir.z][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, ix, iy, m,
                                                   std::abs((int) dir_lhs[iy].dir.z - (int) dir_rhs[iy].dir.z),
                                                   (int) dir_lhs[iy].dir.z, (int) dir_rhs[iy].dir.z);
                        ++fineErrNb;
                    }
                    //}
                } // end for comp dir
            }

            //facet hist
            {
                auto &hist_lhs = facetCounter_lhs.histogram;
                auto &hist_rhs = facetCounter_rhs.histogram;
                for (size_t hIndex = 0; hIndex < hist_lhs.nbHitsHistogram.size(); ++hIndex) {
                    if (std::sqrt(std::max(1.0, std::min(hist_lhs.nbHitsHistogram[hIndex],
                                                         hist_rhs.nbHitsHistogram[hIndex]))) < 80) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual(hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Bounces][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.nbHitsHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));
                        ++fineErrNb;
                    }
                }

                for (size_t hIndex = 0; hIndex < hist_lhs.distanceHistogram.size(); ++hIndex) {
                    if (std::sqrt(std::max(1.0, std::min(hist_lhs.distanceHistogram[hIndex],
                                                         hist_rhs.distanceHistogram[hIndex]))) < 80) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual(hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Dist][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.distanceHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));

                        ++fineErrNb;
                    }
                }

                for (size_t hIndex = 0; hIndex < hist_lhs.timeHistogram.size(); ++hIndex) {
                    if (std::sqrt(
                            std::max(1.0, std::min(hist_lhs.timeHistogram[hIndex], hist_rhs.timeHistogram[hIndex]))) <
                        80) {
                        // Sample size not large enough
                        continue;
                    }
                    if (!IsEqual(hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit),
                                 hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit), locThreshold)) {
                        cmpFileFine += fmt::format("[Facet][{}][Hist][Time][Ind={}][{}] has large difference: "
                                                   "{} : {} - {}\n",
                                                   facetId, hIndex, m,
                                                   std::abs(hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit) -
                                                            hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit)),
                                                   hist_lhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_lhs.hits.nbMCHit), hist_rhs.timeHistogram[hIndex] / static_cast<double>(facetCounter_rhs.hits.nbMCHit));

                        ++fineErrNb;
                    }
                }
            }
        }
    }

    std::string cmp_string;
    int i = 0;
    {
        std::istringstream compStream(cmpFile);
        for (; i < 32 && std::getline(compStream, cmp_string, '\n'); ++i) {
            Log::console_error("%s\n", cmp_string.c_str());
        }
    }
    {
        std::istringstream compStreamFine(cmpFileFine);
        for (; i < 32 && std::getline(compStreamFine, cmp_string, '\n'); ++i) {
            Log::console_msg_master(4, "%s\n", cmp_string.c_str());
        }
    }

    if(i >= 32) {
        Log::console_error("[Warning] List of differences too long: Total = %lu\n", globalErrNb + facetErrNb + fineErrNb);
    }

    return std::make_tuple(globalErrNb, facetErrNb, fineErrNb);
}

/**
* \brief Resize histograms according to sizes in params
* \param params contains data about sizes
*/
/*void FacetHistogramBuffer::Resize(const HistogramParams& params) {
    nbHitsHistogram = std::vector<double>(params.recordBounce ? params.GetBounceHistogramSize() : 0);
    distanceHistogram = std::vector<double>(params.recordDistance ? params.GetDistanceHistogramSize() : 0);
    timeHistogram = std::vector<double>(params.recordTime ? params.GetTimeHistogramSize() : 0);
}*/

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetHistogramBuffer& FacetHistogramBuffer::operator+=(const FacetHistogramBuffer & rhs) {
    // if (model.wp.globalHistogramParams.recordBounce)
    this->nbHitsHistogram += rhs.nbHitsHistogram;
    this->distanceHistogram += rhs.distanceHistogram;
    this->timeHistogram += rhs.timeHistogram;
    return *this;
}

FacetMomentSnapshot::FacetMomentSnapshot() : hits(), histogram(){

}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetMomentSnapshot& FacetMomentSnapshot::operator+=(const FacetMomentSnapshot & rhs) {
    this->hits += rhs.hits;
    this->profile += rhs.profile;
    this->texture += rhs.texture;
    this->direction += rhs.direction;
    this->histogram += rhs.histogram;
    return *this;
}

/**
* \brief + operator, simply calls implemented +=
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetMomentSnapshot& FacetMomentSnapshot::operator+(const FacetMomentSnapshot & rhs) {
    *this += rhs;
    return *this;
}

/**
* \brief += operator, with simple += of underlying structures
* \param rhs reference object on the right hand
* \return address of this (lhs)
*/
FacetState& FacetState::operator+=(const FacetState & rhs) {
    // Check in case simulation pdf is empty (record==false) but global pdf is not (hasRecorded==true)
    if(this->recordedAngleMapPdf.size() == rhs.recordedAngleMapPdf.size())
        this->recordedAngleMapPdf += rhs.recordedAngleMapPdf;
    this->momentResults += rhs.momentResults;
    return *this;
}