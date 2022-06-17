/*
Program:     MolFlow+ / Molflow+
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

#include "MolflowSimFacet.h"
#include <fmt/core.h>

MolflowSimFacet::MolflowSimFacet(const MolflowSimFacet& cpy)  : SimulationFacet(cpy) {
    *this = cpy;
}

MolflowSimFacet::MolflowSimFacet(MolflowSimFacet&& cpy) noexcept : SimulationFacet(cpy){
    *this = std::move(cpy);
}

MolflowSimFacet& MolflowSimFacet::operator=(const MolflowSimFacet& cpy){
    SimulationFacet::operator=(cpy);
    this->angleMap = cpy.angleMap;
    this->ogMap = cpy.ogMap;

    return *this;
}

MolflowSimFacet& MolflowSimFacet::operator=(MolflowSimFacet&& cpy) noexcept {
    SimulationFacet::operator=(cpy);
    this->angleMap = cpy.angleMap;
    this->ogMap = cpy.ogMap;

    return *this;
}

/**
* \brief Initialise various properties for a single facet
* \param id global facet id for the application
* \param regions vector containing the regions and their parameters used to fetch the min/max energy for a spectrum
* \return boolean that is true on successful init
*/
bool MolflowSimFacet::InitializeOnLoad(const size_t &id, const size_t &nbMoments) {
    globalId = id;
    if (!InitializeLinkAndVolatile(id)) return false;
    InitializeOutgassingMap();

    if(InitializeAngleMap() < 0)
        return false;

    InitializeTexture(nbMoments);
    InitializeProfile(nbMoments);
    InitializeDirectionTexture(nbMoments);
    InitializeHistogram(nbMoments);

    return true;
}

int MolflowSimFacet::InitializeAngleMap()
{
    //Incident angle map
    int angleMapSize = 0;
    if (sh.desorbType == DES_ANGLEMAP) { //Use mode
        if (angleMap.pdf.empty()) // check if a valid recorded angle map exists, otherwise we can't desorb
            throw Error(fmt::format("Facet {}: should generate by angle map but has none recorded.", globalId + 1).c_str());

        //Construct CDFs
        try {
            angleMap.phi_pdfsums_lowerTheta.resize(sh.anglemapParams.thetaLowerRes);
            angleMap.phi_pdfsums_higherTheta.resize(sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (phi CDF line sums)");
        }
        try {
            angleMap.theta_CDF_lower.resize(sh.anglemapParams.thetaLowerRes);
            angleMap.theta_CDF_higher.resize(sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (line sums, CDF)");
        }
        try {
            angleMap.phi_CDFs_lowerTheta.resize(sh.anglemapParams.phiWidth * sh.anglemapParams.thetaLowerRes);
            angleMap.phi_CDFs_higherTheta.resize(sh.anglemapParams.phiWidth * sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (CDF)");
        }
        try {
            angleMap.phi_pdfs_lowerTheta.resize(sh.anglemapParams.phiWidth * sh.anglemapParams.thetaLowerRes);
            angleMap.phi_pdfs_higherTheta.resize(sh.anglemapParams.phiWidth * sh.anglemapParams.thetaHigherRes);
        }
        catch (...) {
            throw std::runtime_error("Not enough memory to load incident angle map (CDF)");
        }

        //First pass: determine phi line sums and total map sum
        //Lower part
        angleMap.theta_CDFsum_lower = 0;
        memset(angleMap.phi_pdfsums_lowerTheta.data(), 0, sizeof(size_t) * sh.anglemapParams.thetaLowerRes);
        for (size_t thetaIndex = 0; thetaIndex < sh.anglemapParams.thetaLowerRes; thetaIndex++) {
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                angleMap.phi_pdfsums_lowerTheta[thetaIndex] += angleMap.pdf[thetaIndex*sh.anglemapParams.phiWidth + phiIndex]; //phi line sum
            }
            angleMap.theta_CDFsum_lower += angleMap.phi_pdfsums_lowerTheta[thetaIndex]; //total lower map sum
        }
        //Higher part
        angleMap.theta_CDFsum_higher = angleMap.theta_CDFsum_lower; //higher includes lower part
        memset(angleMap.phi_pdfsums_higherTheta.data(), 0, sizeof(size_t) * sh.anglemapParams.thetaHigherRes);
        for (size_t thetaIndex = sh.anglemapParams.thetaLowerRes; thetaIndex < (sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes); thetaIndex++) {
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                angleMap.phi_pdfsums_higherTheta[thetaIndex-sh.anglemapParams.thetaLowerRes] += angleMap.pdf[thetaIndex*sh.anglemapParams.phiWidth + phiIndex]; //phi line sum
            }
            angleMap.theta_CDFsum_higher += angleMap.phi_pdfsums_higherTheta[thetaIndex - sh.anglemapParams.thetaLowerRes]; //total map sum
        }
        if (angleMap.theta_CDFsum_higher == 0) {
            auto err = fmt::format("Facet {} has all-zero recorded angle map, but is being used for desorption.", globalId + 1);
            throw std::runtime_error(err.c_str());
        }
        angleMap.thetaLowerRatio = (double)angleMap.theta_CDFsum_lower / (double)angleMap.theta_CDFsum_higher; //higher includes lower

        //Second pass: construct midpoint CDFs, normalized phi pdfs
        //We use midpoint because user expects linear interpolation between PDF midpoints, not between bin endpoints
        double thetaNormalizingFactor = 1.0 / (double)angleMap.theta_CDFsum_higher;
        //lower part
        for (size_t thetaIndex = 0; thetaIndex < sh.anglemapParams.thetaLowerRes; thetaIndex++) {
            if (thetaIndex == 0) {
                //First CDF value: sums from theta=0 to midpoint of first theta bin
                angleMap.theta_CDF_lower[thetaIndex] = 0.5 * (double)angleMap.phi_pdfsums_lowerTheta[0] * thetaNormalizingFactor;
            }
            else {
                //value summing second half of previous theta bin and first half of current theta bin
                angleMap.theta_CDF_lower[thetaIndex] = angleMap.theta_CDF_lower[thetaIndex - 1] + (double)(angleMap.phi_pdfsums_lowerTheta[thetaIndex - 1] + angleMap.phi_pdfsums_lowerTheta[thetaIndex])*0.5*thetaNormalizingFactor;
            }
            double phiNormalizingFactor = 1.0 / (double)angleMap.phi_pdfsums_lowerTheta[thetaIndex];
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                size_t index = sh.anglemapParams.phiWidth * thetaIndex + phiIndex;
                if (angleMap.phi_pdfsums_lowerTheta[thetaIndex] == 0) { //no hits in this line, create phi CDF of uniform distr.
                    angleMap.phi_CDFs_lowerTheta[index] = (0.5 + (double)phiIndex) / (double)sh.anglemapParams.phiWidth;
                    angleMap.phi_pdfs_lowerTheta[index] = 1.0 / (double)sh.anglemapParams.phiWidth;
                }
                else {  //construct regular CDF based on midpoints
                    if (phiIndex == 0) {
                        //First CDF value, sums to midpoint of first phi bin
                        angleMap.phi_CDFs_lowerTheta[index] = 0.5 * (double)angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex] * phiNormalizingFactor; //sh.anglemapParams.phiWidth * thetaIndex is the first pdf bin of this phi line
                    }
                    else {
                        //value covering second half of last phi bin and first of current phi bin
                        angleMap.phi_CDFs_lowerTheta[index] = angleMap.phi_CDFs_lowerTheta[sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + (double)(angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex])*0.5*phiNormalizingFactor;
                    }
                    angleMap.phi_pdfs_lowerTheta[index] = angleMap.pdf[index] * phiNormalizingFactor;
                }
            }
        }
        //second pass, higher part, pay attention to index shift: thetaIndex upshifted by thetaLowerRes
        for (size_t thetaIndex = sh.anglemapParams.thetaLowerRes; thetaIndex < sh.anglemapParams.thetaLowerRes + sh.anglemapParams.thetaHigherRes; thetaIndex++) {
            if (thetaIndex == sh.anglemapParams.thetaLowerRes) {
                //First CDF value: sums from theta=limit to midpoint of first higher theta bin
                angleMap.theta_CDF_higher[thetaIndex-sh.anglemapParams.thetaLowerRes] = angleMap.thetaLowerRatio + 0.5 * (double)angleMap.phi_pdfsums_higherTheta[0] * thetaNormalizingFactor;
            }
            else {
                //value summing second half of previous theta bin and first half of current theta bin
                angleMap.theta_CDF_higher[thetaIndex-sh.anglemapParams.thetaLowerRes] = angleMap.theta_CDF_higher[thetaIndex - sh.anglemapParams.thetaLowerRes - 1] + (double)(angleMap.phi_pdfsums_higherTheta[thetaIndex - sh.anglemapParams.thetaLowerRes - 1] + angleMap.phi_pdfsums_higherTheta[thetaIndex - sh.anglemapParams.thetaLowerRes])*0.5*thetaNormalizingFactor;
            }
            double phiNormalizingFactor = 1.0 / (double)angleMap.phi_pdfsums_higherTheta[thetaIndex - sh.anglemapParams.thetaLowerRes];
            for (size_t phiIndex = 0; phiIndex < sh.anglemapParams.phiWidth; phiIndex++) {
                size_t index = sh.anglemapParams.phiWidth * (thetaIndex - sh.anglemapParams.thetaLowerRes) + phiIndex;
                if (angleMap.phi_pdfsums_higherTheta[thetaIndex - sh.anglemapParams.thetaLowerRes] == 0) { //no hits in this line, create phi CDF of uniform distr.
                    angleMap.phi_CDFs_higherTheta[index] = (0.5 + (double)phiIndex) / (double)sh.anglemapParams.phiWidth;
                    angleMap.phi_pdfs_higherTheta[index] = 1.0 / (double)sh.anglemapParams.phiWidth;
                }
                else {
                    if (phiIndex == 0) {
                        //First CDF value, sums to midpoint of first phi bin
                        angleMap.phi_CDFs_higherTheta[index] = 0.5 * (double)angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex] * phiNormalizingFactor; //sh.anglemapParams.phiWidth * thetaIndex is the first pdf bin of this phi line
                    }
                    else {
                        //value covering second half of last phi bin and first of current phi bin
                        angleMap.phi_CDFs_higherTheta[index] = angleMap.phi_CDFs_higherTheta[sh.anglemapParams.phiWidth * (thetaIndex - sh.anglemapParams.thetaLowerRes) + phiIndex - 1] + (double)(angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex - 1] + angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex])*0.5*phiNormalizingFactor;
                    }
                    angleMap.phi_pdfs_higherTheta[index] = angleMap.pdf[sh.anglemapParams.phiWidth * thetaIndex + phiIndex] * phiNormalizingFactor;
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

void MolflowSimFacet::InitializeOutgassingMap()
{
    if (sh.useOutgassingFile) {
        //Precalc actual outgassing map width and height for faster generation:
        ogMap.outgassingMapWidth_precise = sh.U.Length() * ogMap.outgassingFileRatioU;
        ogMap.outgassingMapHeight_precise = sh.V.Length() * ogMap.outgassingFileRatioV;
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

size_t MolflowSimFacet::InitializeHistogram(const size_t &nbMoments) const
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

bool MolflowSimFacet::InitializeLinkAndVolatile(const size_t & id)
{
    SimulationFacet::InitializeLinkAndVolatile(id);
    if (sh.superDest || sh.isVolatile) {
        sh.opacity_paramId = -1;
        sh.sticking_paramId = -1;
    }
    return true;
}

/**
* \brief Calculates the hits size for a single facet
* \param nbMoments amount of moments (TD simulations are currently not used for Molflow)
* \return calculated size of the facet hits
*/
size_t MolflowSimFacet::GetHitsSize(size_t nbMoments) const {
    // Init with base size
    size_t hits_size = SimulationFacet::GetHitsSize(nbMoments);

    return hits_size
    + (sh.anglemapParams.record ? (sh.anglemapParams.GetRecordedDataSize()) : 0);
}

size_t MolflowSimFacet::GetMemSize() const {

    // Init with base size
    size_t mem_size = SimulationFacet::GetMemSize();

    mem_size += sizeof (double) * ogMap.outgassingMap.capacity();
    mem_size += angleMap.GetMemSize();
    return mem_size;
}