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

#ifndef MOLFLOW_PROJ_MOLFLOWSIMFACET_H
#define MOLFLOW_PROJ_MOLFLOWSIMFACET_H

#include "SimulationFacet.h"
#include "MolflowTypes.h"

struct Anglemap {
public:
    Anglemap(){
        theta_CDFsum_lower=theta_CDFsum_higher=0;
    };
    std::vector<size_t> pdf;          // Incident angle distribution, large array of phi and theta, not normalized by solid angle or to 1 (simply number of hits in each bin). Used either for recording or for 2nd order interpolation. Unites lower and higher theta parts
    std::vector<double> phi_CDFs_lowerTheta;    // A table containing cumulative phi distributions, normalized to 1, summed up to the phi bin midpoint (!), for each theta, starting from 0 for every line (1 line = 1 theta bin). For speed, one big array of (theta_lowerRes*phi_size) instead of vector of vectors
    std::vector<double> phi_CDFs_higherTheta;
    std::vector<size_t> phi_pdfsums_lowerTheta; // since phi_CDFs sums only to the middle of the last phi bin, for each theta a line sum (of the raw i.e. not normalized angle map pdf) is stored here. Also a pdf for theta, as it contains the total possibility, summed over all phi angles, for that theta bin.
    std::vector<size_t> phi_pdfsums_higherTheta;
    std::vector<double> phi_pdfs_lowerTheta; //Normalized pdf of each phi line. For speed, one big array of (theta_lowerRes*phi_size) instead of vector of vectors
    std::vector<double> phi_pdfs_higherTheta;
    std::vector<double> theta_CDF_lower;      // Theta CDF to each theta bin midpoint, normalized to 1. nth value is the CDF at the midpoint of theta bin n.)
    std::vector<double> theta_CDF_higher;
    size_t theta_CDFsum_lower;  // since theta CDF only sums till the midpoint of the last segment, the total map sum is here
    size_t theta_CDFsum_higher; // theta_CDFsum_higher>=theta_CDFsum_lower as it includes lower part. Also total number of hits in raw pdf
    double thetaLowerRatio; // ratio of angle map below theta limit, to decide which side to look up in

    [[nodiscard]] size_t GetMemSize() const {
        size_t sum = 0;
        sum += sizeof(Anglemap);
        sum += sizeof(size_t) * pdf.capacity();
        sum += sizeof(double) * phi_CDFs_lowerTheta.capacity();
        sum += sizeof(double) * phi_CDFs_higherTheta.capacity();
        sum += sizeof(double) * phi_pdfs_lowerTheta.capacity();
        sum += sizeof(double) * phi_pdfs_higherTheta.capacity();
        sum += sizeof(size_t) * phi_pdfsums_lowerTheta.capacity();
        sum += sizeof(size_t) * phi_pdfsums_higherTheta.capacity();
        sum += sizeof(double) * theta_CDF_lower.capacity();
        sum += sizeof(double) * theta_CDF_higher.capacity();
        return sum;
    }
};

struct MolflowSimFacet : public SimulationFacet {
    MolflowSimFacet() : SimulationFacet() {};
    explicit MolflowSimFacet(size_t nbIndex) : SimulationFacet(nbIndex) {};
    MolflowSimFacet(const MolflowSimFacet &o);
    MolflowSimFacet(MolflowSimFacet &&cpy) noexcept;
    MolflowSimFacet &operator=(const MolflowSimFacet &o);
    MolflowSimFacet &operator=(MolflowSimFacet &&o) noexcept;

    bool InitializeOnLoad(const size_t &id, const size_t &nbMoments);

    size_t InitializeHistogram(const size_t &nbMoments) const;

    size_t InitializeDirectionTexture(const size_t &nbMoments){
        return SimulationFacet::InitializeDirectionTexture();
    };

    size_t InitializeProfile(const size_t &nbMoments){
        return SimulationFacet::InitializeProfile();
    };

    size_t InitializeTexture(const size_t &nbMoments) {
        return SimulationFacet::InitializeTexture();
    };

    bool InitializeLinkAndVolatile(const size_t & id) override;

    int InitializeAngleMap();

    void InitializeOutgassingMap();

    [[nodiscard]] size_t GetHitsSize(size_t nbMoments) const override;
    size_t GetMemSize() const override;

    OutgassingMap ogMap;
    Anglemap angleMap;
};

#endif //MOLFLOW_PROJ_MOLFLOWSIMFACET_H