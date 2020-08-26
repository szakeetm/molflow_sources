//
// Created by Pascal Baehr on 22.04.20.
//

#include "MolflowBuffer.h"

void FacetHistogramBuffer::Resize(const HistogramParams& params){
    this->nbHitsHistogram.resize(params.recordBounce ? params.GetBounceHistogramSize() : 0); 
    this->nbHitsHistogram.shrink_to_fit();
    this->distanceHistogram.resize(params.recordDistance ? params.GetDistanceHistogramSize() : 0); 
    this->distanceHistogram.shrink_to_fit();
    this->timeHistogram.resize(params.recordTime ? params.GetTimeHistogramSize() : 0); 
    this->timeHistogram.shrink_to_fit();
}

void FacetHistogramBuffer::Reset(){
    ZEROVECTOR(nbHitsHistogram);
    ZEROVECTOR(distanceHistogram);
    ZEROVECTOR(timeHistogram);
}