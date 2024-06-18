#include <cstdio>
#include <tuple> //std::tie
#include <cstring>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include "MolflowSimulation.h"
#include "MolflowSimFacet.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "Helper/MathTools.h"

#include "Parameter.h"
#include <omp.h>

double MolflowSimulationModel::GetStickingAt(const MolflowSimFacet *f, const double time) const {
    if (f->sticking_paramId == -1) //constant sticking
        return f->sh.sticking;
    else {
        auto &par = tdParams.parameters[f->sticking_paramId];
        return InterpolateY(time, par.GetValues(), par.logXinterp, par.logYinterp, false);
    }
}

double MolflowSimulationModel::GetOpacityAt(const MolflowSimFacet *f, const double time) const {
    if (f->opacity_paramId == -1) //constant opacity
        return f->sh.opacity;
    else {
        auto &par = tdParams.parameters[f->opacity_paramId];
        return InterpolateY(time, par.GetValues(), par.logXinterp, par.logYinterp, false);
    }
}

double MolflowSimulationModel::GetTemperatureAt(const MolflowSimFacet *f, const double time) const {
    if (f->temperature_paramId == -1) //constant temp
        return f->sh.temperature;
    else {
        auto &par = tdParams.parameters[f->temperature_paramId];
        return InterpolateY(time, par.GetValues(), par.logXinterp, par.logYinterp, false);
    }
}