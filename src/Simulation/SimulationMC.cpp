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
#include <cmath>
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