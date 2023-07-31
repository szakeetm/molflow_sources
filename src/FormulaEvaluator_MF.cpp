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

#include <Helper/StringHelper.h>
#include <Helper/MathTools.h>
#include <numeric>
#include "FormulaEvaluator_MF.h"
#include "Worker.h"
#include "MolflowGeometry.h"
#include "Facet_shared.h"

FormulaEvaluator_MF::FormulaEvaluator_MF(Worker* w, MolflowGeometry* geom, std::vector<SelectionGroup>* sel){
    worker = w;
    geometry = geom;
    selections = sel;
}

bool FormulaEvaluator_MF::EvaluateVariable(std::list<Variable>::iterator v, const std::vector <std::pair<std::string, std::optional<double>>>& aboveFormulaValues) {
    bool ok = true;
    Geometry* geom = worker->GetGeometry();
    int nbFacet = geom->GetNbFacet();
    int idx;


    if ((idx = GetFacetIndex(v->varName, "A")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.nbAbsEquiv;
    }
    else if ((idx = GetFacetIndex(v->varName, "D")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed;
    }
    else if ((idx = GetFacetIndex(v->varName, "MCH")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbMCHit;
    }
    else if ((idx = GetFacetIndex(v->varName, "H")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv;
    }
    else if ((idx = GetFacetIndex(v->varName, "P")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.sum_v_ort *
                           worker->GetMoleculesPerTP(worker->displayedMoment)*1E4 / geom->GetFacet(idx - 1)->GetArea() * (worker->model->wp.gasMass / 1000 / 6E23)*0.0100;
    }
    else if ((idx = GetFacetIndex(v->varName, "DEN")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) {
            InterfaceFacet *f = geom->GetFacet(idx - 1);
            v->value = f->DensityCorrection() * f->facetHitCache.sum_1_per_ort_velocity /
                       f->GetArea() *
                       worker->GetMoleculesPerTP(worker->displayedMoment)*1E4;
        }
    }
    else if ((idx = GetFacetIndex(v->varName, "Z")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv /
                           geom->GetFacet(idx - 1)->GetArea() *
                           worker->GetMoleculesPerTP(worker->displayedMoment)*1E4;
    }
    else if ((idx = GetFacetIndex(v->varName, "V")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) /*v->value = 4.0*(double)(geom->GetFacet(idx - 1)->facetHitCache.nbMCHit + geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed) /
			geom->GetFacet(idx - 1)->facetHitCache.sum_1_per_ort_velocity;*/
            v->value = (geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv + static_cast<double>(geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed)) / geom->GetFacet(idx - 1)->facetHitCache.sum_1_per_velocity;
    }
    else if ((idx = GetFacetIndex(v->varName, "T")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->sh.temperature;
    }
    else if ((idx = GetFacetIndex(v->varName, "AR")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->sh.area;
    }
    else if ((idx = GetFacetIndex(v->varName, "Force")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) {
            InterfaceFacet* f = geom->GetFacet(idx - 1);
            auto forceN = f->facetHitCache.impulse.Norme() * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
            v->value = forceN;
        }
    }
    else if ((idx = GetFacetIndex(v->varName, "ForceX")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) {
            InterfaceFacet* f = geom->GetFacet(idx - 1);
            auto forceX = f->facetHitCache.impulse.x * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
            v->value = forceX;
        }
    }
	else if ((idx = GetFacetIndex(v->varName, "ForceY")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto forceY = f->facetHitCache.impulse.y * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
			v->value = forceY;
		}
	}
	else if ((idx = GetFacetIndex(v->varName, "ForceZ")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto forceZ = f->facetHitCache.impulse.z * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
			v->value = forceZ;
		}
	}
    else if ((idx = GetFacetIndex(v->varName, "ForceSqr")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) {
            InterfaceFacet* f = geom->GetFacet(idx - 1);
            auto force_sqrN = f->facetHitCache.impulse_square.Norme() * worker->GetMoleculesPerTP(worker->displayedMoment)
                * Square(worker->model->wp.gasMass / 1000 / 6E23);
            v->value = force_sqrN;
            if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
        }
    }
    else if ((idx = GetFacetIndex(v->varName, "ForceSqrX")) > 0) {
        ok = (idx <= nbFacet);
        if (ok) {
            InterfaceFacet* f = geom->GetFacet(idx - 1);
            auto force_sqrX = f->facetHitCache.impulse_square.x * worker->GetMoleculesPerTP(worker->displayedMoment) * Square(worker->model->wp.gasMass / 1000 / 6E23);
            v->value = force_sqrX;
            if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
        }
    }
	else if ((idx = GetFacetIndex(v->varName, "ForceSqrY")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto force_sqrY = f->facetHitCache.impulse_square.y * worker->GetMoleculesPerTP(worker->displayedMoment) * Square(worker->model->wp.gasMass / 1000 / 6E23);
			v->value = force_sqrY;
            if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
		}
	}
	else if ((idx = GetFacetIndex(v->varName, "ForceSqrZ")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto force_sqrZ = f->facetHitCache.impulse_square.z * worker->GetMoleculesPerTP(worker->displayedMoment) * Square(worker->model->wp.gasMass / 1000 / 6E23);
			v->value = force_sqrZ;
            if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
		}
	}
	else if ((idx = GetFacetIndex(v->varName, "Torque")) > 0) {
	    ok = (idx <= nbFacet);
	    if (ok) {
		    InterfaceFacet* f = geom->GetFacet(idx - 1);
		    auto torqueN = f->facetHitCache.impulse_momentum.Norme() * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //0.01: N*cm to Nm
		    v->value = torqueN;
	    }
	}
    else if ((idx = GetFacetIndex(v->varName, "TorqueX")) > 0) {
    ok = (idx <= nbFacet);
    if (ok) {
        InterfaceFacet* f = geom->GetFacet(idx - 1);
        auto torqueX = f->facetHitCache.impulse_momentum.x * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //0.01: N*cm to Nm
        v->value = torqueX;
    }
    }
	else if ((idx = GetFacetIndex(v->varName, "TorqueY")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto torqueY = f->facetHitCache.impulse_momentum.y * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //0.01: N*cm to Nm
			v->value = torqueY;
		}
	}
	else if ((idx = GetFacetIndex(v->varName, "TorqueZ")) > 0) {
		ok = (idx <= nbFacet);
		if (ok) {
			InterfaceFacet* f = geom->GetFacet(idx - 1);
			auto torqueZ = f->facetHitCache.impulse_momentum.z * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //0.01: N*cm to Nm
			v->value = torqueZ;
		}
	}
	else if (iequals(v->varName, "SUMDES")) {
        v->value = (double)worker->globalStatCache.globalHits.nbDesorbed;
    }
    else if (iequals(v->varName, "SUMABS")) {
        v->value = worker->globalStatCache.globalHits.nbAbsEquiv;
    }
    else if (iequals(v->varName, "SUMMCHIT")) {
        v->value = (double)worker->globalStatCache.globalHits.nbMCHit;
    }
    else if (iequals(v->varName, "SUMHIT")) {
        v->value = worker->globalStatCache.globalHits.nbHitEquiv;
    }
    else if (iequals(v->varName, "MPP")) {
        v->value = worker->globalStatCache.distTraveled_total / (double)worker->globalStatCache.globalHits.nbDesorbed;
    }
    else if (iequals(v->varName, "MFP")) {
        v->value = worker->globalStatCache.distTraveledTotal_fullHitsOnly / worker->globalStatCache.globalHits.nbHitEquiv;
    }
    else if (iequals(v->varName, "DESAR")) {
        double sumArea = 0.0;
        for (size_t i2 = 0; i2 < geom->GetNbFacet(); i2++) {
            InterfaceFacet *f_tmp = geom->GetFacet(i2);
            if (f_tmp->sh.desorbType) sumArea += f_tmp->GetArea();
        }
        v->value = sumArea;
    }
    else if (iequals(v->varName, "ABSAR")) {
        double sumArea = 0.0;

        for (size_t i2 = 0; i2 < geom->GetNbFacet(); i2++) {
            InterfaceFacet *f_tmp = geom->GetFacet(i2);
            if (f_tmp->sh.sticking > 0.0) sumArea += f_tmp->GetArea()*f_tmp->sh.opacity;
        }
        v->value = sumArea;
    }
    else if (iequals(v->varName, "QCONST")) {
        v->value = worker->model->wp.finalOutgassingRate_Pa_m3_sec*10.00; //10: Pa*m3/sec -> mbar*l/s
    }
    else if (iequals(v->varName, "QCONST_N")) {
        v->value = worker->model->wp.finalOutgassingRate;
    }
    else if (iequals(v->varName, "NTOT")) {
        v->value = worker->model->wp.totalDesorbedMolecules;
    }
    else if (iequals(v->varName, "GASMASS")) {
        v->value = worker->model->wp.gasMass;
    }
    else if (iequals(v->varName, "KB")) {
        v->value = 1.3806504e-23;
    }
    else if (iequals(v->varName, "R")) {
        v->value = 8.314472;
    }
    else if (iequals(v->varName, "Na")) {
        v->value = 6.02214179e23;
    }
    else if ((idx = GetFacetIndex(v->varName,"Formula"))>0) { //Refer to previously evaluated formula
        if (idx > aboveFormulaValues.size()) {
            throw Error(fmt::format("Formula {} should be defined above to refer to its value", idx));
        }
        if (!aboveFormulaValues[idx - 1].second.has_value()) {
            throw Error(fmt::format("Formula {} is not yet evaluated", idx));
        }
        ok = true;
        v->value = aboveFormulaValues[idx - 1].second.value();
    }
    else if (v->varName.length() >= 2 && v->varName[0] == '"' && v->varName[v->varName.length() - 1] == '"') { //"formula name"
        bool found = false;
        auto it = aboveFormulaValues.begin();
        while (!found && it != aboveFormulaValues.end()) {
            found = ('"' + it->first + '"') == v->varName;
            if (!found) {
                ++it;
            }
        }
        if (!found) throw Error(fmt::format("Formula name {} not found above", v->varName));
        if (!it->second.has_value()) throw Error(fmt::format("Formula {} is not yet evaluated", v->varName));
        ok = true;
        v->value = it->second.value();
    }
    else if ((beginsWith(uppercase(v->varName), "SUM(") || beginsWith(uppercase(v->varName), "AVG(")) && endsWith(v->varName, ")")) {
        bool avgMode = (beginsWith(uppercase(v->varName), "AVG(")); //else SUM mode
        std::string inside = v->varName; inside.erase(0, 4); inside.erase(inside.size() - 1, 1);
        std::vector<std::string> tokens = SplitString(inside,',');
        if (!Contains({ 2,3 }, tokens.size()))
            return false;
        if (avgMode) {
            if (!iContains({ "P","DEN","Z","Force","ForceX","ForceY","ForceZ","ForceSqr","ForceSqrX","ForceSqrY","ForceSqrZ","Torque","TorqueX","TorqueY","TorqueZ"}, tokens[0]))
                return false;
        }
        else {
            if (!iContains({ "MCH","H","D","A","AR","Force","ForceX","ForceY","ForceZ","ForceSqr","ForceSqrX","ForceSqrY","ForceSqrZ","Torque","TorqueX","TorqueY","TorqueZ" }, tokens[0]))
                return false;
        }
        std::vector<size_t> facetsToSum;
        if (tokens.size() == 3) { // Like SUM(H,3,6) = H3 + H4 + H5 + H6
            size_t startId, endId, pos;
            try {
                startId = std::stol(tokens[1], &pos); if (pos != tokens[1].size() || startId > geom->GetNbFacet() || startId == 0) return false;
                endId = std::stol(tokens[2], &pos); if (pos != tokens[2].size() || endId > geom->GetNbFacet() || endId == 0) return false;
            }
            catch (...) {
                return false;
            }
            if (!(startId < endId)) return false;
            facetsToSum = std::vector<size_t>(endId-startId+1);
            std::iota(facetsToSum.begin(), facetsToSum.end(), startId-1);
        }
        else { //Selection group
            if (!(beginsWith(uppercase(tokens[1]), "S"))) return false;
            std::string selIdString = tokens[1]; selIdString.erase(0, 1);
            if (iContains({ "EL" }, selIdString)) { //Current selections
                facetsToSum = geom->GetSelectedFacets();
            }
            else {
                size_t selGroupId, pos;
                try {
                    selGroupId = std::stol(selIdString, &pos); if (pos != selIdString.size() || selGroupId > selections->size() || selGroupId == 0) return false;
                }
                catch (...) {
                    return false;
                }
                facetsToSum = (*selections)[selGroupId - 1].facetIds;
            }
        }
        size_t sumLL=0;
        double sumD=0.0;
        double sumArea = 0.0; //We average by area
        for (auto& sel : facetsToSum) {
            if (iequals("MCH",tokens[0])) {
                sumLL+=geom->GetFacet(sel)->facetHitCache.nbMCHit;
            }
            else if (Contains({ "H", "h" }, tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.nbHitEquiv;
            }
            else if (Contains({ "D", "d" }, tokens[0])) {
                sumLL+=geom->GetFacet(sel)->facetHitCache.nbDesorbed;
            } else if (Contains({ "A", "a" }, tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.nbAbsEquiv;
            } else if (Contains({ "AR", "ar" }, tokens[0])) {
                sumArea += geom->GetFacet(sel)->GetArea();
            }
            else if (Contains({ "P", "p" }, tokens[0])) {
                sumD+= geom->GetFacet(sel)->facetHitCache.sum_v_ort *
                       (worker->model->wp.gasMass / 1000 / 6E23)*0.0100;
                sumArea += geom->GetFacet(sel)->GetArea();
            } else if (Contains({ "DEN", "den" }, tokens[0])) {
                InterfaceFacet *f = geom->GetFacet(sel);
                sumD += f->DensityCorrection() * f->facetHitCache.sum_1_per_ort_velocity;
                sumArea += geom->GetFacet(sel)->GetArea();
            } else if (Contains({ "Z", "z" }, tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.nbHitEquiv;
                sumArea += geom->GetFacet(sel)->GetArea();
            }
            else if(iequals("Force", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse.Norme();
            }
            else if(iequals("ForceX", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse.x;
            }
            else if (iequals("ForceY", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse.y;
            }
            else if (iequals("ForceZ", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse.z;
            }
            else if (iequals("ForceSqr", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_square.Norme();
            }
            else if (iequals("ForceSqrX", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_square.x;
            }
            else if (iequals("ForceSqrY", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_square.y;
            }
            else if (iequals("ForceSqrZ", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_square.z;
            }
            else if (iequals("Torque", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_momentum.Norme();
            }
            else if (iequals("TorqueX", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_momentum.x;
            }
            else if (iequals("TorqueY", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_momentum.y;
            }
            else if (iequals("TorqueZ", tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.impulse_momentum.z;
            }
            else return false;
        }
        if (avgMode) {
            if (iContains({ "P","DEN","Z" },tokens[0])) v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * 1E4 / sumArea;
            else if (iContains({ "Force","ForceX","ForceY","ForceZ" }, tokens[0])) v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) / (double)facetsToSum.size();
            else if (iContains({ "ForceSqr","ForceSqrX","ForceSqrY","ForceSqrZ" }, tokens[0])) {
                v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * Square(worker->model->wp.gasMass / 1000 / 6E23) / (double)facetsToSum.size();
                if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
            }
            else if (iContains({ "Torque","TorqueX", "TorqueY", "TorqueZ" }, tokens[0])) v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01 / (double)facetsToSum.size(); //Ncm to Nm
        }
        else { //sum mode
            if (iequals("AR" , tokens[0])) v->value = sumArea;
            else if (iContains({ "H", "A" }, tokens[0])) v->value = sumD;
            else if (iContains({ "Force","ForceX","ForceY","ForceZ" }, tokens[0])) v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23);
            else if (iContains({ "ForceSqr","ForceSqrX","ForceSqrY","ForceSqrZ" }, tokens[0])){
                v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * Square(worker->model->wp.gasMass / 1000 / 6E23);
                if (worker->displayedMoment!=0) v->value/=worker->interfaceMomentCache[worker->displayedMoment - 1].window; //force2 divided by dt^2 to get N^2
            }
            else if (iContains({ "Torque","TorqueX", "TorqueY", "TorqueZ" }, tokens[0])) v->value = sumD * worker->GetMoleculesPerTP(worker->displayedMoment) * (worker->model->wp.gasMass / 1000 / 6E23) * 0.01; //Ncm to Nm
            else v->value = static_cast<double>(sumLL); //One long->double conversion at the end (instead of at each summing operation)
        }
     }
    else ok = false;
    return ok;
}