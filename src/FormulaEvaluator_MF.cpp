//
// Created by pbahr on 12/10/2020.
//

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

bool FormulaEvaluator_MF::EvaluateVariable(VLIST *v) {
    bool ok = true;
    Geometry* geom = worker->GetGeometry();
    int nbFacet = geom->GetNbFacet();
    int idx;

    if ((idx = GetVariable(v->name, "A")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.nbAbsEquiv;
    }
    else if ((idx = GetVariable(v->name, "D")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed;
    }
    else if ((idx = GetVariable(v->name, "MCH")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbMCHit;
    }
    else if ((idx = GetVariable(v->name, "H")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = (double)geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv;
    }
    else if ((idx = GetVariable(v->name, "P")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.sum_v_ort *
                           worker->GetMoleculesPerTP(worker->displayedMoment)*1E4 / geom->GetFacet(idx - 1)->GetArea() * (worker->model.wp.gasMass / 1000 / 6E23)*0.0100;
    }
    else if ((idx = GetVariable(v->name, "DEN")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) {
            InterfaceFacet *f = geom->GetFacet(idx - 1);
            v->value = f->DensityCorrection() * f->facetHitCache.sum_1_per_ort_velocity /
                       f->GetArea() *
                       worker->GetMoleculesPerTP(worker->displayedMoment)*1E4;
        }
    }
    else if ((idx = GetVariable(v->name, "Z")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv /
                           geom->GetFacet(idx - 1)->GetArea() *
                           worker->GetMoleculesPerTP(worker->displayedMoment)*1E4;
    }
    else if ((idx = GetVariable(v->name, "V")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) /*v->value = 4.0*(double)(geom->GetFacet(idx - 1)->facetHitCache.nbMCHit + geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed) /
			geom->GetFacet(idx - 1)->facetHitCache.sum_1_per_ort_velocity;*/
            v->value = (geom->GetFacet(idx - 1)->facetHitCache.nbHitEquiv + static_cast<double>(geom->GetFacet(idx - 1)->facetHitCache.nbDesorbed)) / geom->GetFacet(idx - 1)->facetHitCache.sum_1_per_velocity;
    }
    else if ((idx = GetVariable(v->name, "T")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->sh.temperature;
    }
    else if ((idx = GetVariable(v->name, "AR")) > 0) {
        ok = (idx > 0 && idx <= nbFacet);
        if (ok) v->value = geom->GetFacet(idx - 1)->sh.area;
    }
    else if (iequals(v->name, "SUMDES")) {
        v->value = (double)worker->globalHitCache.globalHits.nbDesorbed;
    }
    else if (iequals(v->name, "SUMABS")) {
        v->value = worker->globalHitCache.globalHits.nbAbsEquiv;
    }
    else if (iequals(v->name, "SUMMCHIT")) {
        v->value = (double)worker->globalHitCache.globalHits.nbMCHit;
    }
    else if (iequals(v->name, "SUMHIT")) {
        v->value = worker->globalHitCache.globalHits.nbHitEquiv;
    }
    else if (iequals(v->name, "MPP")) {
        v->value = worker->globalHitCache.distTraveled_total / (double)worker->globalHitCache.globalHits.nbDesorbed;
    }
    else if (iequals(v->name, "MFP")) {
        v->value = worker->globalHitCache.distTraveledTotal_fullHitsOnly / worker->globalHitCache.globalHits.nbHitEquiv;
    }
    else if (iequals(v->name, "DESAR")) {
        double sumArea = 0.0;
        for (size_t i2 = 0; i2 < geom->GetNbFacet(); i2++) {
            InterfaceFacet *f_tmp = geom->GetFacet(i2);
            if (f_tmp->sh.desorbType) sumArea += f_tmp->GetArea();
        }
        v->value = sumArea;
    }
    else if (iequals(v->name, "ABSAR")) {
        double sumArea = 0.0;

        for (size_t i2 = 0; i2 < geom->GetNbFacet(); i2++) {
            InterfaceFacet *f_tmp = geom->GetFacet(i2);
            if (f_tmp->sh.sticking > 0.0) sumArea += f_tmp->GetArea()*f_tmp->sh.opacity;
        }
        v->value = sumArea;
    }
    else if (iequals(v->name, "QCONST")) {
        v->value = worker->model.wp.finalOutgassingRate_Pa_m3_sec*10.00; //10: Pa*m3/sec -> mbar*l/s
    }
    else if (iequals(v->name, "QCONST_N")) {
        v->value = worker->model.wp.finalOutgassingRate;
    }
    else if (iequals(v->name, "NTOT")) {
        v->value = worker->model.wp.totalDesorbedMolecules;
    }
    else if (iequals(v->name, "GASMASS")) {
        v->value = worker->model.wp.gasMass;
    }
    else if (iequals(v->name, "KB")) {
        v->value = 1.3806504e-23;
    }
    else if (iequals(v->name, "R")) {
        v->value = 8.314472;
    }
    else if (iequals(v->name, "Na")) {
        v->value = 6.02214179e23;
    }
    else if ((beginsWith(v->name, "SUM(") || beginsWith(v->name, "sum(") || beginsWith(v->name, "AVG(") || beginsWith(v->name, "avg(")) && endsWith(v->name, ")")) {
        bool avgMode = (beginsWith(v->name, "AVG(") || beginsWith(v->name, "avg(")); //else SUM mode
        std::string inside = v->name; inside.erase(0, 4); inside.erase(inside.size() - 1, 1);
        std::vector<std::string> tokens = SplitString(inside,',');
        if (!Contains({ 2,3 }, tokens.size()))
            return false;
        if (avgMode) {
            if (!Contains({ "P","DEN","Z","p","den","z" }, tokens[0]))
                return false;
        }
        else {
            if (!Contains({ "MCH","H","D","A","AR","mch","h","d","a","ar" }, tokens[0]))
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
            if (!(beginsWith(tokens[1], "S") || beginsWith(tokens[1], "s"))) return false;
            std::string selIdString = tokens[1]; selIdString.erase(0, 1);
            if (Contains({ "EL","el" }, selIdString)) { //Current selections
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
                facetsToSum = (*selections)[selGroupId - 1].selection;
            }
        }
        size_t sumLL=0;
        double sumD=0.0;
        double sumArea = 0.0; //We average by area
        for (auto& sel : facetsToSum) {
            if (Contains({"MCH", "mch"},tokens[0])) {
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
                       (worker->model.wp.gasMass / 1000 / 6E23)*0.0100;
                sumArea += geom->GetFacet(sel)->GetArea();
            } else if (Contains({ "DEN", "den" }, tokens[0])) {
                InterfaceFacet *f = geom->GetFacet(sel);
                sumD += f->DensityCorrection() * f->facetHitCache.sum_1_per_ort_velocity;
                sumArea += geom->GetFacet(sel)->GetArea();
            } else if (Contains({ "Z", "z" }, tokens[0])) {
                sumD += geom->GetFacet(sel)->facetHitCache.nbHitEquiv;
                sumArea += geom->GetFacet(sel)->GetArea();
            } else return false;
        }
        if (avgMode) v->value=sumD * worker->GetMoleculesPerTP(worker->displayedMoment)*1E4 / sumArea;
        else if (Contains({ "AR", "ar" }, tokens[0])) v->value = sumArea;
        else if (Contains({ "H", "h", "A", "a" }, tokens[0])) v->value = sumD;
        else v->value = static_cast<double>(sumLL); //Only one conversion at the end (instead of at each summing operation)
    }
    else ok = false;
    return ok;
}