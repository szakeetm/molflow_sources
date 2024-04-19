

#ifndef MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H
#define MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H

#include "FormulaEvaluator.h"
#include <Interface/GeometryViewer.h>

class Worker;
class MolflowGeometry;

class FormulaEvaluator_MF : public FormulaEvaluator {
public:
    FormulaEvaluator_MF(Worker* w, MolflowGeometry* interfGeom, std::vector<SelectionGroup>* sel);
    bool EvaluateVariable(std::list<Variable>::iterator v, const std::vector <std::pair<std::string, std::optional<double>>>& previousFormulaValues) override;

    Worker* worker;
    MolflowGeometry* geometry;
    std::vector<SelectionGroup>* selections;
};


#endif //MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H
