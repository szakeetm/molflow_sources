//
// Created by pbahr on 12/10/2020.
//

#ifndef MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H
#define MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H

#include "FormulaEvaluator.h"
#include <Interface/GeometryViewer.h>

class Worker;
class MolflowGeometry;

class FormulaEvaluator_MF : FormulaEvaluator {
public:
    FormulaEvaluator_MF(Worker* w, MolflowGeometry* geom, std::vector<SelectionGroup>* sel);
    bool EvaluateVariable(VLIST *v) override;

    Worker* worker;
    MolflowGeometry* geometry;
    std::vector<SelectionGroup>* selections;
};


#endif //MOLFLOW_PROJ_FORMULAEVALUATOR_MF_H
