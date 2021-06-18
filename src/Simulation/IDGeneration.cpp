//
// Created by pascal on 2/9/21.
//

#include <Helper/MathTools.h>
#include "IDGeneration.h"

namespace IDGeneration {

    /**
* \brief Get ID (if it exists) of the integrated desorption (ID) function for a particular paramId
* \param paramId parameter ID
* \return Id of the integrated desorption function
*/
    int GetIDId(const std::set<size_t> &desorptionParameterIDs, int paramId) {
        if(!desorptionParameterIDs.empty()) {
            auto lowerBound = std::lower_bound(desorptionParameterIDs.begin(), desorptionParameterIDs.end(), paramId);
            if(lowerBound == desorptionParameterIDs.begin())
                return -1;
            --lowerBound; //even temperatureList.end() can be a bound

            if (paramId == *lowerBound) {
                return std::distance(desorptionParameterIDs.begin(), lowerBound);
            }
        }
        return -1;
    }

/**
* \brief Generate a new ID (integrated desorption) for desorption parameter for time-dependent simulations
* \param paramId parameter ID
* \return Previous size of IDs vector, which determines new id in the vector
*/
    std::pair<int, std::vector<ID_p>> GenerateNewID( std::set<size_t>& desorptionParameterIDs, int paramId, SimulationModel* model) {
        size_t i = desorptionParameterIDs.size();
        desorptionParameterIDs.insert(paramId);

        std::vector<ID_p> id_v = Generate_ID(paramId,model);
        return std::make_pair((int) i, id_v);
    }

/**
* \brief Generate integrated desorption (ID) function
* \param paramId parameter identifier
* \return ID as a Vector containing a pair of double values (x value = moment, y value = desorption value)
*/
    std::vector<std::pair<double, double>> Generate_ID(int paramId, SimulationModel *model) {
        std::vector<std::pair<double, double>> newID;
        //First, let's check at which index is the latest moment
        size_t indexBeforeLastMoment;
        for (indexBeforeLastMoment = 0; indexBeforeLastMoment < model->tdParams.parameters[paramId].GetSize() &&
                                        (model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment) <
                                         model->wp.latestMoment); indexBeforeLastMoment++);
        if (indexBeforeLastMoment >= model->tdParams.parameters[paramId].GetSize())
            indexBeforeLastMoment = model->tdParams.parameters[paramId].GetSize() - 1; //not found, set as last moment

//Construct integral from 0 to latest moment
//Zero
        newID.emplace_back(0.0, 0.0);

        //First moment
        newID.emplace_back(model->tdParams.parameters[paramId].GetX(0),
                           model->tdParams.parameters[paramId].GetX(0) * model->tdParams.parameters[paramId].GetY(0) *
                           0.100); //for the first moment (0.1: mbar*l/s -> Pa*m3/s)

        //Intermediate moments
        for (size_t pos = 1; pos <= indexBeforeLastMoment; pos++) {
            if (IsEqual(model->tdParams.parameters[paramId].GetY(pos),
                        model->tdParams.parameters[paramId].GetY(pos - 1))) //two equal values follow, simple integration by multiplying
                newID.emplace_back(model->tdParams.parameters[paramId].GetX(pos),
                                   newID.back().second +
                                   (model->tdParams.parameters[paramId].GetX(pos) - model->tdParams.parameters[paramId].GetX(pos - 1)) *
                                   model->tdParams.parameters[paramId].GetY(pos) * 0.100);
            else { //difficult case, we'll integrate by dividing to 20 equal sections
                for (double delta = 0.05; delta < 1.0001; delta += 0.05) {
                    double delta_t = model->tdParams.parameters[paramId].GetX(pos) - model->tdParams.parameters[paramId].GetX(pos - 1);
                    double time = model->tdParams.parameters[paramId].GetX(pos - 1) + delta * delta_t;
                    double avg_value = (model->tdParams.parameters[paramId].InterpolateY(time - 0.05 * delta_t, false) +
                                        model->tdParams.parameters[paramId].InterpolateY(time, false)) * 0.100 / 2.0;
                    newID.emplace_back(time, newID.back().second + 0.05 * delta_t * avg_value);
                }
            }
        }

        //wp.latestMoment
        double valueAtLatestMoment = model->tdParams.parameters[paramId].InterpolateY(model->wp.latestMoment, false);
        if (IsEqual(valueAtLatestMoment, model->tdParams.parameters[paramId].GetY(
                indexBeforeLastMoment))) //two equal values follow, simple integration by multiplying
            newID.emplace_back(model->wp.latestMoment,
                               newID.back().second +
                               (model->wp.latestMoment - model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment)) *
                               model->tdParams.parameters[paramId].GetY(indexBeforeLastMoment) * 0.100);
        else { //difficult case, we'll integrate by dividing two 5equal sections
            for (double delta = 0.0; delta < 1.0001; delta += 0.05) {
                double delta_t = model->wp.latestMoment - model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment);
                double time = model->tdParams.parameters[paramId].GetX(indexBeforeLastMoment) + delta * delta_t;
                double avg_value = (model->tdParams.parameters[paramId].GetY(indexBeforeLastMoment) * 0.100 +
                                    model->tdParams.parameters[paramId].InterpolateY(time, false) * 0.100) / 2.0;
                newID.emplace_back(time, newID.back().second + 0.05 * delta_t * avg_value);
            }
        }

        return newID;
    }
    
};