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
        auto pos = desorptionParameterIDs.find(paramId);
        if(pos == desorptionParameterIDs.end())
            return -1;
        //--lowerBound; //even temperatureList.end() can be a bound

        if (paramId == *pos) {
            return std::distance(desorptionParameterIDs.begin(), pos);
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
    std::vector<std::pair<double, double>> ID; //time-cumulative desorption pairs. Can have more points than the corresponding outgassing (some sections are divided to subsections)

    //We need to slightly modify the original outgassing:
    //Beginning: Add a point at t=0, with the first outgassing value (assuming constant outgassing from 0 to t1 - const. extrap.)
    //End: if latestMoment is after the last user-defined moment, copy last user value to latestMoment (constant extrapolation)
    //     if latestMoment is before, create a new point at latestMoment with interpolated outgassing value and throw away the rest (don't generate particles after latestMoment)
    std::vector<std::pair<double, double>> myOutgassing; //modified parameter

    //First, let's check at which index is the latest moment
    size_t indexAfterLatestMoment;
    auto &par = model->tdParams.parameters[paramId]; //we'll reference it a lot
    for (indexAfterLatestMoment = 0; indexAfterLatestMoment < par.GetSize() &&
                                     (par.GetX(indexAfterLatestMoment) <
                                      model->wp.latestMoment); indexAfterLatestMoment++); //loop exits after first index after latestMoment
    bool lastUserMomentBeforeLatestMoment;
    if (indexAfterLatestMoment >= par.GetSize()) {
        indexAfterLatestMoment = par.GetSize() - 1; //not found, set as last moment
        lastUserMomentBeforeLatestMoment = true;
    } else {
        lastUserMomentBeforeLatestMoment = false;
    }

//Construct integral from 0 to the simulation's latest moment
    //First point: t=0, Q(0)=Q(t0)
    ID.emplace_back(0.0, 0.0); //At t=0 no particles have desorbed yet
    if (par.GetX(0) > 0.0) { //If the user starts later than t=0, copy first user value to t=0(const extrapolation)
        myOutgassing.emplace_back(0.0, par.GetY(0));
    } else { //user has set an outgassing at t=0, so just copy it to myOutgassing
        myOutgassing.push_back(par.GetValues()[0]);
    }
    //Consecutive points: user-defined points that are before latestMoment
    //We throw away user-defined moments after latestMoment:
    //Example: we sample the system until t=10s but outgassing is defined until t=1000s -> ignore values after 10s
    {
        const auto &valuesCopy = par.GetValues();
        if (lastUserMomentBeforeLatestMoment) {
            myOutgassing.insert(myOutgassing.end(), valuesCopy.begin(), valuesCopy.end()); //copy all values...
            myOutgassing.emplace_back(model->wp.latestMoment,
                                      myOutgassing.back().second); //...and create last point equal to last outgassing (const extrapolation)
        } else {
            if (indexAfterLatestMoment > 0) {
                myOutgassing.insert(myOutgassing.end(), valuesCopy.begin(),
                                    valuesCopy.begin() + indexAfterLatestMoment -
                                    1); //copy values that are before latestMoment
            }
            if (!IsEqual(myOutgassing.back().first, model->wp.latestMoment)) { //if interpolation is needed
                myOutgassing.emplace_back(model->wp.latestMoment,
                                          InterpolateY(model->wp.latestMoment, valuesCopy, par.logXinterp,
                                                       par.logYinterp)); //interpolate outgassing value to t=latestMoment
            }
        }

    } //valuesCopy goes out of scope

    //Intermediate moments, from first to t=latestMoment
    for (size_t i = 1; i < myOutgassing.size(); i++) { //myOutgassing[0] is always at t=0, skipping
        if (IsEqual(myOutgassing[i].second, myOutgassing[i - 1].second)) {
            //easy case of two equal y0=y1 values, simple integration by multiplying, reducing number of points
            ID.emplace_back(myOutgassing[i].first,
                            ID.back().second +
                            (myOutgassing[i].first - myOutgassing[i - 1].first) *
                            myOutgassing[i].second * MBARLS_TO_PAM3S); //integral = y1*(x1-x0)
        } else { //we need to split the user-defined section to 10 subsections, and integrate each
            //(terminology: section is between two consecutive user-defined time-value pairs)
            const int nbSteps = 10; //change here
            double sectionStartTime = myOutgassing[i - 1].first;
            double sectionEndTime = myOutgassing[i].first;
            double sectionTimeInterval = sectionEndTime - sectionStartTime;
            double sectionDelta = 1.0 / nbSteps;
            double subsectionTimeInterval, subsectionLogTimeInterval;
            double logSectionStartTime, logSectionEndTime, logSectionTimeInterval;

            subsectionTimeInterval = subsectionLogTimeInterval =
            logSectionStartTime = logSectionEndTime = logSectionTimeInterval = 0.0;

            if (par.logXinterp) { //time is logarithmic
                if (sectionStartTime == 0) logSectionStartTime = -99;
                else logSectionStartTime = log10(sectionStartTime);
                logSectionEndTime = log10(sectionEndTime);
                logSectionTimeInterval = logSectionEndTime - logSectionStartTime;
                subsectionLogTimeInterval = logSectionTimeInterval * sectionDelta;
            } else {
                subsectionTimeInterval = sectionTimeInterval * sectionDelta;
            }
            double previousSubsectionValue = myOutgassing[i -
                                                          1].second; //Points to start of section, will be updated to point to start of subsections
            double previousSubsectionTime = sectionStartTime;

            double sectionFraction = sectionDelta;
            while (sectionFraction < 1.0001) {
                //for (double sectionFraction = sectionDelta; sectionFraction < 1.0001; sectionFraction += sectionDelta) { //sectionFraction: where we are within the section [0..1]
                double subSectionEndTime;
                if (!par.logXinterp) {
                    //linX: Distribute sampled points evenly
                    subSectionEndTime = sectionStartTime + sectionFraction * sectionTimeInterval;
                } else {
                    //logX: Distribute sampled points logarithmically
                    subSectionEndTime = Pow10(logSectionStartTime + sectionFraction * logSectionTimeInterval);
                }
                double subSectionEndValue = InterpolateY(subSectionEndTime, myOutgassing, par.logXinterp,
                                                         par.logYinterp);
                double subsectionDesorbedGas; //desorbed gas in this subsection, in mbar*l/s, will be converted to Pa*m3/s when adding to ID
                if (!par.logXinterp && !par.logYinterp) { //lin-lin interpolation
                    //Area under a straight section from (x0,y0) to (x1,y1) on a lin-lin plot: I = (x1-x0) * (y0+y1)/2
                    subsectionDesorbedGas =
                            subsectionTimeInterval * 0.5 * (previousSubsectionValue + subSectionEndValue);
                } else if (par.logXinterp &&
                           !par.logYinterp) { //log-lin: time (X) is logarithmic, outgassing (Y) is linear
                    //From Mathematica/WolframAlpha: integral of a straight section from (x0,y0) to (x1,y1) on a log10-lin plot: I = (y1-y0)(x0-x1)/ln(x1/x0) - x0y0 + x1y1
                    double x0 = previousSubsectionTime;
                    double x1 = subSectionEndTime;
                    double y0 = previousSubsectionValue;
                    double y1 = subSectionEndValue;
                    subsectionDesorbedGas = (y1 - y0) * (x0 - x1) / log(x1 / x0) - x0 * y0 + x1 * y1;
                } else if (!par.logXinterp &&
                           par.logYinterp) { //lin-log: time (X) is linear, outgassing (Y) is logarithmic
                    //Area under a straight section from (x0,y0) to (x1,y1) on a lin-log plot: I = 1/m * (y1-y0) where m=log10(y1/y0)/(x1-x0)
                    double logSubSectionEndValue = (subSectionEndValue > 0.0) ? log10(subSectionEndValue) : -99;
                    double logPreviousSubSectionValue = (previousSubsectionTime > 0.0) ? log10(
                            previousSubsectionValue)
                                                                                       : -99;
                    //I = (x2-x1)*(y2-y1)*log10(exp(1))/(log10(y2)-log10(y1)) from https://fr.mathworks.com/matlabcentral/answers/404930-finding-the-area-under-a-semilogy-plot-with-straight-line-segments
                    subsectionDesorbedGas = subsectionTimeInterval * (subSectionEndValue - previousSubsectionValue)
                                            * log10(exp(1)) / (logSubSectionEndValue - logPreviousSubSectionValue);
                } else { //log-log
                    //Area under a straight section from (x0,y0) to (x1,y1) on a log-log plot:
                    //I = (y0/(x0^m))/(m+1) * (x1^(m+1)-x0^(m+1)) if m!=-1
                    //I = x0*y0*log10(x1/x0) if m==-1
                    //where m= log10(y1/y0) / log10(x1/x0)
                    double logSubSectionEndValue = (subSectionEndValue > 0.0) ? log10(subSectionEndValue) : -99;
                    double logPreviousSubSectionValue = (previousSubsectionTime > 0.0) ? log10(
                            previousSubsectionValue)
                                                                                       : -99;
                    double m = (logSubSectionEndValue - logPreviousSubSectionValue) /
                               subsectionLogTimeInterval; //slope
                    if (m != -1.0) {
                        subsectionDesorbedGas =
                                previousSubsectionValue / pow(previousSubsectionTime, m) / (m + 1.0) *
                                (pow(subSectionEndTime, m + 1.0) -
                                 pow(previousSubsectionTime, m + 1.0));
                    } else { //m==-1
                        subsectionDesorbedGas =
                                previousSubsectionTime * previousSubsectionValue * subsectionLogTimeInterval;
                    }
                }

                ID.emplace_back(subSectionEndTime, ID.back().second + subsectionDesorbedGas * MBARLS_TO_PAM3S);

                //Cache end values for next iteration
                previousSubsectionValue = subSectionEndValue;
                previousSubsectionTime = subSectionEndTime;

                sectionFraction += sectionDelta;
            } //end subsection loop
        }
    } //end section loop

    return ID;
}
    
};