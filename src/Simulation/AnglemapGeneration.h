//
// Created by pascal on 2/5/21.
//

#ifndef MOLFLOW_PROJ_ANGLEMAPGENERATION_H
#define MOLFLOW_PROJ_ANGLEMAPGENERATION_H


#include <MolflowTypes.h>
#include <GeometrySimu.h>

class AnglemapGeneration {
public:
    static double GetTheta(const double &thetaIndex, const AnglemapParams &anglemapParams);

    static double GetPhi(const double &phiIndex, const AnglemapParams &anglemapParams);

    static double GetPhipdfValue(const double &thetaIndex, const int &phiLowerIndex,
                                 const AnglemapParams &anglemapParams, const std::vector<size_t> &angleMapPDF);

    static double GetPhiCDFValue(const double &thetaIndex, const int &phiLowerIndex,
                                 const AnglemapParams &anglemapParams, const Anglemap &anglemap);

    static double GetPhiCDFSum(const double &thetaIndex, const AnglemapParams &anglemapParams,
                               const Anglemap &anglemap);

    static std::tuple<double, int, double>
    GenerateThetaFromAngleMap(const AnglemapParams &anglemapParams, const Anglemap &anglemap,
                              const double lookupValue);

    static double GeneratePhiFromAngleMap(const int &thetaLowerIndex, const double &thetaOvershoot,
                                          const AnglemapParams &anglemapParams, Anglemap &anglemap,
                                          const std::vector<size_t> &angleMapPDF,
                                          double lookupValue);
};


#endif //MOLFLOW_PROJ_ANGLEMAPGENERATION_H
