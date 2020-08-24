/*
Program:     MolFlow+ / Synrad+
Description: Monte Carlo simulator for ultra-high vacuum and synchrotron radiation
Authors:     Jean-Luc PONS / Roberto KERSEVAN / Marton ADY
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
#include "Simulation.h"
#include "IntersectAABB_shared.h"
#include "Random.h"
#include "GLApp/MathTools.h"

#include "Parameter.h"
#include <omp.h>

//extern Simulation *sHandle; //delcared in molflowSub.cpp

class AnglemapGeneration {
public:
    static double GetTheta(const double& thetaIndex, const AnglemapParams& anglemapParams);
    static double GetPhi(const double& phiIndex, const AnglemapParams& anglemapParams);
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

/*bool Simulation::UpdateMCHits(Dataport *dpHit, int prIdx, size_t nbMoments, DWORD timeout) {

    BYTE *buffer;
    GlobalHitBuffer *gHits;
    TEXTURE_MIN_MAX texture_limits_old[3];
    int i, j, s, x, y;

    //SetState(PROCESS_STARTING, "Waiting for 'hits' dataport access...", false, true);
    bool lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
    if (!lastHitUpdateOK) return false; //Timeout, will try again later
    //SetState(PROCESS_STARTING, "Updating MC hits...", false, true);
#if defined(_DEBUG)
    double t0, t1;
    t0 = GetTick();
#endif

    buffer = (BYTE *) dpHit->buff;
    gHits = (GlobalHitBuffer *) buffer;

    // Global hits and leaks: adding local hits to shared memory
    int loopInd = 0;
    for(auto& tmpSimuState : tmpGlobalResults) {
        auto& tmpGlobalResult = tmpSimuState.globalHits;
        //std::cout << "["<<loopInd<<"] "<< gHits->globalHits.hit.nbMCHit << " += "<< tmpGlobalResult.globalHits.hit.nbMCHit<<std::endl;
        gHits->globalHits += tmpGlobalResult.globalHits;

        *//*gHits->globalHits.hit.nbMCHit += tmpGlobalResult.globalHits.hit.nbMCHit;
        gHits->globalHits.hit.nbHitEquiv += tmpGlobalResult.globalHits.hit.nbHitEquiv;
        gHits->globalHits.hit.nbAbsEquiv += tmpGlobalResult.globalHits.hit.nbAbsEquiv;
        gHits->globalHits.hit.nbDesorbed += tmpGlobalResult.globalHits.hit.nbDesorbed;*//*
        gHits->distTraveled_total += tmpGlobalResult.distTraveled_total;
        gHits->distTraveledTotal_fullHitsOnly += tmpGlobalResult.distTraveledTotal_fullHitsOnly;

        //Memorize current limits, then do a min/max search
        for (i = 0; i < 3; i++) {
            texture_limits_old[i] = gHits->texture_limits[i];
            gHits->texture_limits[i].min.all = gHits->texture_limits[i].min.moments_only = HITMAX;
            gHits->texture_limits[i].max.all = gHits->texture_limits[i].max.moments_only = 0;
        }

        //model.wp.sMode = MC_MODE;
        //for(i=0;i<BOUNCEMAX;i++) gHits->wallHits[i] += wallHits[i];

        // Leak
        for (size_t leakIndex = 0; leakIndex < tmpGlobalResult.leakCacheSize; leakIndex++)
            gHits->leakCache[(leakIndex + gHits->lastLeakIndex) %
                             LEAKCACHESIZE] = tmpGlobalResult.leakCache[leakIndex];
        gHits->nbLeakTotal += tmpGlobalResult.nbLeakTotal;
        gHits->lastLeakIndex = (gHits->lastLeakIndex + tmpGlobalResult.leakCacheSize) % LEAKCACHESIZE;
        gHits->leakCacheSize = Min(LEAKCACHESIZE, gHits->leakCacheSize + tmpGlobalResult.leakCacheSize);

        // HHit (Only prIdx 0)
        if (prIdx == 0) {
            for (size_t hitIndex = 0; hitIndex < tmpGlobalResult.hitCacheSize; hitIndex++)
                gHits->hitCache[(hitIndex + gHits->lastHitIndex) %
                                HITCACHESIZE] = tmpGlobalResult.hitCache[hitIndex];

            if (tmpGlobalResult.hitCacheSize > 0) {
                gHits->lastHitIndex = (gHits->lastHitIndex + tmpGlobalResult.hitCacheSize) % HITCACHESIZE;
                gHits->hitCache[gHits->lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
                gHits->hitCacheSize = Min(HITCACHESIZE, gHits->hitCacheSize + tmpGlobalResult.hitCacheSize);
            }
        }

        //Global histograms

        for (int m = 0; m < (1 + nbMoments); m++) {
            BYTE *histCurrentMoment =
                    buffer + sizeof(GlobalHitBuffer) + m * model.wp.globalHistogramParams.GetDataSize();
            auto& tmpHistogram = tmpSimuState.globalHistograms[m];
            if (model.wp.globalHistogramParams.recordBounce) {
                double *nbHitsHistogram = (double *) histCurrentMoment;
                for (size_t i = 0; i < model.wp.globalHistogramParams.GetBounceHistogramSize(); i++) {
                    nbHitsHistogram[i] += tmpHistogram.nbHitsHistogram[i];
                }
            }
            if (model.wp.globalHistogramParams.recordDistance) {
                double *distanceHistogram = (double *) (histCurrentMoment +
                                                        model.wp.globalHistogramParams.GetBouncesDataSize());
                for (size_t i = 0; i < (model.wp.globalHistogramParams.GetDistanceHistogramSize()); i++) {
                    distanceHistogram[i] += tmpHistogram.distanceHistogram[i];
                }
            }
            if (model.wp.globalHistogramParams.recordTime) {
                double *timeHistogram = (double *) (histCurrentMoment +
                                                    model.wp.globalHistogramParams.GetBouncesDataSize() +
                                                    model.wp.globalHistogramParams.GetDistanceDataSize());
                for (size_t i = 0; i < (model.wp.globalHistogramParams.GetTimeHistogramSize()); i++) {
                    timeHistogram[i] += tmpHistogram.timeHistogram[i];
                }
            }
        }


        size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);
        // Facets
        for (s = 0; s < model.sh.nbSuper; s++) {
            for (SubprocessFacet &f : model.structures[s].facets) {
                if (f.isHit) {
                    BYTE *bufferStart = buffer + f.sh.hitOffset;
                    for (int m = 0; m < (1 + nbMoments); m++) {
                        FacetHitBuffer *facetHitBuffer = (FacetHitBuffer *) (bufferStart + m * sizeof(FacetHitBuffer));
                        (*facetHitBuffer) += f.tmpCounter[m];
                    }

                    if (f.sh.isProfile) {
                        bufferStart = buffer + f.sh.hitOffset + facetHitsSize;
                        for (int m = 0; m < (1 + nbMoments); m++) {
                            ProfileSlice *shProfile = (ProfileSlice *) (bufferStart + m * f.profileSize);
                            for (j = 0; j < PROFILE_SIZE; j++) {
                                shProfile[j] += f.profile[m][j];
                            }
                        }
                    }

                    if (f.sh.isTextured) {
                        bufferStart = buffer + (f.sh.hitOffset + facetHitsSize + f.profileSize * (1 + nbMoments));
                        for (int m = 0; m < (1 + nbMoments); m++) {
                            TextureCell *shTexture = (TextureCell *) (bufferStart + (m * f.textureSize));
                            //double dCoef = gHits->globalHits.hit.nbDesorbed * 1E4 * model.wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
                            double timeCorrection =
                                    m == 0 ? model.wp.finalOutgassingRate : (model.wp.totalDesorbedMolecules) /
                                                                            model.tdParams.moments[m - 1].second;
                            //model.wp.timeWindowSize;
                            //Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)

                            for (y = 0; y < f.sh.texHeight; y++) {
                                for (x = 0; x < f.sh.texWidth; x++) {
                                    size_t add = x + y * f.sh.texWidth;

                                    //Add temporary hit counts
                                    shTexture[add] += f.texture[m][add];

                                    double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)

                                    val[0] = shTexture[add].sum_v_ort_per_area *
                                             timeCorrection; //pressure without dCoef_pressure
                                    val[1] = shTexture[add].countEquiv * f.textureCellIncrements[add] *
                                             timeCorrection; //imp.rate without dCoef
                                    val[2] = f.textureCellIncrements[add] * shTexture[add].sum_1_per_ort_velocity *
                                             timeCorrection; //particle density without dCoef

                                    //Global autoscale
                                    for (int v = 0; v < 3; v++) {
                                        if (val[v] > gHits->texture_limits[v].max.all && f.largeEnough[add])
                                            gHits->texture_limits[v].max.all = val[v];

                                        if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.all &&
                                            f.largeEnough[add])
                                            gHits->texture_limits[v].min.all = val[v];

                                        //Autoscale ignoring constant flow (moments only)
                                        if (m != 0) {
                                            if (val[v] > gHits->texture_limits[v].max.moments_only &&
                                                f.largeEnough[add])
                                                gHits->texture_limits[v].max.moments_only = val[v];

                                            if (val[v] > 0.0 && val[v] < gHits->texture_limits[v].min.moments_only &&
                                                f.largeEnough[add])
                                                gHits->texture_limits[v].min.moments_only = val[v];
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (f.sh.countDirection) {
                        bufferStart = buffer + (f.sh.hitOffset + facetHitsSize +
                                                f.profileSize * (1 + nbMoments) +
                                                f.textureSize * (1 + nbMoments));
                        for (int m = 0; m < (1 + nbMoments); m++) {
                            DirectionCell *shDir = (DirectionCell *) (bufferStart + f.directionSize * m);
                            for (y = 0; y < f.sh.texHeight; y++) {
                                for (x = 0; x < f.sh.texWidth; x++) {
                                    size_t add = x + y * f.sh.texWidth;
                                    shDir[add].dir.x += f.direction[m][add].dir.x;
                                    shDir[add].dir.y += f.direction[m][add].dir.y;
                                    shDir[add].dir.z += f.direction[m][add].dir.z;
                                    //shDir[add].sumSpeed += f.direction[m][add].sumSpeed;
                                    shDir[add].count += f.direction[m][add].count;
                                }
                            }
                        }
                    }

                    if (f.sh.anglemapParams.record) {
                        size_t *shAngleMap = (size_t *) (buffer + f.sh.hitOffset + facetHitsSize +
                                                         f.profileSize * (1 + nbMoments) +
                                                         f.textureSize * (1 + nbMoments) +
                                                         f.directionSize * (1 + nbMoments));
                        for (y = 0; y < (f.sh.anglemapParams.thetaLowerRes + f.sh.anglemapParams.thetaHigherRes); y++) {
                            for (x = 0; x < f.sh.anglemapParams.phiWidth; x++) {
                                size_t add = x + y * f.sh.anglemapParams.phiWidth;
                                shAngleMap[add] += f.angleMap.pdf[add];
                            }
                        }
                    }

                    //Facet histograms

                    bufferStart = buffer + f.sh.hitOffset + facetHitsSize + f.profileSize * (1 + nbMoments) +
                                  f.textureSize * (1 + nbMoments) + f.directionSize * (1 + nbMoments);
                    for (int m = 0; m < (1 + nbMoments); m++) {
                        size_t angleMapRecordedDataSize = sizeof(size_t) * (f.sh.anglemapParams.phiWidth *
                                                                            (f.sh.anglemapParams.thetaLowerRes +
                                                                             f.sh.anglemapParams.thetaHigherRes));
                        BYTE *histCurrentMoment =
                                bufferStart + angleMapRecordedDataSize + m * f.sh.facetHistogramParams.GetDataSize();
                        if (f.sh.facetHistogramParams.recordBounce) {
                            double *nbHitsHistogram = (double *) histCurrentMoment;
                            for (size_t i = 0; i < f.sh.facetHistogramParams.GetBounceHistogramSize(); i++) {
                                nbHitsHistogram[i] += f.tmpHistograms[m].nbHitsHistogram[i];
                            }
                        }
                        if (f.sh.facetHistogramParams.recordDistance) {
                            double *distanceHistogram = (double *) (histCurrentMoment +
                                                                    f.sh.facetHistogramParams.GetBouncesDataSize());
                            for (size_t i = 0; i < (f.sh.facetHistogramParams.GetDistanceHistogramSize()); i++) {
                                distanceHistogram[i] += f.tmpHistograms[m].distanceHistogram[i];
                            }
                        }
                        if (f.sh.facetHistogramParams.recordTime) {
                            double *timeHistogram = (double *) (histCurrentMoment +
                                                                f.sh.facetHistogramParams.GetBouncesDataSize() +
                                                                f.sh.facetHistogramParams.GetDistanceDataSize());
                            for (size_t i = 0; i < (f.sh.facetHistogramParams.GetTimeHistogramSize()); i++) {
                                timeHistogram[i] += f.tmpHistograms[m].timeHistogram[i];
                            }
                        }
                    }

                } // End if(isHit)
            } // End nbFacet
        } // End nbSuper
    }
    //if there were no textures:
    for (int v = 0; v < 3; v++) {
        if (gHits->texture_limits[v].min.all == HITMAX)
            gHits->texture_limits[v].min.all = texture_limits_old[v].min.all;
        if (gHits->texture_limits[v].min.moments_only == HITMAX)
            gHits->texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
        if (gHits->texture_limits[v].max.all == 0.0)
            gHits->texture_limits[v].max.all = texture_limits_old[v].max.all;
        if (gHits->texture_limits[v].max.moments_only == 0.0)
            gHits->texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
    }


    ReleaseDataport(dpHit);

    //extern char *GetSimuStatus();
    //SetState(PROCESS_STARTING, GetSimuStatus(), false, true);

#if defined(_DEBUG)
    t1 = GetTick();
    printf("Update hits: %f us\n", (t1 - t0) * 1000000.0);
#endif

    return true;
}*/

bool Simulation::UpdateMCHits(GlobalSimuState& globState, int prIdx, size_t nbMoments, DWORD timeout) {

    TEXTURE_MIN_MAX texture_limits_old[3];
    int i, j, s, x, y;
#if defined(_DEBUG)
    double t0, t1;
    t0 = GetTick();
#endif
    //SetState(PROCESS_STARTING, "Waiting for 'hits' dataport access...", false, true);
    // TODO: Mutex
    //bool lastHitUpdateOK = LockMutex(worker->results.mutex, timeout);
    //if (!lastHitUpdateOK) return false; //Timeout, will try again later
    //SetState(PROCESS_STARTING, "Updating MC hits...", false, true);

    //buffer = (BYTE *) dpHit->buff;
    //gHits = (GlobalHitBuffer *) buffer;

    //Memorize current limits, then do a min/max search
    for (i = 0; i < 3; i++) {
        texture_limits_old[i] = globState.globalHits.texture_limits[i];
        globState.globalHits.texture_limits[i].min.all = globState.globalHits.texture_limits[i].min.moments_only = HITMAX;
        globState.globalHits.texture_limits[i].max.all = globState.globalHits.texture_limits[i].max.moments_only = 0;
    }

    // Global hits and leaks: adding local hits to shared memory
    /*gHits->globalHits += tmpGlobalResult.globalHits;
    gHits->distTraveled_total += tmpGlobalResult.distTraveled_total;
    gHits->distTraveledTotal_fullHitsOnly += tmpGlobalResult.distTraveledTotal_fullHitsOnly;*/
    for(auto globIter = tmpGlobalResults.begin()+1;globIter!=tmpGlobalResults.end();++globIter){
    //for(auto& tmpResults : tmpGlobalResults) {
        auto& tmpResults = *globIter;
        globState.globalHits.globalHits += tmpResults.globalHits.globalHits;
        globState.globalHits.distTraveled_total += tmpResults.globalHits.distTraveled_total;
        globState.globalHits.distTraveledTotal_fullHitsOnly += tmpResults.globalHits.distTraveledTotal_fullHitsOnly;

        /*gHits->globalHits.hit.nbMCHit += tmpGlobalResult.globalHits.hit.nbMCHit;
        gHits->globalHits.hit.nbHitEquiv += tmpGlobalResult.globalHits.hit.nbHitEquiv;
        gHits->globalHits.hit.nbAbsEquiv += tmpGlobalResult.globalHits.hit.nbAbsEquiv;
        gHits->globalHits.hit.nbDesorbed += tmpGlobalResult.globalHits.hit.nbDesorbed;*/

        //model.wp.sMode = MC_MODE;
        //for(i=0;i<BOUNCEMAX;i++) globState.globalHits.wallHits[i] += wallHits[i];

        // Leak
        for (size_t leakIndex = 0; leakIndex < tmpResults.globalHits.leakCacheSize; leakIndex++)
            globState.globalHits.leakCache[(leakIndex + globState.globalHits.lastLeakIndex) %
                                           LEAKCACHESIZE] = tmpResults.globalHits.leakCache[leakIndex];
        globState.globalHits.nbLeakTotal += tmpResults.globalHits.nbLeakTotal;
        globState.globalHits.lastLeakIndex =
                (globState.globalHits.lastLeakIndex + tmpResults.globalHits.leakCacheSize) % LEAKCACHESIZE;
        globState.globalHits.leakCacheSize = Min(LEAKCACHESIZE, globState.globalHits.leakCacheSize +
                                                                tmpResults.globalHits.leakCacheSize);

        // HHit (Only prIdx 0)
        if (prIdx == 0) {
            for (size_t hitIndex = 0; hitIndex < tmpResults.globalHits.hitCacheSize; hitIndex++)
                globState.globalHits.hitCache[(hitIndex + globState.globalHits.lastHitIndex) %
                                              HITCACHESIZE] = globState.globalHits.hitCache[hitIndex];

            if (tmpResults.globalHits.hitCacheSize > 0) {
                globState.globalHits.lastHitIndex =
                        (globState.globalHits.lastHitIndex + tmpResults.globalHits.hitCacheSize) % HITCACHESIZE;
                globState.globalHits.hitCache[globState.globalHits.lastHitIndex].type = HIT_LAST; //Penup (border between blocks of consecutive hits in the hit cache)
                globState.globalHits.hitCacheSize = Min(HITCACHESIZE, globState.globalHits.hitCacheSize +
                                                                      tmpResults.globalHits.hitCacheSize);
            }
        }

        //Global histograms
        globState.globalHistograms += tmpResults.globalHistograms;

        // Facets
        globState.facetStates += tmpResults.facetStates;
    }

    for (s = 0; s < model.sh.nbSuper; s++) {
        for (SubprocessFacet &f : model.structures[s].facets) {
            if (f.isHit) {
                if (f.sh.isTextured) {
                    for (int m = 0; m < (1 + nbMoments); m++) {
                        //double dCoef = globState.globalHits.globalHits.hit.nbDesorbed * 1E4 * model.wp.gasMass / 1000 / 6E23 * MAGIC_CORRECTION_FACTOR;  //1E4 is conversion from m2 to cm2
                        const double timeCorrection =
                                m == 0 ? model.wp.finalOutgassingRate : (model.wp.totalDesorbedMolecules) /
                                                                        model.tdParams.moments[m - 1].second;
                        //model.wp.timeWindowSize;
                        //Timecorrection is required to compare constant flow texture values with moment values (for autoscaling)
                        const auto &texture = globState.facetStates[f.globalId].momentResults[m].texture;
                        const size_t textureSize = texture.size();
                        for (size_t t = 0; t < textureSize; t++) {
                            double val[3];  //pre-calculated autoscaling values (Pressure, imp.rate, density)
                            val[0] = texture[t].sum_v_ort_per_area *
                                     timeCorrection; //pressure without dCoef_pressure
                            val[1] = texture[t].countEquiv * f.textureCellIncrements[t] *
                                     timeCorrection; //imp.rate without dCoef
                            val[2] = f.textureCellIncrements[t] * texture[t].sum_1_per_ort_velocity *
                                     timeCorrection; //particle density without dCoef

                            //Global autoscale
                            for (int v = 0; v < 3; v++) {
                                if (val[v] > globState.globalHits.texture_limits[v].max.all && f.largeEnough[t])
                                    globState.globalHits.texture_limits[v].max.all = val[v];

                                if (val[v] > 0.0 && val[v] < globState.globalHits.texture_limits[v].min.all &&
                                    f.largeEnough[t])
                                    globState.globalHits.texture_limits[v].min.all = val[v];

                                //Autoscale ignoring constant flow (moments only)
                                if (m != 0) {
                                    if (val[v] > globState.globalHits.texture_limits[v].max.moments_only &&
                                        f.largeEnough[t])
                                        globState.globalHits.texture_limits[v].max.moments_only = val[v];

                                    if (val[v] > 0.0 &&
                                        val[v] < globState.globalHits.texture_limits[v].min.moments_only &&
                                        f.largeEnough[t])
                                        globState.globalHits.texture_limits[v].min.moments_only = val[v];
                                }
                            }
                        }
                    }
                }

                //Facet histograms -> with FacetState+=

            } // End if(isHit)
        } // End nbFacet
    } // End nbSuper

    //if there were no textures:
    for (int v = 0; v < 3; v++) {
        if (globState.globalHits.texture_limits[v].min.all == HITMAX)
            globState.globalHits.texture_limits[v].min.all = texture_limits_old[v].min.all;
        if (globState.globalHits.texture_limits[v].min.moments_only == HITMAX)
            globState.globalHits.texture_limits[v].min.moments_only = texture_limits_old[v].min.moments_only;
        if (globState.globalHits.texture_limits[v].max.all == 0.0)
            globState.globalHits.texture_limits[v].max.all = texture_limits_old[v].max.all;
        if (globState.globalHits.texture_limits[v].max.moments_only == 0.0)
            globState.globalHits.texture_limits[v].max.moments_only = texture_limits_old[v].max.moments_only;
    }

    //ReleaseDataport(dpHit);
    //TODO: //ReleaseMutex(worker->results.mutex);

    //extern char *GetSimuStatus();
    //SetState(PROCESS_STARTING, GetSimuStatus(), false, true);

#if defined(_DEBUG)
    t1 = GetTick();
    printf("Update hits (glob): %f us\n", (t1 - t0) * 1000000.0);
#endif

    return true;
}

/*
bool Simulation::UploadHits(Dataport *dpHit, Dataport* dpLog, int prIdx, DWORD timeout) {
    // Prepare
#if defined(_DEBUG)
    double t0, t1;
    t0 = GetTick();
#endif
    bool lastHitUpdateOK = AccessDataportTimed(dpHit, timeout);
    if (!lastHitUpdateOK) return false; //Timeout, will try again later

    // Prepare data
    std::ostringstream result;
    {
        cereal::BinaryOutputArchive outputArchive(result);

        // Reduce to first threads storage and transfer
        GlobalSimuState& tmpResults = tmpGlobalResults[0];
        for(size_t ind=1;ind < tmpGlobalResults.size();++ind) {
            tmpResults += tmpGlobalResults[ind];
        }

        outputArchive(
                cereal::make_nvp("GlobalHits", tmpResults.globalHits),
                cereal::make_nvp("GlobalHistograms", tmpResults.globalHistograms),
                cereal::make_nvp("FacetStates", tmpResults.facetStates)
        );
    }

    // Upload
    if(dpHit){
        BYTE* buffer = (BYTE*)dpHit->buff;
        auto resString = result.str();
        BYTE* source = (BYTE*)(resString.data());
        auto resultSize =  resString.size();
        std::copy(source,source + resultSize, buffer);
    }

    // Finalize
    ReleaseDataport(dpHit);
#if defined(_DEBUG)
    t1 = GetTick();
    printf("Upload hits: %f us\n", (t1 - t0) * 1000000.0);
#endif

    return true;
}
*/

/*void Simulation::UpdateLog(Dataport *dpLog, DWORD timeout) {
    if (tmpParticleLog.size()) {
#if defined(_DEBUG)
        double t0, t1;
        t0 = GetTick();
#endif
        //SetState(PROCESS_STARTING, "Waiting for 'dpLog' dataport access...", false, true);
        lastLogUpdateOK = AccessDataportTimed(dpLog, timeout);
        //SetState(PROCESS_STARTING, "Updating Log...", false, true);
        if (!lastLogUpdateOK) return;

        size_t *logBuff = (size_t *) dpLog->buff;
        size_t recordedLogSize = *logBuff;
        ParticleLoggerItem *logBuff2 = (ParticleLoggerItem *) (logBuff + 1);

        size_t writeNb;
        if (recordedLogSize > model.otfParams.logLimit) writeNb = 0;
        else writeNb = Min(tmpParticleLog.size(), model.otfParams.logLimit - recordedLogSize);
        memcpy(&logBuff2[recordedLogSize], &tmpParticleLog[0],
               writeNb * sizeof(ParticleLoggerItem)); //Knowing that vector memories are contigious
        (*logBuff) += writeNb;
        ReleaseDataport(dpLog);
        tmpParticleLog.clear();
        //extern char *GetSimuStatus();
        //SetState(PROCESS_STARTING, GetSimuStatus(), false, true);

#if defined(_DEBUG)
        t1 = GetTick();
        printf("Update log: %f us\n", (t1 - t0) * 1000000.0);
#endif
    }
}*/

// Compute particle teleport

void Simulation::PerformTeleport(SubprocessFacet *iFacet, CurrentParticleStatus &currentParticle) {
    
    //Search destination
    SubprocessFacet *destination;
    bool found = false;
    bool revert = false;
    int destIndex;
    if (iFacet->sh.teleportDest == -1) {
        destIndex = currentParticle.teleportedFrom;
        if (destIndex == -1) {
            /*char err[128];
            sprintf(err, "Facet %d tried to teleport to the facet where the particle came from, but there is no such facet.", iFacet->globalId + 1);
            SetErrorSub(err);*/
            RecordHit(HIT_REF, currentParticle);
            currentParticle.lastHitFacet = iFacet;
            return; //LEAK
        }
    } else destIndex = iFacet->sh.teleportDest - 1;

    //Look in which superstructure is the destination facet:
    for (size_t i = 0; i < model.sh.nbSuper && (!found); i++) {
        for (size_t j = 0; j < model.structures[i].facets.size() && (!found); j++) {
            if (destIndex == model.structures[i].facets[j].globalId) {
                destination = &(model.structures[i].facets[j]);
                if (destination->sh.superIdx != -1) {
                    currentParticle.structureId = destination->sh.superIdx; //change current superstructure, unless the target is a universal facet
                }
                currentParticle.teleportedFrom = (int) iFacet->globalId; //memorize where the particle came from
                found = true;
            }
        }
    }
    if (!found) {
        /*char err[128];
        sprintf(err, "Teleport destination of facet %d not found (facet %d does not exist)", iFacet->globalId + 1, iFacet->sh.teleportDest);
        SetErrorSub(err);*/
        RecordHit(HIT_REF, currentParticle);
        currentParticle.lastHitFacet = iFacet;
        return; //LEAK
    }
    // Count this hit as a transparent pass
    RecordHit(HIT_TELEPORTSOURCE, currentParticle);
    if (/*iFacet->texture && */iFacet->sh.countTrans)
        RecordHitOnTexture(iFacet, currentParticle.flightTime, true, 2.0, 2.0, currentParticle);
    if (/*iFacet->direction && */iFacet->sh.countDirection)
        RecordDirectionVector(iFacet, currentParticle.flightTime, currentParticle);
    ProfileFacet(iFacet, currentParticle.flightTime, true, 2.0, 2.0, currentParticle);
    LogHit(iFacet, currentParticle);
    if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet, currentParticle);

    // Relaunch particle from new facet
    auto[inTheta, inPhi] = CartesianToPolar(currentParticle.direction, iFacet->sh.nU, iFacet->sh.nV,
                                            iFacet->sh.N);
    currentParticle.direction = PolarToCartesian(destination, inTheta, inPhi, false);
    // Move particle to teleport destination point
    double u = iFacet->colU;
    double v = iFacet->colV;
    currentParticle.position = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
    RecordHit(HIT_TELEPORTDEST, currentParticle);
    int nbTry = 0;
    if (!IsInFacet(*destination, u, v)) { //source and destination facets not the same shape, would generate leak
        // Choose a new starting point
        RecordHit(HIT_ABS, currentParticle);
        bool found = false;
        while (!found && nbTry < 1000) {
            u = currentParticle.randomGenerator.rnd();
            v = currentParticle.randomGenerator.rnd();
            if (IsInFacet(*destination, u, v)) {
                found = true;
                currentParticle.position = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
                RecordHit(HIT_DES, currentParticle);
            }
        }
        nbTry++;
    }

    currentParticle.lastHitFacet = destination;

    //Count hits on teleport facets
    /*iFacet->sh.tmpCounter.hit.nbAbsEquiv++;
    destination->sh.tmpCounter.hit.nbDesorbed++;*/

    double ortVelocity =
            currentParticle.velocity * std::abs(Dot(currentParticle.direction, iFacet->sh.N));
    //We count a teleport as a local hit, but not as a global one since that would affect the MFP calculation
    /*iFacet->sh.tmpCounter.hit.nbMCHit++;
    iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity;
    iFacet->sh.tmpCounter.hit.sum_v_ort += 2.0*(model.wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
    IncreaseFacetCounter(iFacet, currentParticle.flightTime, 1, 0, 0, 2.0 / ortVelocity,
                         2.0 * (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity, currentParticle);
    iFacet->isHit = true;
    /*destination->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / currentParticle.velocity;
    destination->sh.tmpCounter.hit.sum_v_ort += currentParticle.velocity*abs(DOT3(
    currentParticle.direction.x, currentParticle.direction.y, currentParticle.direction.z,
    destination->sh.N.x, destination->sh.N.y, destination->sh.N.z));*/
}

// Perform nbStep simulation steps (a step is a bounce)

bool Simulation::SimulationMCStep(size_t nbStep) {

    // Check end of simulation
    if (model.otfParams.desorptionLimit > 0) {
        if (totalDesorbed >= model.otfParams.desorptionLimit / model.otfParams.nbProcess) {
            //currentParticle.lastHitFacet = nullptr; // reset full particle status or go on from where we left
            return false;
        }
    }

    // Perform simulation steps
    int returnVal = true;
    int allQuit = 0;
#pragma omp parallel num_threads(nbThreads) default(none) firstprivate(nbStep) shared( returnVal, allQuit)
    {

#if defined(_DEBUG)
        double t0, t1;
        t0 = GetTick();
#endif

        int index = omp_get_thread_num();
        CurrentParticleStatus& currentParticle = currentParticles[index];

        size_t i;
        for (i = 0; i < nbStep && allQuit <= 0; i++) {
            if (!currentParticle.lastHitFacet)
                if (!StartFromSource(currentParticle)) {
                    returnVal = false; // desorp limit reached
                    break;
                }
            //return (currentParticle.lastHitFacet != nullptr);

            //Prepare output values
            auto[found, collidedFacet, d] = Intersect(this, currentParticle.position,
                                                      currentParticle.direction, currentParticle);

            if (found) {

                // Move particle to intersection point
                currentParticle.position =
                        currentParticle.position + d * currentParticle.direction;
                //currentParticle.distanceTraveled += d;

                double lastFlightTime = currentParticle.flightTime; //memorize for partial hits
                currentParticle.flightTime +=
                        d / 100.0 / currentParticle.velocity; //conversion from cm to m

                if ((!model.wp.calcConstantFlow && (currentParticle.flightTime > model.wp.latestMoment))
                    || (model.wp.enableDecay &&
                        (currentParticle.expectedDecayMoment < currentParticle.flightTime))) {
                    //hit time over the measured period - we create a new particle
                    //OR particle has decayed
                    double remainderFlightPath = currentParticle.velocity * 100.0 *
                                                 Min(model.wp.latestMoment - lastFlightTime,
                                                     currentParticle.expectedDecayMoment -
                                                     lastFlightTime); //distance until the point in space where the particle decayed
                    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
                    tmpResults.globalHits.distTraveled_total += remainderFlightPath * currentParticle.oriRatio;
                    RecordHit(HIT_LAST, currentParticle);
                    //distTraveledSinceUpdate += currentParticle.distanceTraveled;
                    if (!StartFromSource(currentParticle)) {
                        // desorptionLimit reached
                        returnVal = false;
                        break;
                    }
                } else { //hit within measured time, particle still alive
                    if (collidedFacet->sh.teleportDest != 0) { //Teleport
                        IncreaseDistanceCounters(d * currentParticle.oriRatio, currentParticle);
                        PerformTeleport(collidedFacet, currentParticle);
                    }
                        /*else if ((GetOpacityAt(collidedFacet, currentParticle.flightTime) < 1.0) && (randomGenerator.rnd() > GetOpacityAt(collidedFacet, currentParticle.flightTime))) {
                            //Transparent pass
                            tmpResults.globalHits.distTraveled_total += d;
                            PerformTransparentPass(collidedFacet);
                        }*/
                    else { //Not teleport
                        IncreaseDistanceCounters(d * currentParticle.oriRatio, currentParticle);
                        double stickingProbability = GetStickingAt(collidedFacet, currentParticle.flightTime);
                        if (!model.otfParams.lowFluxMode) { //Regular stick or bounce
                            if (stickingProbability == 1.0 ||
                                ((stickingProbability > 0.0) && (currentParticle.randomGenerator.rnd() < (stickingProbability)))) {
                                //Absorbed
                                RecordAbsorb(collidedFacet, currentParticle);
                                //distTraveledSinceUpdate += currentParticle.distanceTraveled;
                                if (!StartFromSource(currentParticle)) {
                                    // desorptionLimit reached
                                    returnVal = false;
                                    break;
                                }
                            } else {
                                //Reflected
                                PerformBounce(collidedFacet, currentParticle);
                            }
                        } else { //Low flux mode
                            if (stickingProbability > 0.0) {
                                double oriRatioBeforeCollision = currentParticle.oriRatio; //Local copy
                                currentParticle.oriRatio *= (stickingProbability); //Sticking part
                                RecordAbsorb(collidedFacet, currentParticle);
                                currentParticle.oriRatio =
                                        oriRatioBeforeCollision * (1.0 - stickingProbability); //Reflected part
                            } else
                                currentParticle.oriRatio *= (1.0 - stickingProbability);
                            if (currentParticle.oriRatio > model.otfParams.lowFluxCutoff) {
                                PerformBounce(collidedFacet, currentParticle);
                            } else { //eliminate remainder and create new particle
                                if (!StartFromSource(currentParticle)) {
                                    // desorptionLimit reached
                                    returnVal = false;
                                    break;
                                }
                            }
                        }
                    }
                } //end hit within measured time
            } //end intersection found
            else {
                // No intersection found: Leak
                GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
                tmpResults.globalHits.nbLeakTotal++;
                RecordLeakPos(currentParticle);
                if (!StartFromSource(currentParticle)) {
                    // desorptionLimit reached
                    returnVal = false;
                    break;
                }
            }
        }

#pragma omp critical
            ++allQuit;
/*#if defined(_DEBUG)
        t1 = GetTick();
        printf("[%d] done after %zu steps: %f s\n", index, i,(t1 - t0) * 1.0);
#endif*/
    } // omp parallel
    return returnVal;
}

void Simulation::IncreaseDistanceCounters(double distanceIncrement, CurrentParticleStatus &currentParticle) {
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    tmpResults.globalHits.distTraveled_total += distanceIncrement;
    tmpResults.globalHits.distTraveledTotal_fullHitsOnly += distanceIncrement;
    currentParticle.distanceTraveled += distanceIncrement;
}

// Launch a ray from a source facet. The ray 
// direction is chosen according to the desorption type.

bool Simulation::StartFromSource(CurrentParticleStatus &currentParticle) {
    bool found = false;
    bool foundInMap = false;
    bool reverse;
    size_t mapPositionW, mapPositionH;
    SubprocessFacet *src = nullptr;
    double srcRnd;
    double sumA = 0.0;
    int i = 0, j = 0;
    int nbTry = 0;

    // Check end of simulation
    if (model.otfParams.desorptionLimit > 0) {
        if (totalDesorbed >= model.otfParams.desorptionLimit / model.otfParams.nbProcess) {
            //currentParticle.lastHitFacet = nullptr; // reset full particle status or go on from where we left
            return false;
        }
    }

    // Select source
    srcRnd = currentParticle.randomGenerator.rnd() * model.wp.totalDesorbedMolecules;

    while (!found && j < model.sh.nbSuper) { //Go through superstructures
        i = 0;
        while (!found && i < model.structures[j].facets.size()) { //Go through facets in a structure
            SubprocessFacet &f = model.structures[j].facets[i];
            if (f.sh.desorbType != DES_NONE) { //there is some kind of outgassing
                if (f.sh.useOutgassingFile) { //Using SynRad-generated outgassing map
                    if (f.sh.totalOutgassing > 0.0) {
                        found = (srcRnd >= sumA) && (srcRnd < (sumA + model.wp.latestMoment * f.sh.totalOutgassing /
                                                                      (1.38E-23 * f.sh.temperature)));
                        if (found) {
                            //look for exact position in map
                            double rndRemainder = (srcRnd - sumA) / model.wp.latestMoment * (1.38E-23 *
                                                                                                f.sh.temperature); //remainder, should be less than f.sh.totalOutgassing
                            /*double sumB = 0.0;
                            for (w = 0; w < f.sh.outgassingMapWidth && !foundInMap; w++) {
                                for (h = 0; h < f.sh.outgassingMapHeight && !foundInMap; h++) {
                                    double cellOutgassing = f.outgassingMap[h*f.sh.outgassingMapWidth + w];
                                    if (cellOutgassing > 0.0) {
                                        foundInMap = (rndRemainder >= sumB) && (rndRemainder < (sumB + cellOutgassing));
                                        if (foundInMap) mapPositionW = w; mapPositionH = h;
                                        sumB += cellOutgassing;
                                    }
                                }
                            }*/
                            double lookupValue = rndRemainder;
                            int outgLowerIndex = my_lower_bound(lookupValue,
                                                                f.outgassingMap); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
                            outgLowerIndex++;
                            mapPositionH = (size_t)((double) outgLowerIndex / (double) f.sh.outgassingMapWidth);
                            mapPositionW = (size_t) outgLowerIndex - mapPositionH * f.sh.outgassingMapWidth;
                            foundInMap = true;
                            /*if (!foundInMap) {
                                SetErrorSub("Starting point not found in imported desorption map");
                                return false;
                            }*/
                        }
                        sumA += model.wp.latestMoment * f.sh.totalOutgassing / (1.38E-23 * f.sh.temperature);
                    }
                } //end outgassing file block
                else { //constant or time-dependent outgassing
                    double facetOutgassing =
                            (f.sh.outgassing_paramId >= 0)
                            ? model.tdParams.IDs[f.sh.IDid].back().second / (1.38E-23 * f.sh.temperature)
                            : model.wp.latestMoment * f.sh.outgassing / (1.38E-23 * f.sh.temperature);
                    found = (srcRnd >= sumA) && (srcRnd < (sumA + facetOutgassing));
                    sumA += facetOutgassing;
                } //end constant or time-dependent outgassing block
            } //end 'there is some kind of outgassing'
            if (!found) i++;
            if (f.sh.is2sided) reverse = currentParticle.randomGenerator.rnd() > 0.5;
            else reverse = false;
        }
        if (!found) j++;
    }
    if (!found) {
        std::cerr << "No starting point, aborting" << std::endl;
        //SetErrorSub("No starting point, aborting");
        return false;
    }
    src = &(model.structures[j].facets[i]);

    currentParticle.lastHitFacet = src;
    //currentParticle.distanceTraveled = 0.0;  //for mean free path calculations
    //currentParticle.flightTime = desorptionStartTime + (desorptionStopTime - desorptionStartTime)*currentParticle.randomGenerator.rnd();
    currentParticle.flightTime = GenerateDesorptionTime(src, currentParticle.randomGenerator.rnd());
    currentParticle.lastMomentIndex = 0;
    if (model.wp.useMaxwellDistribution) currentParticle.velocity = GenerateRandomVelocity(src->sh.CDFid, currentParticle.randomGenerator.rnd());
    else
        currentParticle.velocity =
                145.469 * std::sqrt(src->sh.temperature / model.wp.gasMass);  //sqrt(8*R/PI/1000)=145.47
    currentParticle.oriRatio = 1.0;
    if (model.wp.enableDecay) { //decaying gas
        currentParticle.expectedDecayMoment =
                currentParticle.flightTime + model.wp.halfLife * 1.44269 * -log(currentParticle.randomGenerator.rnd()); //1.44269=1/ln2
        //Exponential distribution PDF: probability of 't' life = 1/TAU*exp(-t/TAU) where TAU = half_life/ln2
        //Exponential distribution CDF: probability of life shorter than 't" = 1-exp(-t/TAU)
        //Equation: currentParticle.randomGenerator.rnd()=1-exp(-t/TAU)
        //Solution: t=TAU*-log(1-currentParticle.randomGenerator.rnd()) and 1-currentParticle.randomGenerator.rnd()=currentParticle.randomGenerator.rnd() therefore t=half_life/ln2*-log(currentParticle.randomGenerator.rnd())
    } else {
        currentParticle.expectedDecayMoment = 1e100; //never decay
    }
    //temperature = src->sh.temperature; //Thermalize particle
    currentParticle.nbBounces = 0;
    currentParticle.distanceTraveled = 0;

    found = false; //Starting point within facet

    // Choose a starting point
    while (!found && nbTry < 1000) {
        double u, v;

        if (foundInMap) {
            if (mapPositionW < (src->sh.outgassingMapWidth - 1)) {
                //Somewhere in the middle of the facet
                u = ((double) mapPositionW + currentParticle.randomGenerator.rnd()) / src->outgassingMapWidthD;
            } else {
                //Last element, prevent from going out of facet
                u = ((double) mapPositionW + currentParticle.randomGenerator.rnd() * (src->outgassingMapWidthD - (src->sh.outgassingMapWidth - 1))) /
                    src->outgassingMapWidthD;
            }
            if (mapPositionH < (src->sh.outgassingMapHeight - 1)) {
                //Somewhere in the middle of the facet
                v = ((double) mapPositionH + currentParticle.randomGenerator.rnd()) / src->outgassingMapHeightD;
            } else {
                //Last element, prevent from going out of facet
                v = ((double) mapPositionH + currentParticle.randomGenerator.rnd() * (src->outgassingMapHeightD - (src->sh.outgassingMapHeight - 1))) /
                    src->outgassingMapHeightD;
            }
        } else {
            u = currentParticle.randomGenerator.rnd();
            v = currentParticle.randomGenerator.rnd();
        }
        if (IsInFacet(*src, u, v)) {

            // (U,V) -> (x,y,z)
            currentParticle.position = src->sh.O + u * src->sh.U + v * src->sh.V;
            src->colU = u;
            src->colV = v;
            found = true;

        }
        nbTry++;
    }

    if (!found) {
        // Get the center, if the center is not included in the facet, a leak is generated.
        if (foundInMap) {
            //double uLength = sqrt(pow(src->sh.U.x, 2) + pow(src->sh.U.y, 2) + pow(src->sh.U.z, 2));
            //double vLength = sqrt(pow(src->sh.V.x, 2) + pow(src->sh.V.y, 2) + pow(src->sh.V.z, 2));
            double u = ((double) mapPositionW + 0.5) / src->outgassingMapWidthD;
            double v = ((double) mapPositionH + 0.5) / src->outgassingMapHeightD;
            currentParticle.position = src->sh.O + u * src->sh.U + v * src->sh.V;
            src->colU = u;
            src->colV = v;
        } else {
            src->colU = 0.5;
            src->colV = 0.5;
            currentParticle.position = model.structures[j].facets[i].sh.center;
        }

    }

    if (src->sh.isMoving && model.wp.motionType) RecordHit(HIT_MOVING, currentParticle);
    else RecordHit(HIT_DES, currentParticle); //create blue hit point for created particle

    //See docs/theta_gen.png for further details on angular distribution generation
    switch (src->sh.desorbType) {
        case DES_UNIFORM:
            currentParticle.direction = PolarToCartesian(src, std::acos(currentParticle.randomGenerator.rnd()), currentParticle.randomGenerator.rnd() * 2.0 * PI, reverse);
            break;
        case DES_NONE: //for file-based
        case DES_COSINE:
            currentParticle.direction = PolarToCartesian(src, std::acos(std::sqrt(currentParticle.randomGenerator.rnd())), currentParticle.randomGenerator.rnd() * 2.0 * PI,
                                                          reverse);
            break;
        case DES_COSINE_N:
            currentParticle.direction = PolarToCartesian(src, std::acos(
                    std::pow(currentParticle.randomGenerator.rnd(), 1.0 / (src->sh.desorbTypeN + 1.0))), currentParticle.randomGenerator.rnd() * 2.0 * PI, reverse);
            break;
        case DES_ANGLEMAP: {
            auto[theta, thetaLowerIndex, thetaOvershoot] = AnglemapGeneration::GenerateThetaFromAngleMap(
                    src->sh.anglemapParams, src->angleMap, currentParticle.randomGenerator.rnd());
            GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
            auto phi = AnglemapGeneration::GeneratePhiFromAngleMap(thetaLowerIndex, thetaOvershoot,
                                                                   src->sh.anglemapParams, src->angleMap,
                                                                   tmpResults.facetStates[src->globalId].recordedAngleMapPdf,
                                                                   currentParticle.randomGenerator.rnd());
            currentParticle.direction = PolarToCartesian(src, PI - theta, phi, false); //angle map contains incident angle (between N and source dir) and theta is dir (between N and dest dir)

        }
    }

    // Current structure
    if (src->sh.superIdx == -1) {
        std::ostringstream out;
        out << "Facet " << (src->globalId + 1) << " is in all structures, it shouldn't desorb.";
        //SetErrorSub(out.str().c_str());
        std::cerr << out.str() << std::endl;

        return false;
    }
    currentParticle.structureId = src->sh.superIdx;
    currentParticle.teleportedFrom = -1;

    // Count

    src->isHit = true;
    totalDesorbed++;
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    tmpResults.globalHits.globalHits.hit.nbDesorbed++;
    //nbPHit = 0;

    if (src->sh.isMoving) {
        TreatMovingFacet(currentParticle);
    }

    double ortVelocity =
            currentParticle.velocity * std::abs(Dot(currentParticle.direction, src->sh.N));
    /*src->sh.tmpCounter.hit.nbDesorbed++;
    src->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / ortVelocity; //was 2.0 / ortV
    src->sh.tmpCounter.hit.sum_v_ort += (model.wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
    IncreaseFacetCounter(src, currentParticle.flightTime, 0, 1, 0, 2.0 / ortVelocity,
                         (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity, currentParticle);
    //Desorption doesn't contribute to angular profiles, nor to angle maps
    ProfileFacet(src, currentParticle.flightTime, false, 2.0, 1.0, currentParticle); //was 2.0, 1.0
    LogHit(src, currentParticle);
    if (/*src->texture && */src->sh.countDes)
        RecordHitOnTexture(src, currentParticle.flightTime, true, 2.0, 1.0, currentParticle); //was 2.0, 1.0
    //if (src->direction && src->sh.countDirection) RecordDirectionVector(src, currentParticle.flightTime);

    // Reset volatile state
    if (hasVolatile) {
        for (auto &s : model.structures) {
            for (auto &f : s.facets) {
                f.isReady = true;
            }
        }
    }

    found = false;
    return true;
}

std::tuple<double, int, double>
AnglemapGeneration::GenerateThetaFromAngleMap(const AnglemapParams &anglemapParams, const Anglemap &anglemap,
                                              const double lookupValue) {
    //double lookupValue = currentParticle.randomGenerator.rnd();
    int thetaLowerIndex = my_lower_bound(lookupValue,
                                         anglemap.theta_CDF); //returns line number AFTER WHICH LINE lookup value resides in ( -1 .. size-2 )
    double theta, thetaOvershoot;

    if (thetaLowerIndex == -1) { //first half section
        thetaOvershoot = 0.5 + 0.5 * lookupValue / anglemap.theta_CDF[0]; //between 0.5 and 1
        theta = GetTheta((double) thetaLowerIndex + 0.5 + thetaOvershoot,
                         anglemapParams); //between 0 and the first section end
        return {theta, thetaLowerIndex, thetaOvershoot};
    } else if (thetaLowerIndex == (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes -
                                   1)) { //last half section //can this happen?
        thetaOvershoot = 0.5 * (lookupValue - anglemap.theta_CDF[thetaLowerIndex])
                         / (1.0 - anglemap.theta_CDF[thetaLowerIndex]); //between 0 and 0.5
        theta = GetTheta((double) thetaLowerIndex + 0.5 + thetaOvershoot,
                         anglemapParams); //between 0 and the first section end
        return {theta, thetaLowerIndex, thetaOvershoot};
    } else { //regular section
        if (/*true || */anglemap.phi_CDFsums[thetaLowerIndex] == anglemap.phi_CDFsums[thetaLowerIndex + 1]) {
            //The pdf's slope is 0, linear interpolation
            thetaOvershoot = (lookupValue - anglemap.theta_CDF[thetaLowerIndex]) /
                    (anglemap.theta_CDF[thetaLowerIndex + 1] - anglemap.theta_CDF[thetaLowerIndex]);
            theta = GetTheta((double) thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
        } else {
            //2nd degree interpolation
            // y(x) = ax^2 + bx + c
            // c: CDF value at lower index
            // b: pdf value at lower index
            // a: pdf slope at lower index / 2
            // dy := y - c
            // dx := x - [x at lower index]
            // dy = ax^2 + bx
            // dx = ( -b + sqrt(b^2 +4*a*dy) ) / (2a)
            double thetaStep = GetTheta((double) thetaLowerIndex + 1.5, anglemapParams) -
                               GetTheta((double) thetaLowerIndex + 0.5, anglemapParams);
            double c = anglemap.theta_CDF[thetaLowerIndex]; //CDF value at lower index
            double b = (double) anglemap.phi_CDFsums[thetaLowerIndex] / (double) anglemap.theta_CDFsum /
                       thetaStep; //pdf value at lower index
            double a = 0.5 * ((double) (anglemap.phi_CDFsums[thetaLowerIndex + 1]) - (double) anglemap.phi_CDFsums[thetaLowerIndex]) /
                       (double) anglemap.theta_CDFsum / Sqr(thetaStep); //pdf slope at lower index
            double dy = lookupValue - c;

            double dx = (-b + sqrt(Sqr(b) + 4 * a * dy)) /
                        (2 * a); //Since b>=0 it's the + branch of the +- that is valid for us

            thetaOvershoot = dx / thetaStep;
            theta = GetTheta((double) thetaLowerIndex + 0.5 + thetaOvershoot, anglemapParams);
        }
    }
    return {theta, thetaLowerIndex, thetaOvershoot};
}

/**
* \brief Generates phi angle (azimuth) from angle map
* \param thetaLowerIndex lower index of theta angle of bin in CDF (cummulative distribution function)
* \param thetaOvershoot corresponding to a weight of the previous and next lines
* \param anglemapParams parameters of the angle map
* \param randomGenerator reference to the random number generator (Mersenne Twister)
* \return phi angle
*/

double AnglemapGeneration::GeneratePhiFromAngleMap(const int &thetaLowerIndex, const double &thetaOvershoot,
                                                   const AnglemapParams &anglemapParams, Anglemap &anglemap,
                                                   const std::vector<size_t> &angleMapPDF,
                                                   double lookupValue) {
    //double lookupValue = currentParticle.randomGenerator.rnd();
    if (anglemapParams.phiWidth == 1) return -PI + 2.0 * PI * lookupValue; //special case, uniform phi distribution
    int phiLowerIndex;
    double weigh; //0: take previous theta line, 1: take next theta line, 0..1: interpolate in-between
    if (thetaLowerIndex == -1) { //first theta half section
        lookupValue += anglemap.phi_CDFs[0]; //periodic BCs over -PI...PI, can be larger than 1
        phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs[0],
                                       anglemapParams.phiWidth); //take entirely the phi ditro belonging to first theta
        weigh = thetaOvershoot; // [0.5 - 1], will subtract 0.5 when evaluating thetaIndex
    } else if (thetaLowerIndex ==
               (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)) { //last theta half section
        lookupValue += anglemap.phi_CDFs[thetaLowerIndex *
                                anglemapParams.phiWidth]; //periodic BCs over -PI...PI, can be larger than 1
        phiLowerIndex = my_lower_bound(lookupValue, &anglemap.phi_CDFs[thetaLowerIndex * anglemapParams.phiWidth],
                                       anglemapParams.phiWidth); //take entirely the phi ditro belonging to latest theta
        weigh = thetaOvershoot; // [0 - 0.5], will add 0.5 when evaluating thetaIndex
    } else {
        //Here we do a weighing both by the hit sum of the previous and next lines (w1 and w2) and also the weighs of the two lines based on thetaOvershoot (w3 and w4)
        // w1: sum of hits in previous line
        // w2: sum of hits in next line
        // w3: weigh of previous line (1 - thetaOvershoot)
        // w4: weigh of next line     (thetaOvershoot)
        // result: previous value weight: w1*w3 / (w1*w3 + w2*w4)
        //         next     value weight: w2*w4 / (w1*w3 + w2*w4) <- this will be the input for weighed_lower_bound

        double div;
        div = ((double) anglemap.phi_CDFsums[thetaLowerIndex] * (1.0 - thetaOvershoot) +
               (double) anglemap.phi_CDFsums[thetaLowerIndex + 1] * thetaOvershoot); // (w1*w3 + w2*w4)
        if (div > 0.0) {
            weigh = (thetaOvershoot * (double) anglemap.phi_CDFsums[thetaLowerIndex + 1]) /
                    div;    //      w2*w4 / (w1*w3 + w2*w4)
        } else {
            weigh = thetaOvershoot;
        }
        lookupValue += Weigh((double) anglemap.phi_CDFs[thetaLowerIndex * anglemapParams.phiWidth],
                             (double) anglemap.phi_CDFs[(thetaLowerIndex + 1) * anglemapParams.phiWidth], weigh);
        phiLowerIndex = weighed_lower_bound_X(lookupValue, weigh, &anglemap.phi_CDFs[thetaLowerIndex * anglemapParams.phiWidth],
                                              &anglemap.phi_CDFs[(thetaLowerIndex + 1) * anglemapParams.phiWidth],
                                              anglemapParams.phiWidth);
    }

    double phi, phiOvershoot;
    double thetaIndex = (double) thetaLowerIndex + 0.5 + weigh;
    if (phiLowerIndex == -1) { //first half section
        DEBUG_BREAK; //should not happen since we shifted the lookup value with first value
        phiOvershoot = 0.5 + 0.5 * lookupValue / GetPhiCDFValue(thetaIndex, 0, anglemapParams, anglemap); //between 0.5 and 1
        phi = GetPhi((double) phiLowerIndex + 0.5 + phiOvershoot, anglemapParams); //between 0 and the first section end
    }
        /*else if (phiLowerIndex == (anglemapParams.phiWidth - 1)) { //last half section
            phiOvershoot = 0.5 * (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams) )
                / (1.0 - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams)); //between 0 and 0.5
            phi = GetPhi((double)phiLowerIndex + 0.5 + phiOvershoot, anglemapParams); //between 0 and the first section end
        }*/
    else { //regular or last section
        if (/*true ||*/ GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams, angleMapPDF) ==
                        GetPhipdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams, angleMapPDF)) {
            //The pdf's slope is 0, linear interpolation
            phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap))
                           / (GetPhiCDFValue(thetaIndex, phiLowerIndex + 1, anglemapParams, anglemap) -
                              GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap));
            phi = GetPhi((double) phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
        } else {

            //2nd degree interpolation
            // y(x) = ax^2 + bx + c
            // c: CDF value at lower index
            // b: pdf value at lower index
            // a: pdf slope at lower index / 2
            // dy := y - c
            // dx := x - [x at lower index]
            // dy = ax^2 + bx
            // dx = ( -b + sqrt(b^2 +4*a*dy) ) / (2a)
            double phiStep = 2.0 * PI / (double) anglemapParams.phiWidth;
            double c = GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap); //CDF value at lower index
            double b = GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams, angleMapPDF) /
                    GetPhiCDFSum(thetaIndex, anglemapParams, anglemap) / phiStep; //pdf value at lower index
            double a = 0.5 * (GetPhipdfValue(thetaIndex, phiLowerIndex + 1, anglemapParams, angleMapPDF) -
                              GetPhipdfValue(thetaIndex, phiLowerIndex, anglemapParams, angleMapPDF)) /
                    GetPhiCDFSum(thetaIndex, anglemapParams, anglemap) / Sqr(phiStep); //pdf slope at lower index
            double dy = lookupValue - c;

            double D = Sqr(b) + 4 * a *
                                dy; //Discriminant. In rare cases it might be slightly negative, then fall back to linear interpolation:
            if (D < 0) {
                phiOvershoot = (lookupValue - GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams,
                                                             anglemap))
                               / (GetPhiCDFValue(thetaIndex, (int) IDX(phiLowerIndex + 1, anglemapParams.phiWidth),
                                                 anglemapParams, anglemap) -
                        GetPhiCDFValue(thetaIndex, phiLowerIndex, anglemapParams, anglemap));
            } else {
                double dx = (-b + sqrt(Sqr(b) + 4 * a * dy)) /
                            (2 * a); //Since b>=0 it's the + branch of the +- that is valid for us
                phiOvershoot = dx / phiStep;
            }
            phi = GetPhi((double) phiLowerIndex + 0.5 + phiOvershoot, anglemapParams);
        }
    }
    assert(phi > -PI && phi < PI);
    return phi;
}

/**
* \brief Converts from index to theta value
* \param thetaIndex theta index 
* \param anglemapParams parameters of the angle map
* \return theta angle
*/

double AnglemapGeneration::GetTheta(const double &thetaIndex, const AnglemapParams &anglemapParams) {
    if ((size_t) (thetaIndex) < anglemapParams.thetaLowerRes) { // 0 < theta < limit
        return anglemapParams.thetaLimit * (thetaIndex) / (double) anglemapParams.thetaLowerRes;
    } else { // limit < theta < PI/2
        return anglemapParams.thetaLimit +
               (PI / 2.0 - anglemapParams.thetaLimit) * (thetaIndex - (double) anglemapParams.thetaLowerRes) /
               (double) anglemapParams.thetaHigherRes;
    }
}

/**
* \brief makes phiIndex circular and converts from index to -pi...pi
* \param phiIndex phi index
* \param anglemapParams parameters of the angle map
* \return phi angle
*/

double AnglemapGeneration::GetPhi(const double &phiIndex, const AnglemapParams &anglemapParams)
//makes phiIndex circular and converts from index to -pi...pi
{
    double width = (double) anglemapParams.phiWidth;
    double correctedIndex = (phiIndex < width) ? phiIndex : phiIndex - width;
    return -PI + 2.0 * PI * correctedIndex / width;
}

double
AnglemapGeneration::GetPhipdfValue(const double &thetaIndex, const int &phiLowerIndex,
                                   const AnglemapParams &anglemapParams, const std::vector<size_t> &angleMapPDF)
//phiLowerIndex is circularized
{
    if (thetaIndex < 0.5) {
        return (double) angleMapPDF[IDX(phiLowerIndex, anglemapParams.phiWidth)];
    } else if (thetaIndex > (double) (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
        return (double) angleMapPDF[
                anglemapParams.phiWidth * (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1) +
                IDX(phiLowerIndex, anglemapParams.phiWidth)];
    } else {
        size_t thetaLowerIndex = (size_t) (thetaIndex - 0.5);
        double thetaOvershoot = thetaIndex - 0.5 - (double) thetaLowerIndex;
        double valueFromLowerpdf = (double) angleMapPDF[anglemapParams.phiWidth * thetaLowerIndex +
                                                            IDX(phiLowerIndex, anglemapParams.phiWidth)];
        double valueFromHigherpdf = (double) angleMapPDF[anglemapParams.phiWidth * (thetaLowerIndex + 1) +
                                                             IDX(phiLowerIndex, anglemapParams.phiWidth)];
        return Weigh(valueFromLowerpdf, valueFromHigherpdf, thetaOvershoot);
    }
}

/**
* \brief Get phi value from cummulative density function
* \param thetaIndex theta index
* \param phiLowerIndex lower index of the bin
* \param anglemapParams parameters of the angle map
* \return phi cdf value
*/

double
AnglemapGeneration::GetPhiCDFValue(const double &thetaIndex, const int &phiLowerIndex,
                                   const AnglemapParams &anglemapParams, const Anglemap &anglemap) {
    if (thetaIndex < 0.5) {
        return (phiLowerIndex < anglemapParams.phiWidth) ? anglemap.phi_CDFs[phiLowerIndex] : 1.0 + anglemap.phi_CDFs[0];
    } else if (thetaIndex > (double) (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
        return (phiLowerIndex < anglemapParams.phiWidth) ? anglemap.phi_CDFs[
                anglemapParams.phiWidth * (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1) +
                phiLowerIndex] : 1.0 + anglemap.phi_CDFs[anglemapParams.phiWidth *
                                                (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1)];
    } else {
        size_t thetaLowerIndex = (size_t) (thetaIndex - 0.5);
        double thetaOvershoot = thetaIndex - 0.5 - (double) thetaLowerIndex;
        double valueFromLowerCDF = (phiLowerIndex < anglemapParams.phiWidth) ? anglemap.phi_CDFs[
                anglemapParams.phiWidth * thetaLowerIndex + phiLowerIndex] : 1.0 + anglemap.phi_CDFs[anglemapParams.phiWidth *
                                                                                            (thetaLowerIndex)];
        double valueFromHigherCDF = (phiLowerIndex < anglemapParams.phiWidth) ? anglemap.phi_CDFs[
                anglemapParams.phiWidth * (thetaLowerIndex + 1) + phiLowerIndex] : 1.0 +
                anglemap.phi_CDFs[anglemapParams.phiWidth *
                                                                                            (thetaLowerIndex + 1)];
        return Weigh(valueFromLowerCDF, valueFromHigherCDF, thetaOvershoot);
    }

}

/**
* \brief Get phi value from cummulative density function
* \param thetaIndex theta index
* \param anglemapParams parameters of the angle map
* \return phi cdf summed value
*/

double AnglemapGeneration::GetPhiCDFSum(const double &thetaIndex, const AnglemapParams &anglemapParams,
                                        const Anglemap &anglemap) {
    if (thetaIndex < 0.5) {
        return (double) anglemap.phi_CDFsums[0];
    } else if (thetaIndex > (double) (anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes) - 0.5) {
        return (double) anglemap.phi_CDFsums[anglemapParams.thetaLowerRes + anglemapParams.thetaHigherRes - 1];
    } else {
        size_t thetaLowerIndex = (size_t) (thetaIndex - 0.5);
        double thetaOvershoot = thetaIndex - 0.5 - (double) thetaLowerIndex;
        double valueFromLowerSum = (double) anglemap.phi_CDFsums[thetaLowerIndex];
        double valueFromHigherSum = (double) anglemap.phi_CDFsums[thetaLowerIndex + 1];
        return Weigh(valueFromLowerSum, valueFromHigherSum, thetaOvershoot);
    }
}

/**
* \brief Perform a bounce from a facet by logging the hit and sometimes relaunching it
* \param iFacet facet corresponding to the bounce event
*/

void Simulation::PerformBounce(SubprocessFacet *iFacet, CurrentParticleStatus &currentParticle) {

    bool revert = false;
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    tmpResults.globalHits.globalHits.hit.nbMCHit++; //global
    tmpResults.globalHits.globalHits.hit.nbHitEquiv += currentParticle.oriRatio;

    // Handle super structure link facet. Can be
    if (iFacet->sh.superDest) {
        IncreaseFacetCounter(iFacet, currentParticle.flightTime, 1, 0, 0, 0, 0, currentParticle);
        currentParticle.structureId = iFacet->sh.superDest - 1;
        if (iFacet->sh.isMoving) { //A very special case where link facets can be used as transparent but moving facets
            RecordHit(HIT_MOVING, currentParticle);
            TreatMovingFacet(currentParticle);
        } else {
            // Count this hit as a transparent pass
            RecordHit(HIT_TRANS, currentParticle);
        }
        LogHit(iFacet, currentParticle);
        ProfileFacet(iFacet, currentParticle.flightTime, true, 2.0, 2.0, currentParticle);
        if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet, currentParticle);
        if (/*iFacet->texture &&*/ iFacet->sh.countTrans)
            RecordHitOnTexture(iFacet, currentParticle.flightTime, true, 2.0, 2.0, currentParticle);
        if (/*iFacet->direction &&*/ iFacet->sh.countDirection)
            RecordDirectionVector(iFacet, currentParticle.flightTime, currentParticle);

        return;

    }

    // Handle volatile facet
    if (iFacet->sh.isVolatile) {

        if (iFacet->isReady) {
            IncreaseFacetCounter(iFacet, currentParticle.flightTime, 0, 0, 1, 0, 0, currentParticle);
            iFacet->isReady = false;
            LogHit(iFacet, currentParticle);
            ProfileFacet(iFacet, currentParticle.flightTime, true, 2.0, 1.0, currentParticle);
            if (/*iFacet->texture && */iFacet->sh.countAbs)
                RecordHitOnTexture(iFacet, currentParticle.flightTime, true, 2.0, 1.0, currentParticle);
            if (/*iFacet->direction && */iFacet->sh.countDirection)
                RecordDirectionVector(iFacet, currentParticle.flightTime, currentParticle);
        }
        return;

    }

    if (iFacet->sh.is2sided) {
        // We may need to revert normal in case of 2 sided hit
        revert = Dot(currentParticle.direction, iFacet->sh.N) > 0.0;
    }

    //Texture/Profile incoming hit


    //Register (orthogonal) velocity
    double ortVelocity =
            currentParticle.velocity * std::abs(Dot(currentParticle.direction, iFacet->sh.N));

    /*iFacet->sh.tmpCounter.hit.nbMCHit++; //hit facet
    iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
    iFacet->sh.tmpCounter.hit.sum_v_ort += (model.wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/

    IncreaseFacetCounter(iFacet, currentParticle.flightTime, 1, 0, 0, 1.0 / ortVelocity,
                         (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity, currentParticle);
    currentParticle.nbBounces++;
    if (/*iFacet->texture &&*/ iFacet->sh.countRefl)
        RecordHitOnTexture(iFacet, currentParticle.flightTime, true, 1.0, 1.0, currentParticle);
    if (/*iFacet->direction &&*/ iFacet->sh.countDirection)
        RecordDirectionVector(iFacet, currentParticle.flightTime, currentParticle);
    LogHit(iFacet, currentParticle);
    ProfileFacet(iFacet, currentParticle.flightTime, true, 1.0, 1.0, currentParticle);
    if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet, currentParticle);

    // Relaunch particle
    UpdateVelocity(iFacet, currentParticle);
    //Sojourn time
    if (iFacet->sh.enableSojournTime) {
        double A = exp(-iFacet->sh.sojournE / (8.31 * iFacet->sh.temperature));
        currentParticle.flightTime += -log(currentParticle.randomGenerator.rnd()) / (A * iFacet->sh.sojournFreq);
    }

    if (iFacet->sh.reflection.diffusePart > 0.999999) { //Speedup branch for most common, diffuse case
        currentParticle.direction = PolarToCartesian(iFacet, std::acos(std::sqrt(currentParticle.randomGenerator.rnd())), currentParticle.randomGenerator.rnd() * 2.0 * PI,
                                                      revert);
    } else {
        double reflTypeRnd = currentParticle.randomGenerator.rnd();
        if (reflTypeRnd < iFacet->sh.reflection.diffusePart) {
            //diffuse reflection
            //See docs/theta_gen.png for further details on angular distribution generation
            currentParticle.direction = PolarToCartesian(iFacet, std::acos(std::sqrt(currentParticle.randomGenerator.rnd())), currentParticle.randomGenerator.rnd() * 2.0 * PI,
                                                          revert);
        } else if (reflTypeRnd < (iFacet->sh.reflection.diffusePart + iFacet->sh.reflection.specularPart)) {
            //specular reflection
            auto[inTheta, inPhi] = CartesianToPolar(currentParticle.direction, iFacet->sh.nU, iFacet->sh.nV,
                                                    iFacet->sh.N);
            currentParticle.direction = PolarToCartesian(iFacet, PI - inTheta, inPhi, false);

        } else {
            //Cos^N reflection
            currentParticle.direction = PolarToCartesian(iFacet, std::acos(
                    std::pow(currentParticle.randomGenerator.rnd(), 1.0 / (iFacet->sh.reflection.cosineExponent + 1.0))), currentParticle.randomGenerator.rnd() * 2.0 * PI, revert);
        }
    }

    if (iFacet->sh.isMoving) {
        TreatMovingFacet(currentParticle);
    }

    //Texture/Profile outgoing particle
    //Register outgoing velocity
    ortVelocity = currentParticle.velocity * std::abs(Dot(currentParticle.direction, iFacet->sh.N));

    /*iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 1.0 / ortVelocity;
    iFacet->sh.tmpCounter.hit.sum_v_ort += (model.wp.useMaxwellDistribution ? 1.0 : 1.1781)*ortVelocity;*/
    IncreaseFacetCounter(iFacet, currentParticle.flightTime, 0, 0, 0, 1.0 / ortVelocity,
                         (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity, currentParticle);
    if (/*iFacet->texture &&*/ iFacet->sh.countRefl)
        RecordHitOnTexture(iFacet, currentParticle.flightTime, false, 1.0,
                           1.0, currentParticle); //count again for outward velocity
    ProfileFacet(iFacet, currentParticle.flightTime, false, 1.0, 1.0, currentParticle);
    //no direction count on outgoing, neither angle map

    if (iFacet->sh.isMoving && model.wp.motionType) RecordHit(HIT_MOVING, currentParticle);
    else RecordHit(HIT_REF, currentParticle);
    currentParticle.lastHitFacet = iFacet;
    //nbPHit++;
}

void Simulation::PerformTransparentPass(SubprocessFacet *iFacet) { //disabled, caused finding hits with the same facet
    /*double directionFactor = abs(DOT3(
        currentParticle.direction.x, currentParticle.direction.y, currentParticle.direction.z,
        iFacet->sh.N.x, iFacet->sh.N.y, iFacet->sh.N.z));
    iFacet->sh.tmpCounter.hit.nbMCHit++;
    iFacet->sh.tmpCounter.hit.sum_1_per_ort_velocity += 2.0 / (currentParticle.velocity*directionFactor);
    iFacet->sh.tmpCounter.hit.sum_v_ort += 2.0*(model.wp.useMaxwellDistribution ? 1.0 : 1.1781)*currentParticle.velocity*directionFactor;
    iFacet->isHit = true;
    if (iFacet->texture && iFacet->sh.countTrans) RecordHitOnTexture(iFacet, currentParticle.flightTime + iFacet->colDist / 100.0 / currentParticle.velocity,
        true, 2.0, 2.0);
    if (iFacet->direction && iFacet->sh.countDirection) RecordDirectionVector(iFacet, currentParticle.flightTime + iFacet->colDist / 100.0 / currentParticle.velocity);
    ProfileFacet(iFacet, currentParticle.flightTime + iFacet->colDist / 100.0 / currentParticle.velocity,
        true, 2.0, 2.0);
    RecordHit(HIT_TRANS);
    lastHit = iFacet;*/
}

void Simulation::RecordAbsorb(SubprocessFacet *iFacet, CurrentParticleStatus &currentParticle) {
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    tmpResults.globalHits.globalHits.hit.nbMCHit++; //global
    tmpResults.globalHits.globalHits.hit.nbHitEquiv += currentParticle.oriRatio;
    tmpResults.globalHits.globalHits.hit.nbAbsEquiv += currentParticle.oriRatio;

    RecordHistograms(iFacet, currentParticle);

    RecordHit(HIT_ABS, currentParticle);
    double ortVelocity =
            currentParticle.velocity * std::abs(Dot(currentParticle.direction, iFacet->sh.N));
    IncreaseFacetCounter(iFacet, currentParticle.flightTime, 1, 0, 1, 2.0 / ortVelocity,
                         (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity, currentParticle);
    LogHit(iFacet, currentParticle);
    ProfileFacet(iFacet, currentParticle.flightTime, true, 2.0, 1.0, currentParticle); //was 2.0, 1.0
    if (iFacet->sh.anglemapParams.record) RecordAngleMap(iFacet, currentParticle);
    if (/*iFacet->texture &&*/ iFacet->sh.countAbs)
        RecordHitOnTexture(iFacet, currentParticle.flightTime, true, 2.0, 1.0, currentParticle); //was 2.0, 1.0
    if (/*iFacet->direction &&*/ iFacet->sh.countDirection)
        RecordDirectionVector(iFacet, currentParticle.flightTime, currentParticle);
}

void Simulation::RecordHistograms(SubprocessFacet *iFacet, CurrentParticleStatus &currentParticle) {
    //Record in global and facet histograms
    size_t binIndex;
    auto& tmpGlobalHistograms = tmpGlobalResults[omp_get_thread_num()].globalHistograms;

    {
        auto& facetHistogram = tmpGlobalResults[omp_get_thread_num()].facetStates[iFacet->globalId].momentResults[0].histogram;

        if (model.wp.globalHistogramParams.recordBounce) {
            binIndex = Min(currentParticle.nbBounces / model.wp.globalHistogramParams.nbBounceBinsize,
                           model.wp.globalHistogramParams.GetBounceHistogramSize() - 1);
            tmpGlobalHistograms[0].nbHitsHistogram[binIndex] += currentParticle.oriRatio;
        }
        if (model.wp.globalHistogramParams.recordDistance) {
            binIndex = Min(static_cast<size_t>(currentParticle.distanceTraveled /
                                               model.wp.globalHistogramParams.distanceBinsize),
                           model.wp.globalHistogramParams.GetDistanceHistogramSize() - 1);
            tmpGlobalHistograms[0].distanceHistogram[binIndex] += currentParticle.oriRatio;
        }
        if (model.wp.globalHistogramParams.recordTime) {
            binIndex = Min(static_cast<size_t>(currentParticle.flightTime /
                                               model.wp.globalHistogramParams.timeBinsize),
                           model.wp.globalHistogramParams.GetTimeHistogramSize() - 1);
            tmpGlobalHistograms[0].timeHistogram[binIndex] += currentParticle.oriRatio;
        }
        if (iFacet->sh.facetHistogramParams.recordBounce) {
            binIndex = Min(currentParticle.nbBounces / iFacet->sh.facetHistogramParams.nbBounceBinsize,
                           iFacet->sh.facetHistogramParams.GetBounceHistogramSize() - 1);
            facetHistogram.nbHitsHistogram[binIndex] += currentParticle.oriRatio;
        }
        if (iFacet->sh.facetHistogramParams.recordDistance) {
            binIndex = Min(static_cast<size_t>(currentParticle.distanceTraveled /
                                               iFacet->sh.facetHistogramParams.distanceBinsize),
                           iFacet->sh.facetHistogramParams.GetDistanceHistogramSize() - 1);
            facetHistogram.distanceHistogram[binIndex] += currentParticle.oriRatio;
        }
        if (iFacet->sh.facetHistogramParams.recordTime) {
            binIndex = Min(static_cast<size_t>(currentParticle.flightTime /
                                               iFacet->sh.facetHistogramParams.timeBinsize),
                           iFacet->sh.facetHistogramParams.GetTimeHistogramSize() - 1);
            facetHistogram.timeHistogram[binIndex] += currentParticle.oriRatio;
        }
    }

    {

        int m = -1;
        if ((m = LookupMomentIndex(currentParticle.flightTime, model.tdParams.moments,
                                   currentParticle.lastMomentIndex)) >= 0) {
            currentParticle.lastMomentIndex = m;
            auto& facetHistogram = tmpGlobalResults[omp_get_thread_num()].facetStates[iFacet->globalId].momentResults[m].histogram;

            if (model.wp.globalHistogramParams.recordBounce) {
                binIndex = Min(currentParticle.nbBounces / model.wp.globalHistogramParams.nbBounceBinsize,
                               model.wp.globalHistogramParams.GetBounceHistogramSize() - 1);
                tmpGlobalHistograms[m].nbHitsHistogram[binIndex] += currentParticle.oriRatio;
            }
            if (model.wp.globalHistogramParams.recordDistance) {
                binIndex = Min(static_cast<size_t>(currentParticle.distanceTraveled /
                                                   model.wp.globalHistogramParams.distanceBinsize),
                               model.wp.globalHistogramParams.GetDistanceHistogramSize() - 1);
                tmpGlobalHistograms[m].distanceHistogram[binIndex] += currentParticle.oriRatio;
            }
            if (model.wp.globalHistogramParams.recordTime) {
                binIndex = Min(static_cast<size_t>(currentParticle.flightTime /
                                                   model.wp.globalHistogramParams.timeBinsize),
                               model.wp.globalHistogramParams.GetTimeHistogramSize() - 1);
                tmpGlobalHistograms[m].timeHistogram[binIndex] += currentParticle.oriRatio;
            }
            if (iFacet->sh.facetHistogramParams.recordBounce) {
                binIndex = Min(currentParticle.nbBounces / iFacet->sh.facetHistogramParams.nbBounceBinsize,
                               iFacet->sh.facetHistogramParams.GetBounceHistogramSize() - 1);
                facetHistogram.nbHitsHistogram[binIndex] += currentParticle.oriRatio;
            }
            if (iFacet->sh.facetHistogramParams.recordDistance) {
                binIndex = Min(static_cast<size_t>(currentParticle.distanceTraveled /
                                                   iFacet->sh.facetHistogramParams.distanceBinsize),
                               iFacet->sh.facetHistogramParams.GetDistanceHistogramSize() - 1);
                facetHistogram.distanceHistogram[binIndex] += currentParticle.oriRatio;
            }
            if (iFacet->sh.facetHistogramParams.recordTime) {
                binIndex = Min(static_cast<size_t>(currentParticle.flightTime /
                                                   iFacet->sh.facetHistogramParams.timeBinsize),
                               iFacet->sh.facetHistogramParams.GetTimeHistogramSize() - 1);
                facetHistogram.timeHistogram[binIndex] += currentParticle.oriRatio;
            }
        }
    }
}

void Simulation::RecordHitOnTexture(SubprocessFacet *f, double time, bool countHit, double velocity_factor,
                                    double ortSpeedFactor, CurrentParticleStatus &currentParticle) {

    size_t tu = (size_t)(f->colU * f->sh.texWidthD);
    size_t tv = (size_t)(f->colV * f->sh.texHeightD);
    size_t add = tu + tv * (f->sh.texWidth);
    double ortVelocity = (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * currentParticle.velocity *
                         std::abs(Dot(currentParticle.direction,
                                      f->sh.N)); //surface-orthogonal velocity component

    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    {
        TextureCell &texture = tmpResults.facetStates[f->globalId].momentResults[0].texture[add];
        if (countHit) texture.countEquiv += currentParticle.oriRatio;
        texture.sum_1_per_ort_velocity +=
                currentParticle.oriRatio * velocity_factor / ortVelocity;
        texture.sum_v_ort_per_area += currentParticle.oriRatio * ortSpeedFactor * ortVelocity *
                                      f->textureCellIncrements[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
    }
    int m = -1;
    if((m = LookupMomentIndex(time, model.tdParams.moments, currentParticle.lastMomentIndex)) >= 0){
        currentParticle.lastMomentIndex = m;
        TextureCell& texture = tmpResults.facetStates[f->globalId].momentResults[m].texture[add];
        if (countHit) texture.countEquiv += currentParticle.oriRatio;
        texture.sum_1_per_ort_velocity +=
                currentParticle.oriRatio * velocity_factor / ortVelocity;
        texture.sum_v_ort_per_area += currentParticle.oriRatio * ortSpeedFactor * ortVelocity *
                                      f->textureCellIncrements[add]; // sum ortho_velocity[m/s] / cell_area[cm2]
    }
}

void Simulation::RecordDirectionVector(SubprocessFacet *f, double time, CurrentParticleStatus &currentParticle) {
    size_t tu = (size_t)(f->colU * f->sh.texWidthD);
    size_t tv = (size_t)(f->colV * f->sh.texHeightD);
    size_t add = tu + tv * (f->sh.texWidth);

    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    {
        DirectionCell &direction = tmpResults.facetStates[f->globalId].momentResults[0].direction[add];
        direction.dir += currentParticle.oriRatio * currentParticle.direction * currentParticle.velocity;
        direction.count++;
    }
    int m = -1;
    if((m = LookupMomentIndex(time, model.tdParams.moments, currentParticle.lastMomentIndex)) >= 0){
        currentParticle.lastMomentIndex = m;
        DirectionCell &direction = tmpResults.facetStates[f->globalId].momentResults[m].direction[add];
        direction.dir += currentParticle.oriRatio * currentParticle.direction * currentParticle.velocity;
        direction.count++;
    }

}

void
Simulation::ProfileFacet(SubprocessFacet *f, double time, bool countHit, double velocity_factor, double ortSpeedFactor,
                         CurrentParticleStatus &currentParticle) {

    size_t nbMoments = model.tdParams.moments.size();
    int m = LookupMomentIndex(time, model.tdParams.moments, currentParticle.lastMomentIndex);
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    if (countHit && f->sh.profileType == PROFILE_ANGULAR) {
        double dot = Dot(f->sh.N, currentParticle.direction);
        double theta = std::acos(std::abs(dot));     // Angle to normal (PI/2 => PI)
        size_t pos = (size_t)(theta / (PI / 2) * ((double) PROFILE_SIZE)); // To Grad
        Saturate(pos, 0, PROFILE_SIZE - 1);

        tmpResults.facetStates[f->globalId].momentResults[0].profile[pos].countEquiv += currentParticle.oriRatio;
        if(m >= 0){
            currentParticle.lastMomentIndex = m;
            tmpResults.facetStates[f->globalId].momentResults[m].profile[pos].countEquiv += currentParticle.oriRatio;
        }
    } else if (f->sh.profileType == PROFILE_U || f->sh.profileType == PROFILE_V) {
        size_t pos = (size_t)((f->sh.profileType == PROFILE_U ? f->colU : f->colV) * (double) PROFILE_SIZE);
        if (pos >= 0 && pos < PROFILE_SIZE) {
            {
                ProfileSlice &profile = tmpResults.facetStates[f->globalId].momentResults[0].profile[pos];
                if (countHit) profile.countEquiv += currentParticle.oriRatio;
                double ortVelocity = currentParticle.velocity *
                                     std::abs(Dot(f->sh.N, currentParticle.direction));
                profile.sum_1_per_ort_velocity +=
                        currentParticle.oriRatio * velocity_factor / ortVelocity;
                profile.sum_v_ort += currentParticle.oriRatio * ortSpeedFactor *
                                     (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity;
            }
            if(m >= 0) {
                currentParticle.lastMomentIndex = m;
                ProfileSlice &profile = tmpResults.facetStates[f->globalId].momentResults[m].profile[pos];
                if (countHit) profile.countEquiv += currentParticle.oriRatio;
                double ortVelocity = currentParticle.velocity *
                                     std::abs(Dot(f->sh.N, currentParticle.direction));
                profile.sum_1_per_ort_velocity +=
                        currentParticle.oriRatio * velocity_factor / ortVelocity;
                profile.sum_v_ort += currentParticle.oriRatio * ortSpeedFactor *
                                     (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * ortVelocity;
            }
        }
    } else if (countHit && (f->sh.profileType == PROFILE_VELOCITY || f->sh.profileType == PROFILE_ORT_VELOCITY ||
                            f->sh.profileType == PROFILE_TAN_VELOCITY)) {
        double dot;
        if (f->sh.profileType == PROFILE_VELOCITY) {
            dot = 1.0;
        } else if (f->sh.profileType == PROFILE_ORT_VELOCITY) {
            dot = std::abs(Dot(f->sh.N, currentParticle.direction));  //cos(theta) as "dot" value
        } else { //Tangential
            dot = std::sqrt(1 - Sqr(std::abs(Dot(f->sh.N, currentParticle.direction))));  //tangential
        }
        size_t pos = (size_t)(dot * currentParticle.velocity / f->sh.maxSpeed *
                              (double) PROFILE_SIZE); //"dot" default value is 1.0
        if (pos >= 0 && pos < PROFILE_SIZE) {
            tmpResults.facetStates[f->globalId].momentResults[0].profile[pos].countEquiv += currentParticle.oriRatio;
            if (m >= 0) {
                currentParticle.lastMomentIndex = m;
                tmpResults.facetStates[f->globalId].momentResults[m].profile[pos].countEquiv += currentParticle.oriRatio;
            }
        }
    }
}

void Simulation::LogHit(SubprocessFacet *f, CurrentParticleStatus &currentParticle) {
    if(omp_get_thread_num() != 0) return; // only let 1 thread update
    if (model.otfParams.enableLogging &&
        model.otfParams.logFacetId == f->globalId &&
        tmpParticleLog.size() < (model.otfParams.logLimit / model.otfParams.nbProcess)) {
        ParticleLoggerItem log;
        log.facetHitPosition = Vector2d(f->colU, f->colV);
        std::tie(log.hitTheta, log.hitPhi) = CartesianToPolar(currentParticle.direction, f->sh.nU, f->sh.nV,
                                                              f->sh.N);
        log.oriRatio = currentParticle.oriRatio;
        log.particleDecayMoment = currentParticle.expectedDecayMoment;
        log.time = currentParticle.flightTime;
        log.velocity = currentParticle.velocity;
        tmpParticleLog.push_back(log);
    }
}

void Simulation::RecordAngleMap(SubprocessFacet *collidedFacet, CurrentParticleStatus &currentParticle) {
    auto[inTheta, inPhi] = CartesianToPolar(currentParticle.direction, collidedFacet->sh.nU,
                                            collidedFacet->sh.nV, collidedFacet->sh.N);
    if (inTheta > PI / 2.0)
        inTheta = std::abs(
                PI - inTheta); //theta is originally respective to N, but we'd like the angle between 0 and PI/2
    bool countTheta = true;
    size_t thetaIndex;
    if (inTheta < collidedFacet->sh.anglemapParams.thetaLimit) {
        if (collidedFacet->sh.anglemapParams.thetaLowerRes > 0) {
            thetaIndex = (size_t)(inTheta / collidedFacet->sh.anglemapParams.thetaLimit *
                                  (double) collidedFacet->sh.anglemapParams.thetaLowerRes);
        } else {
            countTheta = false;
        }
    } else {
        if (collidedFacet->sh.anglemapParams.thetaHigherRes > 0) {
            thetaIndex = collidedFacet->sh.anglemapParams.thetaLowerRes +
                         (size_t)((inTheta - collidedFacet->sh.anglemapParams.thetaLimit)
                                  / (PI / 2.0 - collidedFacet->sh.anglemapParams.thetaLimit) *
                                  (double) collidedFacet->sh.anglemapParams.thetaHigherRes);
        } else {
            countTheta = false;
        }
    }
    if (countTheta) {
        size_t phiIndex = (size_t)((inPhi + 3.1415926) / (2.0 * PI) *
                                   (double) collidedFacet->sh.anglemapParams.phiWidth); //Phi: -PI..PI , and shifting by a number slightly smaller than PI to store on interval [0,2PI[

        GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
        auto& angleMap = tmpResults.facetStates[collidedFacet->globalId].recordedAngleMapPdf;
        angleMap[thetaIndex * collidedFacet->sh.anglemapParams.phiWidth + phiIndex]++;
    }
}

void Simulation::UpdateVelocity(SubprocessFacet *collidedFacet, CurrentParticleStatus &currentParticle) {
    if (collidedFacet->sh.accomodationFactor > 0.9999) { //speedup for the most common case: perfect thermalization
        if (model.wp.useMaxwellDistribution)
            currentParticle.velocity = GenerateRandomVelocity(collidedFacet->sh.CDFid, currentParticle.randomGenerator.rnd());
        else
            currentParticle.velocity =
                    145.469 * std::sqrt(collidedFacet->sh.temperature / model.wp.gasMass);
    } else {
        double oldSpeed2 = pow(currentParticle.velocity, 2);
        double newSpeed2;
        if (model.wp.useMaxwellDistribution) newSpeed2 = pow(GenerateRandomVelocity(collidedFacet->sh.CDFid, currentParticle.randomGenerator.rnd()), 2);
        else newSpeed2 = /*145.469*/ 29369.939 * (collidedFacet->sh.temperature / model.wp.gasMass);
        //sqrt(29369)=171.3766= sqrt(8*R*1000/PI)*3PI/8, that is, the constant part of the v_avg=sqrt(8RT/PI/m/0.001)) found in literature, multiplied by
        //the corrective factor of 3PI/8 that accounts for moving from volumetric speed distribution to wall collision speed distribution
        currentParticle.velocity = std::sqrt(
                oldSpeed2 + (newSpeed2 - oldSpeed2) * collidedFacet->sh.accomodationFactor);
    }
}

double Simulation::GenerateRandomVelocity(int CDFId, const double rndVal) {
    //return FastLookupY(currentParticle.randomGenerator.rnd(),CDFs[CDFId],false);
    //double r = currentParticle.randomGenerator.rnd();
    double v = InterpolateX(rndVal, model.tdParams.CDFs[CDFId], false, true); //Allow extrapolate
    return v;
}

double Simulation::GenerateDesorptionTime(SubprocessFacet *src, const double rndVal) {
    if (src->sh.outgassing_paramId >= 0) { //time-dependent desorption
        return InterpolateX(rndVal * model.tdParams.IDs[src->sh.IDid].back().second, model.tdParams.IDs[src->sh.IDid], false,
                            true); //allow extrapolate
    } else {
        return rndVal * model.wp.latestMoment; //continous desorption between 0 and latestMoment
    }
}

double Simulation::GetStickingAt(SubprocessFacet *f, double time) {
    if (f->sh.sticking_paramId == -1) //constant sticking
        return f->sh.sticking;
    else return model.tdParams.parameters[f->sh.sticking_paramId].InterpolateY(time, false);
}

double Simulation::GetOpacityAt(SubprocessFacet *f, double time) {
    if (f->sh.opacity_paramId == -1) //constant sticking
        return f->sh.opacity;
    else return model.tdParams.parameters[f->sh.opacity_paramId].InterpolateY(time, false);
}

/**
* \brief Updates particle direction and velocity if we are dealing with a moving facet (translated or rotated)
*/


void Simulation::TreatMovingFacet(CurrentParticleStatus &currentParticle) {
    Vector3d localVelocityToAdd;
    if (model.wp.motionType == 1) { //Translation
        localVelocityToAdd = model.wp.motionVector2; //Fixed translational vector
    } else if (model.wp.motionType == 2) { //Rotation
        Vector3d distanceVector = 0.01 * (currentParticle.position -
                                          model.wp.motionVector1); //distance from base, with cm->m conversion, motionVector1 is rotation base point
        localVelocityToAdd = CrossProduct(model.wp.motionVector2, distanceVector); //motionVector2 is rotation axis
    }
    Vector3d oldVelocity, newVelocity;
    oldVelocity = currentParticle.direction * currentParticle.velocity;
    newVelocity = oldVelocity + localVelocityToAdd;
    currentParticle.direction = newVelocity.Normalized();
    currentParticle.velocity = newVelocity.Norme();
}

/**
* \brief Increase facet counter on a hit, pass etc.
* \param f source facet
* \param time simulation time
* \param hit amount of hits to add
* \param desorb amount of desorptions to add
* \param absorb amount of absorptions to add
* \param sum_1_per_v reciprocals of orthogonal speed components to add
* \param sum_v_ort orthogonal momentum change to add
*/
void Simulation::IncreaseFacetCounter(SubprocessFacet *f, double time, size_t hit, size_t desorb, size_t absorb,
                                      double sum_1_per_v, double sum_v_ort, CurrentParticleStatus &currentParticle) {
    const double hitEquiv = static_cast<double>(hit) * currentParticle.oriRatio;
    GlobalSimuState& tmpResults = tmpGlobalResults[omp_get_thread_num()];
    {
        FacetHitBuffer &hits = tmpResults.facetStates[f->globalId].momentResults[0].hits;
        hits.hit.nbMCHit += hit;
        hits.hit.nbHitEquiv += hitEquiv;
        hits.hit.nbDesorbed += desorb;
        hits.hit.nbAbsEquiv += static_cast<double>(absorb) * currentParticle.oriRatio;
        hits.hit.sum_1_per_ort_velocity += currentParticle.oriRatio * sum_1_per_v;
        hits.hit.sum_v_ort += currentParticle.oriRatio * sum_v_ort;
        hits.hit.sum_1_per_velocity += (hitEquiv + static_cast<double>(desorb)) / currentParticle.velocity;
    }
    int m = -1;
    if((m = LookupMomentIndex(time, model.tdParams.moments, currentParticle.lastMomentIndex)) >= 0){
        currentParticle.lastMomentIndex = m;
        FacetHitBuffer& hits = tmpResults.facetStates[f->globalId].momentResults[m].hits;
        hits.hit.nbMCHit += hit;
        hits.hit.nbHitEquiv += hitEquiv;
        hits.hit.nbDesorbed += desorb;
        hits.hit.nbAbsEquiv += static_cast<double>(absorb) * currentParticle.oriRatio;
        hits.hit.sum_1_per_ort_velocity += currentParticle.oriRatio * sum_1_per_v;
        hits.hit.sum_v_ort += currentParticle.oriRatio * sum_v_ort;
        hits.hit.sum_1_per_velocity += (hitEquiv + static_cast<double>(desorb)) / currentParticle.velocity;
    }
}

void Simulation::RegisterTransparentPass(SubprocessFacet *facet, CurrentParticleStatus &currentParticle) {
    double directionFactor = std::abs(Dot(currentParticle.direction, facet->sh.N));
    IncreaseFacetCounter(facet, currentParticle.flightTime +
                                facet->colDist / 100.0 / currentParticle.velocity, 1, 0, 0,
                         2.0 / (currentParticle.velocity * directionFactor),
                         2.0 * (model.wp.useMaxwellDistribution ? 1.0 : 1.1781) * currentParticle.velocity *
                         directionFactor, currentParticle);

    facet->isHit = true;
    if (/*facet->texture &&*/ facet->sh.countTrans) {
        RecordHitOnTexture(facet, currentParticle.flightTime +
                                  facet->colDist / 100.0 / currentParticle.velocity,
                           true, 2.0, 2.0, currentParticle);
    }
    if (/*facet->direction &&*/ facet->sh.countDirection) {
        RecordDirectionVector(facet, currentParticle.flightTime +
                                     facet->colDist / 100.0 / currentParticle.velocity, currentParticle);
    }
    LogHit(facet, currentParticle);
    ProfileFacet(facet, currentParticle.flightTime + facet->colDist / 100.0 / currentParticle.velocity,
                 true, 2.0, 2.0, currentParticle);
    if (facet->sh.anglemapParams.record) RecordAngleMap(facet, currentParticle);
}
