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

#ifndef MOLFLOW_PROJ_MOLFLOWHITCOUNTER_H
#define MOLFLOW_PROJ_MOLFLOWHITCOUNTER_H

#include <vector>
#include <Buffer_shared.h> // FacetHitBuffer

/*struct FacetState {
    std::vector<FacetHitBuffer> facetHits;
    std::vector<TextureCell>     texture;            // Texture hit recording (taking area, temperature, mass into account), 1+nbMoments
    std::vector<ProfileSlice> profile;         // Distribution and hit recording
    std::vector<DirectionCell>     direction;       // Direction field recording (average), 1+nbMoments
    FacetHistogramBuffer tmpHistograms; //1+nbMoment
};*/

/*struct MolflowHit {
    GlobalHitBuffer globalHits;
    // td ---> 1 + nbMoments
    std::vector<FacetState> facetHits;
    // ----
    //    Anglemap angleMap; // 1 per facet


};*/



/*
 size_t hitSize = 0;
hitSize += sizeof(GlobalHitBuffer) + (1 + nbMoments) * model->wp.globalHistogramParams.GetDataSize();
for (int i = 0; i < model->sh.nbFacet; i++) {
    hitSize += loader.loadFacets[i].GetHitsSize(nbMoments);
}
 */

#endif //MOLFLOW_PROJ_MOLFLOWHITCOUNTER_H
