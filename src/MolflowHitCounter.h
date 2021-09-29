//
// Created by Pascal Baehr on 03.08.20.
//

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
