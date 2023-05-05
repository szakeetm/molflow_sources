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
#include "MolflowGeometry.h"
#include "Worker.h"
#include "Facet_shared.h"
//#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#if defined(MOLFLOW)
#include "MolFlow.h"
#endif

#if defined(SYNRAD)
#include "SynRad.h"
#endif
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"

#if defined(MOLFLOW)
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
extern SynRad*mApp;
#endif

/**
* \brief Processes events like button clicks for the advanced facet parameters panel.
* \param results contains all simulation results/states (like hits)
* \param renderRegularTexture bool value
* \param renderDirectionTexture bool value
* \param sMode which simulation mode was used (monte carlo / angular coefficient)
*/
void MolflowGeometry::BuildFacetTextures(GlobalSimuState &globState, bool renderRegularTexture,
                                         bool renderDirectionTexture) {
	
	Worker *w = &(mApp->worker);

	int nbMoments = (int)mApp->worker.moments.size();
	size_t facetHitsSize = (1 + nbMoments) * sizeof(FacetHitBuffer);

	auto prg = GLProgress_GUI("Building texture", "Frame update");
	prg.SetBounds(5, 28, 300, 90);
	int startTime = SDL_GetTicks();

	double dCoef_custom[] = { 1.0, 1.0, 1.0 }; //Three coefficients for pressure, imp.rate, density
											   //Autoscaling limits come from the subprocess corrected by "time factor", which makes constant flow and moment values comparable
											   //Time correction factor in subprocess: MoleculesPerTP * nbDesorbed
	double timeCorrection = 1.0;
	double min, max;

	if (renderRegularTexture) {
        {
		    const GlobalHitBuffer& globHit = globState.globalHits;
			dCoef_custom[0] = 1E4 / (double)globHit.globalHits.nbDesorbed * mApp->worker.model->wp.gasMass / 1000 / 6E23*0.0100; //multiplied by timecorr*sum_v_ort_per_area: pressure
			dCoef_custom[1] = 1E4 / (double)globHit.globalHits.nbDesorbed;
			dCoef_custom[2] = 1E4 / (double)globHit.globalHits.nbDesorbed;
			timeCorrection = (mApp->worker.displayedMoment == 0) ? mApp->worker.model->wp.finalOutgassingRate : mApp->worker.model->wp.totalDesorbedMolecules / mApp->worker.moments[mApp->worker.displayedMoment - 1].second;
		}

		if (!texAutoScale) { //manual values
			min = texture_limits[textureMode].manual.min.steady_state;
			max = texture_limits[textureMode].manual.max.steady_state;
		}
		else { //autoscale
            if(texAutoScaleIncludeConstantFlow == 0){
                min = texture_limits[textureMode].autoscale.min.moments_only;
                max = texture_limits[textureMode].autoscale.max.moments_only;
            } else if(texAutoScaleIncludeConstantFlow == 1){
                min = std::min(texture_limits[textureMode].autoscale.min.steady_state, texture_limits[textureMode].autoscale.min.moments_only);
                max = std::max(texture_limits[textureMode].autoscale.max.steady_state, texture_limits[textureMode].autoscale.max.moments_only);
            } else { // == 2
                min = texture_limits[textureMode].autoscale.min.steady_state;
                max = texture_limits[textureMode].autoscale.max.steady_state;
            }
		}
	}

	for (size_t i = 0; i < sh.nbFacet; i++) {
		int time = SDL_GetTicks();
		if (!prg.IsVisible() && ((time - startTime) > 500)) {
			prg.SetVisible(true);
		}
		prg.SetProgress((double)i / (double)sh.nbFacet);
		InterfaceFacet *f = facets[i];

		size_t profSize = (f->sh.isProfile) ? (PROFILE_SIZE * sizeof(ProfileSlice)) : 0;
		size_t nbElem = f->sh.texWidth*f->sh.texHeight;
		size_t tSize = nbElem * sizeof(TextureCell);

		if (renderRegularTexture && f->sh.isTextured) {

			GLint max_t;
			glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_t);
			if (f->sh.texHeight > max_t || f->sh.texWidth > max_t) {
				if (!f->textureError) {
					char tmp[1024];
					sprintf(tmp, "Facet #%zd has a texture of %zdx%zd cells.\n"
						"Your video card only supports texture dimensions (width or height) up to %d cells.\n"
						"Texture rendering has been disabled on this facet, but you can still read texture values\n"
						"using the Texture Plotter window. Consider using a smaller mesh resolution, or split the facet\n"
						"into smaller parts. (Use Facet/Explode... command)", i + 1, f->sh.texHeight, f->sh.texWidth, max_t);
					GLMessageBox::Display(tmp, "OpenGL Error", GLDLG_OK, GLDLG_ICONWARNING);
				}
				f->textureError = true;
				return;
			}
			else {
				f->textureError = false;
			}

			// Retrieve texture from shared memory (every seconds)
			//TextureCell *hits_local = (TextureCell *)((BYTE *)shGHit + (f->sh.hitOffset + facetHitsSize + profSize*(1 + nbMoments) + tSize*mApp->worker.displayedMoment));
			f->BuildTexture(globState.facetStates[i].momentResults[mApp->worker.displayedMoment].texture, textureMode, min, max, texColormap,
				dCoef_custom[0] * timeCorrection, dCoef_custom[1] * timeCorrection, dCoef_custom[2] * timeCorrection, texLogScale, mApp->worker.displayedMoment);
		}

		if (renderDirectionTexture && f->sh.countDirection && f->dirCache) {
			
			size_t dSize = nbElem * sizeof(DirectionCell);

			/*
			double iDesorbed = 0.0;
			if (shGHit->globalHits.hit.nbDesorbed)
			iDesorbed = 1.0 / (double)shGHit->globalHits.hit.nbDesorbed;
			*/


			const std::vector<DirectionCell>& dirs = globState.facetStates[i].momentResults[mApp->worker.displayedMoment].direction;
			for (size_t j = 0; j < nbElem; j++) {
				double denominator = (dirs[j].count > 0) ? 1.0 / dirs[j].count : 1.0;
				f->dirCache[j].dir = dirs[j].dir * denominator;
				f->dirCache[j].count = dirs[j].count;
			}
		}
	}
}

