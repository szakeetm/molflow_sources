
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
void MolflowGeometry::BuildFacetTextures(const std::shared_ptr<GlobalSimuState> globalState, bool renderRegularTexture,
                                         bool renderDirectionTexture) {
	
	Worker *w = &(mApp->worker);

	int nbMoments = (int)mApp->worker.interfaceMomentCache.size();
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
		    const GlobalHitBuffer& globHit = globalState->globalStats;
			dCoef_custom[0] = 1E4 / (double)globHit.globalHits.nbDesorbed * mApp->worker.model->sp.gasMass / 1000 / 6E23*0.0100; //multiplied by timecorr*sum_v_ort_per_area: pressure
			dCoef_custom[1] = 1E4 / (double)globHit.globalHits.nbDesorbed;
			dCoef_custom[2] = 1E4 / (double)globHit.globalHits.nbDesorbed;
			timeCorrection = (mApp->worker.displayedMoment == 0) ? mApp->worker.model->sp.finalOutgassingRate : mApp->worker.model->sp.totalDesorbedMolecules / mApp->worker.interfaceMomentCache[mApp->worker.displayedMoment - 1].window;
		}

		if (!texAutoScale) { //manual values: Always "steady state" variable used
			min = texture_limits[textureMode].manual.min.steady_state;
			max = texture_limits[textureMode].manual.max.steady_state;
		}
		else { //autoscale
			std::tie(min, max) = GetTextureAutoscaleMinMax();
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
			f->BuildTexture(globalState->facetStates[i].momentResults[mApp->worker.displayedMoment].texture, textureMode, min, max, texColormap,
				dCoef_custom[0] * timeCorrection, dCoef_custom[1] * timeCorrection, dCoef_custom[2] * timeCorrection, texLogScale, mApp->worker.displayedMoment);
		}

		if (renderDirectionTexture && f->sh.countDirection) {
			f->dirCache = globalState->facetStates[i].momentResults[mApp->worker.displayedMoment].direction;
			for (auto& dirCell : f->dirCache) {
				if (dirCell.count > 0) {
					dirCell.dir = dirCell.dir * (1.0 / dirCell.count);
				}
			}
		}
	}
}

