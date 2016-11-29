/*
  File:        GeometryRender.cpp
  Description: Geometry class (OpenGL rendering/selection stuff)
  Program:     MolFlow
  Author:      R. KERSEVAN / J-L PONS / M ADY
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  */

#include "Geometry.h"
#include "Worker.h"
#include "Utils.h"
#include <malloc.h>
#include <string.h>
#include <math.h>
#include "GLApp/GLMatrix.h"
#ifdef MOLFLOW
#include "MolFlow.h"
#endif

#ifdef SYNRAD
#include "SynRad.h"
#endif
#include "GLApp/GLWindowManager.h"
#include "GLApp/GLMessageBox.h"



#ifdef MOLFLOW
extern MolFlow *mApp;
#endif

#ifdef SYNRAD
extern SynRad*mApp;
#endif

void Geometry::BuildFacetTextures(BYTE *hits) {

	SHGHITS *shGHit = (SHGHITS *)hits;

	Worker *w = &(mApp->worker);
	double dCoef = 1E4; //1E4: conversion between m2 and cm2
	double dCoef_custom[] = { 1.0, 1.0, 1.0 }; //Three coefficients for pressure, imp.rate, density

	int nbMoments = (int)mApp->worker.moments.size();
	size_t facetHitsSize = (1 + nbMoments) * sizeof(SHHITS);
	switch (shGHit->mode) {

	case MC_MODE:

		dCoef_custom[0] = dCoef*mApp->worker.gasMass / 1000 / 6E23*0.0100; //pressure
		dCoef_custom[1] = dCoef;
		dCoef_custom[2] = dCoef;

		for (int i = 0; i < 3; i++) {
			//texture limits already corrected by timeFactor in UpdateMCHits()
			texture_limits[i].autoscale.min.moments_only = shGHit->texture_limits[i].min.moments_only*dCoef_custom[i];
			texture_limits[i].autoscale.max.moments_only = shGHit->texture_limits[i].max.moments_only*dCoef_custom[i];
			texture_limits[i].autoscale.min.all = shGHit->texture_limits[i].min.all*dCoef_custom[i];
			texture_limits[i].autoscale.max.all = shGHit->texture_limits[i].max.all*dCoef_custom[i];
		}
		break;
	case AC_MODE:
		texture_limits[0].autoscale.min.all = shGHit->texture_limits[0].min.all;
		texture_limits[0].autoscale.max.all = shGHit->texture_limits[0].max.all;
		break;

	}

	double iDesorbed = 0.0;
	if (shGHit->total.hit.nbDesorbed)
		iDesorbed = 1.0 / (double)shGHit->total.hit.nbDesorbed;

	GLProgress *prg = new GLProgress("Building texture", "Frame update");
	prg->SetBounds(5, 28, 300, 90);
	prg->SetVisible(TRUE);
	for (int i = 0; i < sh.nbFacet; i++) {
		prg->SetProgress((double)i / (double)sh.nbFacet);
		Facet *f = facets[i];
		GLint max_t;
		glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_t);
		if (f->sh.texHeight > max_t || f->sh.texWidth > max_t) {
			if (!f->textureError) {
				char tmp[1024];
				sprintf(tmp, "Facet #%d has a texture of %dx%d cells.\n"
					"Your video card only supports texture dimensions (width or height) up to %d cells.\n"
					"Texture rendering has been disabled on this facet, but you can still read texture values\n"
					"using the Texture Plotter window. Consider using a smaller mesh resolution, or split the facet\n"
					"into smaller parts. (Use Facet/Explode... command)", i + 1, f->sh.texHeight, f->sh.texWidth, max_t);
				GLMessageBox::Display(tmp, "OpenGL Error", GLDLG_OK, GLDLG_ICONWARNING);
			}
			f->textureError = TRUE;
			return;
		}
		else {

			f->textureError = FALSE;
		}

		int profSize = (f->sh.isProfile) ? (PROFILE_SIZE * sizeof(APROFILE)) : 0;


		int nbElem = f->sh.texWidth*f->sh.texHeight;
		int tSize = nbElem * sizeof(AHIT);
		int dSize = nbElem * sizeof(VHIT);

		if (f->sh.isTextured) {

			// Retrieve texture from shared memory (every seconds)
			AHIT *hits_local = (AHIT *)((BYTE *)shGHit + (f->sh.hitOffset + facetHitsSize + profSize*(1 + nbMoments) + tSize*mApp->worker.displayedMoment));

			double min, max;
			if (!texAutoScale) { //manual values
				min = texture_limits[textureMode].manual.min.all;
				max = texture_limits[textureMode].manual.max.all;

			}
			else { //autoscale
				min = texAutoScaleIncludeConstantFlow ?
					texture_limits[textureMode].autoscale.min.all
					: texture_limits[textureMode].autoscale.min.moments_only;
				max = texAutoScaleIncludeConstantFlow ?
					texture_limits[textureMode].autoscale.max.all
					: texture_limits[textureMode].autoscale.max.moments_only;

			}
			f->BuildTexture(hits_local, textureMode, min, max, texColormap,
				dCoef_custom[0] * mApp->worker.GetMoleculesPerTP(), dCoef_custom[1] * mApp->worker.GetMoleculesPerTP(), dCoef_custom[2] * mApp->worker.GetMoleculesPerTP(), texLogScale, mApp->worker.displayedMoment);
		}
		if (f->sh.countDirection && f->dirCache) {
			VHIT *dirs = (VHIT *)((BYTE *)shGHit + (f->sh.hitOffset + facetHitsSize + profSize*(1 + nbMoments) + tSize*(1 + nbMoments) + dSize*mApp->worker.displayedMoment));
			for (int j = 0; j < nbElem; j++) {
				f->dirCache[j].dir.x = dirs[j].dir.x;
				f->dirCache[j].dir.y = dirs[j].dir.y;
				f->dirCache[j].dir.z = dirs[j].dir.z;
				//f->dirCache[j].sumSpeed = dirs[j].sumSpeed;
				f->dirCache[j].count = dirs[j].count;
			}
		}

	}

	prg->SetVisible(FALSE);
	SAFE_DELETE(prg);
}









