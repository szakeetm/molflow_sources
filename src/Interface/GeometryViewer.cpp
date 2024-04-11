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
#include "GLApp/GLApp.h"
#include "GLApp/GLWindowManager.h"
#include "GeometryViewer.h"
#include "GLApp/GLToolkit.h"
#include "GLApp/GLMatrix.h"
#include "Helper/MathTools.h"
#include <math.h>

#include "GLApp/GLButton.h"
#include "GLApp/GLCombo.h"
#include "GLApp/GLLabel.h"
#include "GLApp/GLGradient.h"
//#include <malloc.h>

#if defined(MOLFLOW)
#include "MolFlow.h"
extern MolFlow *mApp;
#endif

#if defined(SYNRAD)
#include "SynRad.h"
extern SynRad *mApp;
#endif

/**
* \brief To set placement of toolbar at the bottom and it's size
* \param x x position of window
* \param y y position of window
* \param w width of window
* \param h height of window
*/
void GeometryViewer::SetBounds(int x, int y, int width, int height) {
	if (this->posX != x || this->posY != y || this->width != width || this->height != height) {
		GLComponent::SetBounds(x, y, width, height);

		frontBtn->SetBounds(posX + 5, posY + height - 22, 33, 19);
		topBtn->SetBounds(posX + 40, posY + height - 22, 33, 19);
		sideBtn->SetBounds(posX + 75, posY + height - 22, 33, 19);
		projCombo->SetBounds(posX + 110, posY + height - 22, 50, 19);
		coordLab->SetBounds(posX + 162, posY + height - 20, 100, 18);
		facetSearchState->SetBounds(posX + 10, posY + 10, 90, 19);
		timeLabel->SetBounds(posX + width - 200, posY + height - 100, 100, 50);

		autoBtn->SetBounds(posX + width - 122, posY + height - 22, 19, 19);
		selBtn->SetBounds(posX + width - 102, posY + height - 22, 19, 19);

		selVxBtn->SetBounds(posX + width - 82, posY + height - 22, 19, 19);
		zoomBtn->SetBounds(posX + width - 62, posY + height - 22, 19, 19);
		handBtn->SetBounds(posX + width - 42, posY + height - 22, 19, 19);
		sysBtn->SetBounds(posX + width - 22, posY + height - 22, 19, 19);

		/*
		tabLabel->SetBounds(posX + 10, posY + height - 87, 0, 19);
		capsLockLabel->SetBounds(posX + 10, posY + height - 67, 0, 19);
		hideLotlabel->SetBounds(posX + 10, posY + height - 47, 0, 19);
		*/ //Will set these dynamically

		toolBack->SetBounds(posX + 1, posY + height - DOWN_MARGIN, width - 2, DOWN_MARGIN - 1);

		//debugLabel->SetBounds(posX + 10, posY + 35, 90, 19); //debug

		SetCurrentView(view);
	}
}

/**
* \brief Draws Lines and Hits into 3D area via OpenGL.
*/
void GeometryViewer::DrawLinesAndHits() {

    auto& hitCache = mApp->worker.globalStatCache;
	// Lines
	if (showLine && hitCache.hitCacheSize) {

		size_t nbLines = std::min(dispNumHits, hitCache.hitCacheSize);

		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_CULL_FACE);
		
		size_t count = 0;
		while (count < nbLines && hitCache.hitCache[count].type != 0) {

			//Regular (green) line color
			if (mApp->whiteBg) { //whitebg
				glColor3f(0.2f, 0.7f, 0.2f);
			}
			else {
				glColor3f(0.5f, 1.0f, 0.5f);
			}

			if (mApp->antiAliasing) {
				glEnable(GL_BLEND);
				glEnable(GL_LINE_SMOOTH);
			}

			glBegin(GL_LINE_STRIP);
			while (count < nbLines && hitCache.hitCache[count].type != HIT_ABS && hitCache.hitCache[count].type != HIT_VOLUME_DECAY) { //While hits are consecutive

				//change color in case of teleport
				if (hitCache.hitCache[count].type == HIT_TELEPORTSOURCE) {
					glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z); //Draw regular line until TP source
					glEnd(); //Finish regular line
					if (showTP && (count+1)< nbLines && hitCache.hitCache[count+1].type == HIT_TELEPORTDEST) {
						//Switch to orange dashed line
						if (!mApp->whiteBg) {
							glColor3f(1.0f, 0.7f, 0.2f);
						}
						else {
							glColor3f(1.0f, 0.0f, 1.0f);
						}
						glPushAttrib(GL_ENABLE_BIT);

						glLineStipple(1, 0x0101);
						glEnable(GL_LINE_STIPPLE);
						glBegin(GL_LINE_STRIP);
						glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z); //source point
						count++;
						glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z);  //teleport dest.
						glEnd();
						glPopAttrib();

						//Switch back to normal (green) color
						if (mApp->whiteBg) { //whitebg
							glColor3f(0.2f, 0.7f, 0.2f);
						}
						else {
							glColor3f(0.5f, 1.0f, 0.5f);
						}
						glBegin(GL_LINE_STRIP);
					}
					else {

						//glVertex3d(hitCache[count].pos.x , hitCache[count].pos.y , hitCache[count].pos.z); //source point
						count++;
						glBegin(GL_LINE_STRIP);
						glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z);  //teleport dest.

					}
				} //end treating teleport 

				if (hitCache.hitCache[count].type == HIT_LAST) { //pen up at cache refresh border
									glEnd();
									count++;
									glBegin(GL_LINE_STRIP);
				}
				else {
					glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z);
					count++;
				}
			}
			//Treat absorption
			if (count < nbLines && hitCache.hitCache[count].type != 0) {
				glVertex3d(hitCache.hitCache[count].pos.x, hitCache.hitCache[count].pos.y, hitCache.hitCache[count].pos.z);
				count++;
			}
			glEnd();
			if (mApp->antiAliasing) {
				glDisable(GL_LINE_SMOOTH);
				glDisable(GL_BLEND);
			}
		}

	}

	// Hits
	if (showHit) {

		size_t nbHits = std::min(dispNumHits, hitCache.hitCacheSize);

		glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
		glDisable(GL_BLEND);
		glDisable(GL_CULL_FACE);

		float pointSize = (bigDots) ? 3.0f : 2.0f;
		glPointSize(pointSize);

		// Define arrays for positions and colors
		std::vector<GLfloat> vertices; vertices.reserve(3*nbHits);
		std::vector<GLfloat> colors; colors.reserve(3 * nbHits);

		for (size_t i = 0; i < nbHits; i++) {
			switch (hitCache.hitCache[i].type) {
			case HIT_DES:
				colors.push_back(0.3f);
				colors.push_back(0.3f);
				colors.push_back(1.0f);
				break;
			case HIT_REF:
				if (mApp->whiteBg) { //whitebg
					colors.push_back(0.2f);
					colors.push_back(0.2f);
					colors.push_back(0.2f);
				}
				else {
					colors.push_back(0.0f);
					colors.push_back(1.0f);
					colors.push_back(0.0f);
				}
				break;
			case HIT_ABS:
				colors.push_back(1.0f);
				colors.push_back(0.0f);
				colors.push_back(0.0f);
				break;
			case HIT_MOVING:
				colors.push_back(1.0f);
				colors.push_back(0.0f);
				colors.push_back(1.0f);
				break;
			case HIT_SCATTER:
				colors.push_back(1.0f);
				colors.push_back(0.3f);
				colors.push_back(1.0f);
				break;
			case HIT_VOLUME_DECAY:
			case HIT_TRANS:
				colors.push_back(0.5f);
				colors.push_back(1.0f);
				colors.push_back(1.0f);
				break;
			case HIT_TELEPORTSOURCE:
			case HIT_TELEPORTDEST:
				if (!mApp->whiteBg) {
					colors.push_back(1.0f);
					colors.push_back(0.7f);
					colors.push_back(0.2f);
				}
				else {
					colors.push_back(1.0f);
					colors.push_back(0.0f);
					colors.push_back(1.0f);
				}
				break;
			}
			if (hitCache.hitCache[i].type!=HIT_LAST) {
				vertices.push_back(hitCache.hitCache[i].pos.x);
				vertices.push_back(hitCache.hitCache[i].pos.y);
				vertices.push_back(hitCache.hitCache[i].pos.z);
			}
		}

		// Enable vertex arrays
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);

		// Set data pointers
		glVertexPointer(3, GL_FLOAT, 0, vertices.data());
		glColorPointer(3, GL_FLOAT, 0, colors.data());

		// Draw the points
		glDrawArrays(GL_POINTS, 0, nbHits);

		// Disable vertex arrays
		glDisableClientState(GL_COLOR_ARRAY);
		glDisableClientState(GL_VERTEX_ARRAY);





	}
}