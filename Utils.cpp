/*
  File:        Utils.cpp
  Description: Various util functions
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

#include "Types.h"
#include "Utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

// -------------------------------------------------------
// Linear algebra
// -------------------------------------------------------

void   Cross(VERTEX3D *result, VERTEX3D *v1, VERTEX3D *v2) {
	result->x = (v1->y)*(v2->z) - (v1->z)*(v2->y);
	result->y = (v1->z)*(v2->x) - (v1->x)*(v2->z);
	result->z = (v1->x)*(v2->y) - (v1->y)*(v2->x);
}

double Dot(VERTEX3D *v1, VERTEX3D *v2) {
	return (v1->x)*(v2->x) + (v1->y)*(v2->y) + (v1->z)*(v2->z);
}

double Norme(VERTEX3D *v) {
	return sqrt(Dot(v, v));
}

int Remainder(int param, int bound) {
	return (bound > 0) ? ((bound + (param) % bound) % bound) : 0;
}

void Normalize(VERTEX3D *v) {
	double length = Norme(v);
	if (length>0.0)
		ScalarMult(v, 1.0 / length);
}

double Norme(VERTEX2D *v1, VERTEX2D *v2) {
	return sqrt((v1->u - v2->u)*(v1->u - v2->u) + (v1->v - v2->v)*(v1->v - v2->v));
}

void ScalarMult(VERTEX3D *result, double r) {
	result->x *= r;
	result->y *= r;
	result->z *= r;
}

void Sub(VERTEX3D *result, VERTEX3D *v1, VERTEX3D *v2) {
	result->x = (v1->x) - (v2->x);
	result->y = (v1->y) - (v2->y);
	result->z = (v1->z) - (v2->z);
}

void Add(VERTEX3D *result, VERTEX3D *v1, VERTEX3D *v2) {
	result->x = (v1->x) + (v2->x);
	result->y = (v1->y) + (v2->y);
	result->z = (v1->z) + (v2->z);
}

int VertexEqual(VERTEX2D *p1, VERTEX2D *p2) {
	return IS_ZERO(p1->u - p2->u) && IS_ZERO(p1->v - p2->v);
}

void Sub2(VERTEX2D *result, VERTEX2D *v1, VERTEX2D *v2) {
	result->u = (v1->u) - (v2->u);
	result->v = (v1->v) - (v2->v);
}

void   ProjectVertex(VERTEX3D *v, VERTEX2D *projected, VERTEX3D U, VERTEX3D V, VERTEX3D origin){
	//Project v on a plane defined by U,V and return the coordinates in base U,V
	VERTEX3D AP, APp, projectionPoint, mult;

	Sub(&AP, v, &origin);
	double t = Dot(&U, &AP) / Dot(&U, &U);
	mult = U;
	ScalarMult(&mult, t);
	Add(&projectionPoint, &origin, &mult);
	Sub(&APp, &projectionPoint, &origin);
	projected->u = Dot(&APp, &U); //u coordinate

	Sub(&AP, v, &origin);
	t = Dot(&V, &AP) / Dot(&V, &V);
	mult = V;
	ScalarMult(&mult, t);
	Add(&projectionPoint, &origin, &mult);
	Sub(&APp, &projectionPoint, &origin);
	projected->v = Dot(&APp, &V); //v coordinate
}

void Mirror(VERTEX3D *P, VERTEX3D P0, VERTEX3D N) {
	VERTEX3D P0_P, P_copy, N_copy;
	P_copy = *P; N_copy = N;
	Normalize(&N_copy);
	Sub(&P0_P, &P_copy, &P0);
	double n_dot_p = Dot(&N_copy, &P0_P);
	ScalarMult(&N_copy, 2.0*n_dot_p);
	Sub(P, &P_copy, &N_copy);
}

void Rotate(VERTEX3D *P, VERTEX3D AXIS_P0, VERTEX3D AXIS_DIR, double theta) {
	theta = theta / 180 * PI; //degree->radians
	Normalize(&AXIS_DIR);
	double x, y, z, a, b, c, u, v, w, costh, sinth, precalc1;
	x = P->x;
	y = P->y;
	z = P->z;
	a = AXIS_P0.x;
	b = AXIS_P0.y;
	c = AXIS_P0.z;
	u = AXIS_DIR.x;
	v = AXIS_DIR.y;
	w = AXIS_DIR.z;
	costh = cos(theta);
	sinth = sin(theta);
	precalc1 = -u*x - v*y - w*z;
	P->x = (a*(v*v + w*w) - u*(b*v + c*w + precalc1))*(1 - costh) + x*costh + (-c*v + b*w - w*y + v*z)*sinth;
	P->y = (b*(u*u + w*w) - v*(a*u + c*w + precalc1))*(1 - costh) + y*costh + (c*u - a*w + w*x - u*z)*sinth;
	P->z = (c*(u*u + v*v) - w*(a*u + b*v + precalc1))*(1 - costh) + z*costh + (-b*u + a*v - v*x + u*y)*sinth;
}

int IsEqual(const double &a, const double &b, double tolerance) {
	return abs(a - b) < tolerance;
}

// ---------------------------------------------------------------
// Polygon geometry stuff
// ---------------------------------------------------------------

// -----------------------------------------------------------

int IsConvex(POLYGON *p, int idx) {

	// Check if p.pts[idx] is a convex vertex (calculate the sign of the oriented angle)

	int i1 = IDX(idx - 1, p->nbPts);
	int i2 = IDX(idx, p->nbPts);
	int i3 = IDX(idx + 1, p->nbPts);

	double d = DET22(p->pts[i1].u - p->pts[i2].u, p->pts[i3].u - p->pts[i2].u,
		p->pts[i1].v - p->pts[i2].v, p->pts[i3].v - p->pts[i2].v);

	return (d*p->sign) >= 0.0;

}

// -----------------------------------------------------------

int IsInsideTri(VERTEX2D *p, VERTEX2D *p1, VERTEX2D *p2, VERTEX2D *p3) {

	// Check if p is inside the triangle p1 p2 p3
	VERTEX2D pts[3];
	pts[0] = *p1;
	pts[1] = *p2;
	pts[2] = *p3;
	return IsInPoly(p->u, p->v, pts, 3);

}

// -----------------------------------------------------------

int ContainsConcave(POLYGON *p, int i1, int i2, int i3) {

	// Determine if the specified triangle contains or not a concave point
	int _i1 = IDX(i1, p->nbPts);
	int _i2 = IDX(i2, p->nbPts);
	int _i3 = IDX(i3, p->nbPts);

	VERTEX2D *p1 = p->pts + _i1;
	VERTEX2D *p2 = p->pts + _i2;
	VERTEX2D *p3 = p->pts + _i3;

	int found = 0;
	int i = 0;
	while (!found && i < p->nbPts) {
		if (i != _i1 && i != _i2 && i != _i3) {
			VERTEX2D *pt = p->pts + i;
			if (IsInsideTri(pt, p1, p2, p3))
				found = !IsConvex(p, i);
		}
		i++;
	}

	return found;

}

// -----------------------------------------------------------

int EmptyTriangle(POLYGON *p, int i1, int i2, int i3, VERTEX2D *center) {

	// Determine if the specified triangle contains or not an other point of the poly
	int _i1 = IDX(i1, p->nbPts);
	int _i2 = IDX(i2, p->nbPts);
	int _i3 = IDX(i3, p->nbPts);

	VERTEX2D *p1 = p->pts + _i1;
	VERTEX2D *p2 = p->pts + _i2;
	VERTEX2D *p3 = p->pts + _i3;

	int found = 0;
	int i = 0;
	while (!found && i < p->nbPts) {
		if (i != _i1 && i != _i2 && i != _i3) {
			VERTEX2D *pt = p->pts + i;
			found = IsInsideTri(pt, p1, p2, p3);
		}
		i++;
	}

	center->u = (p1->u + p2->u + p3->u) / 3.0;
	center->v = (p1->v + p2->v + p3->v) / 3.0;
	return !found;

}

// -----------------------------------------------------------

int IsInPoly(double u, double v, VERTEX2D *pts, int nbPts) {

	// 2D polygon "is inside" solving
	// Using the "Jordan curve theorem" (we intersect in v direction here)

	int n_updown, n_found, j;
	double x1, x2, y1, y2, a, minx, maxx;

	n_updown = 0;
	n_found = 0;

	for (j = 0; j < nbPts - 1; j++) {

		x1 = pts[j].u;
		y1 = pts[j].v;
		x2 = pts[j + 1].u;
		y2 = pts[j + 1].v;

		if (x2>x1) { minx = x1; maxx = x2; }
		else        { minx = x2; maxx = x1; }

		if (u > minx && u <= maxx) {
			a = (y2 - y1) / (x2 - x1);
			if ((a*(u - x1) + y1) < v) {
				n_updown = n_updown + 1;
			}
			else {
				n_updown = n_updown - 1;
			}
			n_found++;
		}

	}

	// Last point
	x1 = pts[j].u;
	y1 = pts[j].v;
	x2 = pts[0].u;
	y2 = pts[0].v;

	if (x2 > x1) { minx = x1; maxx = x2; }
	else        { minx = x2; maxx = x1; }

	if (u > minx && u <= maxx) {
		a = (y2 - y1) / (x2 - x1);
		if ((a*(u - x1) + y1) < v) {
			n_updown = n_updown + 1;
		}
		else {
			n_updown = n_updown - 1;
		}
		n_found++;
	}

	if (n_updown < 0) n_updown = -n_updown;
	return (((n_found / 2) & 1) ^ ((n_updown / 2) & 1));

}

// -----------------------------------------------------------

int Intersect2D(VERTEX2D *p1, VERTEX2D *p2, VERTEX2D *p3, VERTEX2D *p4, VERTEX2D *I) {

	// Computes the intersection between 2 segments
	// Return 1 when an intersection is found, 0 otherwise
	// Solve P1 + t1P1P2 = P3 + t2P3P4

	VERTEX2D p12;
	p12.u = p2->u - p1->u;
	p12.v = p2->v - p1->v;
	VERTEX2D p34;
	p34.u = p4->u - p3->u;
	p34.v = p4->v - p3->v;

	double det = DET22(-p12.u, p34.u,
		-p12.v, p34.v);
	if (IS_ZERO(det)) return 0;

	double idet = 1.0 / det;

	double dt1 = DET22(p1->u - p3->u, p34.u,
		p1->v - p3->v, p34.v);
	double t1 = dt1*idet;
	if (t1<0.0 || t1>1.0) return 0;

	double dt2 = DET22(-p12.u, p1->u - p3->u,
		-p12.v, p1->v - p3->v);
	double t2 = dt2*idet;
	if (t2<0.0 || t2>1.0) return 0;

	// Check coherence (numerical error)
	VERTEX2D I1, I2;
	I1.u = p1->u + t1*p12.u;
	I1.v = p1->v + t1*p12.v;
	I2.u = p3->u + t2*p34.u;
	I2.v = p3->v + t2*p34.v;
	double r = sqrt((I1.u - I2.u)*(I1.u - I2.u) + (I1.v - I2.v)*(I1.v - I2.v));
	if (r > 1e-6) return 0;

	I->u = I1.u;
	I->v = I1.v;
	return 1;

}

// -----------------------------------------------------------

int IsOnEdge(VERTEX2D *p1, VERTEX2D *p2, VERTEX2D *p) {

	// Returns 1 wether p lies in [P1P2], 0 otherwise

	double t = DOT2(p2->u - p1->u, p2->v - p1->v,
		p->u - p1->u, p->v - p1->v);

	double n1 = sqrt(DOT2(p->u - p1->u, p->v - p1->v,
		p->u - p1->u, p->v - p1->v));

	double n2 = sqrt(DOT2(p2->u - p1->u, p2->v - p1->v,
		p2->u - p1->u, p2->v - p1->v));

	double c = t / (n1*n2);

	// Check that p1,p2 and p are aligned
	if (c != 1.0)  return 0;

	// Check wether p is between p1 and p2
	t = n1 / n2;
	return (t >= 0.0 && t <= 1.0);

}

// -----------------------------------------------------------

int GetNode(POLYGRAPH *g, VERTEX2D *p) {

	// Search a node in the polygraph and returns its id (-1 not found)

	int found = 0;
	int i = 0;
	int nb = g->nbNode;

	while (i < nb && !found) {
		found = VertexEqual(&(g->nodes[i].p), p);
		if (!found) i++;
	}

	if (found)
		return i;
	else
		return -1;

}

// -----------------------------------------------------------

int SearchFirst(POLYGRAPH *g, POLYVERTEX **s) {

	// Search a not yet processed starting point
	//(Used by IntersectPoly)

	int found = 0;
	int i = 0;

	while (i < g->nbNode && !found) {
		found = (g->nodes[i].mark == 0) && (g->nodes[i].isStart);
		if (!found) i++;
	}

	if (found)
		*s = &(g->nodes[i]);

	return found;

}

// -----------------------------------------------------------

int AddNode(POLYGRAPH *g, VERTEX2D *p) {

	// Add a node to the polygraph and returns its id

	int i = GetNode(g, p);

	if (i < 0) {
		// New node
		int nb = g->nbNode;
		g->nodes[nb].p = *p;
		g->nodes[nb].VI[0] = -1;
		g->nodes[nb].VI[1] = -1;
		g->nodes[nb].VO[0] = -1;
		g->nodes[nb].VO[1] = -1;
		g->nbNode++;
		return nb;
	}

	return i;

}

// -----------------------------------------------------------

int AddArc(POLYGRAPH *g, int i1, int i2, int source) {

	// Add an arc to the polygraph and returns its id

	int found = 0;
	int i = 0;
	int nb = g->nbArc;

	while (i < nb && !found) {
		found = (g->arcs[i].i1 == i1 && g->arcs[i].i2 == i2);
		if (!found) i++;
	}

	if (!found) {
		// New arc
		g->arcs[nb].i1 = i1;
		g->arcs[nb].i2 = i2;
		g->arcs[i].s = source;
		g->nbArc++;
	}

	return i;

}

// -----------------------------------------------------------

void CutArc(POLYGRAPH *g, int idx, int ni) {

	// Cut the arc idx by inserting ni
	if (g->arcs[idx].i1 != ni && g->arcs[idx].i2 != ni) {
		int tmp = g->arcs[idx].i2;
		g->arcs[idx].i2 = ni;
		AddArc(g, ni, tmp, 1);
	}

}

// -----------------------------------------------------------

void InsertEdge(POLYGRAPH *g, VERTEX2D *p1, VERTEX2D *p2, int a0) {

	// Insert a polygon edge in the polygraph, a0 = first arc to be checked

	if (VertexEqual(p1, p2))
		// Does not add null arc
		return;

	// Insert nodes
	int n1 = AddNode(g, p1);
	int n2 = AddNode(g, p2);

	// Check intersection of the new arc with the arcs of the first polygon.
	int itFound = 0;
	int i = a0;
	while (i < g->nbArc && !itFound) {

		if (g->arcs[i].s == 1) {

			VERTEX2D *e1 = &(g->nodes[g->arcs[i].i1].p);
			VERTEX2D *e2 = &(g->nodes[g->arcs[i].i2].p);
			VERTEX2D I;

			if (Intersect2D(p1, p2, e1, e2, &I)) {
				int ni = AddNode(g, &I);
				InsertEdge(g, p1, &I, i + 1);
				InsertEdge(g, &I, p2, i + 1);
				CutArc(g, i, ni);
				itFound = 1;
			}

		}
		i++;

	}

	if (!itFound) AddArc(g, n1, n2, 2);

}

// -----------------------------------------------------------

void ClearGraph(POLYGRAPH *g) {
	free(g->arcs);
	free(g->nodes);
}

// -----------------------------------------------------------

POLYGON *CopyPoly(POLYGON *p) {

	POLYGON *polys = (POLYGON *)malloc(sizeof(POLYGON));
	polys[0].nbPts = p->nbPts;
	polys[0].pts = (VERTEX2D *)malloc(p->nbPts * sizeof(VERTEX2D));
	memcpy(polys[0].pts, p->pts, p->nbPts * sizeof(VERTEX2D));
	polys[0].sign = p->sign;

	return polys;

}

// -----------------------------------------------------------

double GetOrientedAngle(VERTEX2D *v1, VERTEX2D *v2) {

	// Return oriented angle [0,2PI] (clockwise)
	double n = sqrt((v1->u*v1->u + v1->v*v1->v) * (v2->u*v2->u + v2->v*v2->v));
	double cs = DOT2(v1->u, v1->v, v2->u, v2->v) / n;
	if (cs<-1.0) cs = -1.0;
	if (cs>1.0)  cs = 1.0;
	double a = acos(cs);
	double s = DET22(v1->u, v2->u, v1->v, v2->v);
	if (s < 0.0) a = 2.0*PI - a;
	return 2.0*PI - a;

}

// -----------------------------------------------------------

void CreateGraph(POLYGRAPH *g, POLYGON *inP1, POLYGON *inP2, int *visible2) {

	// Create the polygraph which represent the 2 intersected polygons
	// with their oriented edges.

	int MAXEDGE = inP1->nbPts * inP2->nbPts + 1;

	g->nodes = (POLYVERTEX *)malloc(MAXEDGE * sizeof(POLYVERTEX));
	memset(g->nodes, 0, MAXEDGE * sizeof(POLYVERTEX));
	g->arcs = (POLYARC *)malloc(MAXEDGE * 4 * sizeof(POLYARC));
	memset(g->arcs, 0, MAXEDGE * sizeof(POLYARC));

	// Fill up the graph with the 1st polygon
	g->nbArc = inP1->nbPts;
	g->nbNode = inP1->nbPts;

	for (int i = 0; i < inP1->nbPts; i++) {
		g->nodes[i].p = inP1->pts[i];
		g->nodes[i].VI[0] = -1;
		g->nodes[i].VI[1] = -1;
		g->nodes[i].VO[0] = -1;
		g->nodes[i].VO[1] = -1;
		g->arcs[i].i1 = i;
		g->arcs[i].i2 = IDX(i + 1, inP1->nbPts);
		g->arcs[i].s = 1;
	}

	// Intersect with 2nd polygon
	for (int i = 0; i < inP2->nbPts; i++)  {
		int i2 = IDX(i + 1, inP2->nbPts);
		if ((!visible2 || (visible2 && visible2[i]))) {
			if (inP2->sign < 0.0) {
				InsertEdge(g, inP2->pts + i2, inP2->pts + i, 0);
			}
			else {
				InsertEdge(g, inP2->pts + i, inP2->pts + i2, 0);
			}
		}
	}

	// Remove tangent edge
	for (int i = 0; i<g->nbArc; i++) {
		if ((g->arcs[i].s>0)) {
			int j = i + 1;
			int found = 0;
			while (j<g->nbArc && !found) {
				if ((g->arcs[j].s>0) &&
					(g->arcs[j].i1 == g->arcs[i].i2) &&
					(g->arcs[j].i2 == g->arcs[i].i1))
				{
					g->arcs[i].s = 0;
					g->arcs[j].s = 0;
				}
				if (!found) j++;
			}
		}
	}

	// Fill up successor in the polyvertex array to speed up search
	// of next vertices

	for (int i = 0; i<g->nbArc; i++) {
		if (g->arcs[i].s>0) {
			int idxO = g->arcs[i].i1;
			int idxI = g->arcs[i].i2;
			if (g->nodes[idxI].nbIn < 2) {
				g->nodes[idxI].VI[g->arcs[i].s - 1] = g->arcs[i].i1;
				g->nodes[idxI].nbIn++;
			}
			if (g->nodes[idxO].nbOut < 2) {
				g->nodes[idxO].VO[g->arcs[i].s - 1] = g->arcs[i].i2;
				g->nodes[idxO].nbOut++;
			}
		}
	}

	// Mark starting points (2 outgoing arcs)

	for (int i = 0; i < g->nbNode; i++) {
		if (g->nodes[i].nbOut >= 2) {
			if (g->nodes[i].nbIn >= 2) {

				// Check tangent point
				VERTEX2D vi1, vi2, vo1, vo2;

				// TO DEBUG!!! Causes frequent crashes
				if (g->nodes[i].VI[0] >= 0 && g->nodes[i].VO[0] >= 0 && g->nodes[i].VI[1] >= 0 && g->nodes[i].VO[1] >= 0) {
					Sub2(&vi1, &(g->nodes[g->nodes[i].VI[0]].p), &(g->nodes[i].p));
					Sub2(&vo1, &(g->nodes[g->nodes[i].VO[0]].p), &(g->nodes[i].p));

					Sub2(&vi2, &(g->nodes[g->nodes[i].VI[1]].p), &(g->nodes[i].p));
					Sub2(&vo2, &(g->nodes[g->nodes[i].VO[1]].p), &(g->nodes[i].p));
				}

				double angI = GetOrientedAngle(&vi1, &vo1);
				double angII = GetOrientedAngle(&vi1, &vi2);
				double angIO = GetOrientedAngle(&vi1, &vo2);

				g->nodes[i].isStart = (angII < angI) || (angIO < angI);

			}
			else {
				g->nodes[i].isStart = 1;
			}
		}
	}

}

// -----------------------------------------------------------

int CheckLoop(POLYGRAPH *g) {

	// The grapth is assumed to not contains starting point
	// 0 => More than one loop
	// 1 => tangent node/egde detected

	POLYVERTEX *s, *s0;
	s = &(g->nodes[0]);
	s0 = s;
	int nbVisited = 0;
	int ok = s->nbOut == 1;

	do {
		if (ok) {
			if (s->VO[0] >= 0)
				s = &(g->nodes[s->VO[0]]);
			else
				s = &(g->nodes[s->VO[1]]);
			ok = (s->nbOut == 1);
			nbVisited++;
		}
	} while ((s0 != s) && (ok));

	// !ok                   => only tangent node
	// nbVisited==g->nbNode  => only tangent edge

	//if( !ok || nbVisited==g->nbNode )
	if (nbVisited == g->nbNode)
		return 1;

	return 0;

}

// -----------------------------------------------------------

void FreePolys(POLYGON **polys, int nbPoly) {

	POLYGON *p = *polys;
	for (int i = 0; i < nbPoly; i++) free(p[i].pts);
	free(p);
	*polys = NULL;

}

// -----------------------------------------------------------

int IntersectPoly(POLYGON *inP1, POLYGON *inP2, int *visible2, POLYGON **result) {

	// Computes the polygon intersection between p1 and p2.
	// Operates on simple polygon only. (Hole connection segments 
	// must be marked non visible in visible2, visible2=NULL => all are visible)
	// Return the number of polygon created (0 on null intersection, -1 on failure)

	POLYGRAPH graph;
	POLYGRAPH *g = &graph;
	int MAXEDGE = inP1->nbPts * inP2->nbPts + 1;

	// Create polygraph
	CreateGraph(g, inP1, inP2, visible2);

	// Search a divergent point
	POLYVERTEX *s, *s0;

	if (!SearchFirst(g, &s)) {

		// Check particular cases

		if ((g->nbNode == inP1->nbPts) && (g->nbNode == inP2->nbPts)) {
			// P1 and P2 are equal
			ClearGraph(g);
			*result = CopyPoly(inP1);
			return 1;
		}

		if (CheckLoop(g)) {
			// Only tangent edge/point found => null intersection
			ClearGraph(g);
			*result = NULL;
			return 0;
		}

		ClearGraph(g);
		int i;
		int insideP1 = 0;
		int insideP2 = 0;

		i = 0;
		while (i < inP1->nbPts && !insideP2) {
			insideP2 = IsInPoly(inP1->pts[i].u, inP1->pts[i].v, inP2->pts, inP2->nbPts);
			i++;
		}
		if (insideP2) {
			// P1 is fully inside P2
			*result = CopyPoly(inP1);
			return 1;
		}

		i = 0;
		while (i < inP2->nbPts && !insideP1) {
			insideP1 = IsInPoly(inP2->pts[i].u, inP2->pts[i].v, inP1->pts, inP1->nbPts);
			i++;
		}
		if (insideP1) {
			// P2 is fully inside P1
			*result = CopyPoly(inP2);
			return 1;
		}

		// Null intersection
		*result = NULL;
		return 0;

	}

	// Compute intersection
	POLYGON *polys = (POLYGON *)malloc(sizeof(POLYGON)* 256);
	int nbPoly = 0;
	int eop;
	VERTEX2D n1, n2;
	double sine;

	do {

		// Starts a new polygon
		polys[nbPoly].pts = (VERTEX2D *)malloc(MAXEDGE*sizeof(VERTEX2D));
		polys[nbPoly].sign = 1.0;
		polys[nbPoly].nbPts = 0;
		nbPoly++;
		eop = 0;
		s0 = s;

		while (!eop && polys[nbPoly - 1].nbPts < MAXEDGE) {

			// Add point to the current polygon
			polys[nbPoly - 1].pts[polys[nbPoly - 1].nbPts] = s->p;
			polys[nbPoly - 1].nbPts++;
			s->mark = 1;

			// Go to next point
			switch (s->nbOut) {

			case 1:

				// Next point
				if (s->VO[0] >= 0) {
					// On a P1 edge
					s = &(g->nodes[s->VO[0]]);
				}
				else {
					// On a P2 edge
					s = &(g->nodes[s->VO[1]]);
				}
				break;

			case 2:

				if (s->VO[0] == -1 || s->VO[1] == -1) {
					//Failure!!! (tangent edge not marked)
					FreePolys(&polys, nbPoly);
					ClearGraph(g);
					*result = NULL;
					return -1;
				}

				// We have to turn left
				n1 = g->nodes[s->VO[0]].p;
				n2 = g->nodes[s->VO[1]].p;

				sine = DET22(n1.u - (s->p.u), n2.u - (s->p.u),
					n1.v - (s->p.v), n2.v - (s->p.v));

				if (sine < 0.0)
					// Go to n1
					s = &(g->nodes[s->VO[0]]);
				else
					// Go to n2
					s = &(g->nodes[s->VO[1]]);

				break;

			default:
				//Failure!!! (not ended polygon)
				FreePolys(&polys, nbPoly);
				ClearGraph(g);
				*result = NULL;
				return -1;

			}

			// Reach start point, end of polygon
			eop = (s0 == s);

		}

		if (!eop) {
			//Failure!!! (inner cycle found)
			FreePolys(&polys, nbPoly);
			ClearGraph(g);
			*result = NULL;
			return -1;
		}

	} while (SearchFirst(&graph, &s));

	ClearGraph(g);
	*result = polys;
	return nbPoly;

}

// -------------------------------------------------------

double GetInterArea(POLYGON *inP1, POLYGON *inP2, int *edgeVisible, float *uC, float *vC, int *nbV, double **lList) {

	int nbPoly;
	double A0;
	POLYGON *polys;
	*nbV = 0;
	*lList = NULL;
	*uC = 0.0;
	*vC = 0.0;

	nbPoly = IntersectPoly(inP1, inP2, edgeVisible, &polys);
	if (nbPoly < 0)
		return 0.0;
	if (nbPoly == 0)
		return 0.0;

	// Count number of pts
	int nbE = 0;
	for (int i = 0; i < nbPoly; i++) nbE += polys[i].nbPts;
	*lList = (double *)malloc(nbE * 2 * sizeof(double));
	*nbV = nbE;

	// Area
	nbE = 0;
	double sum = 0.0;
	for (int i = 0; i < nbPoly; i++) {
		double A = 0.0;
		for (int j = 0; j < polys[i].nbPts; j++) {
			int j1 = IDX(j + 1, polys[i].nbPts);
			A += polys[i].pts[j].u*polys[i].pts[j1].v - polys[i].pts[j1].u*polys[i].pts[j].v;
			(*lList)[nbE++] = polys[i].pts[j].u;
			(*lList)[nbE++] = polys[i].pts[j].v;
		}
		if (i == 0) A0 = fabs(0.5 * A);
		sum += fabs(0.5 * A);
	}

	// Centroid (polygon 0)
	double xC = 0.0;
	double yC = 0.0;
	for (int j = 0; j < polys[0].nbPts; j++) {
		int j1 = IDX(j + 1, polys[0].nbPts);
		double d = polys[0].pts[j].u*polys[0].pts[j1].v - polys[0].pts[j1].u*polys[0].pts[j].v;
		xC += (polys[0].pts[j].u + polys[0].pts[j1].u)*d;
		yC += (polys[0].pts[j].v + polys[0].pts[j1].v)*d;
	}
	*uC = (float)(xC / (6.0 * A0));
	*vC = (float)(yC / (6.0 * A0));

	FreePolys(&polys, nbPoly);

	return sum;

}


// -------------------------------------------------------

double GetInterAreaBF(POLYGON *inP1, double u0, double v0, double u1, double v1, float *uC, float *vC) {

	// Compute area of the intersection between the (u0,v0,u1,v1) rectangle 
	// and the inP1 polygon using Jordan theorem.
	// Slow but sure.

	int step = 50;
	int nbHit = 0;
	double ui = (u1 - u0) / (double)step;
	double vi = (v1 - v0) / (double)step;

	for (int i = 0; i < step; i++) {
		double uc = u0 + ui*((double)i + 0.5);
		for (int j = 0; j < step; j++) {
			double vc = v0 + vi*((double)j + 0.5);
			if (IsInPoly(uc, vc, inP1->pts, inP1->nbPts)) {
				nbHit++;
				*uC = (float)uc;
				*vC = (float)vc;
			}
		}
	}

	return (u1 - u0)*(v1 - v0)*((double)nbHit / (double)(step*step));

}

// -------------------------------------------------------
// Various utils function
// -------------------------------------------------------

// Return a power of 2 which is greater or equal than n
int GetPower2(int n) {

	if ((n & (n - 1)) == 0) {
		// already a power of 2
		return n;
	}
	else {
		// Get the power of 2 above
		int p = 0;
		while (n != 0) { n = n >> 1; p++; }
		return 1 << p;
	}

}

// Format a number of byte in KB,MB,...
char *FormatMemory(size_t size) {
	return FormatMemoryLL((llong)size);
}

char *FormatMemoryLL(llong size) {

	static char ret[256];
	const char *suffixStr[] = { "KB", "MB", "GB", "TB", "PB" };
	double dSize = (double)size;
	int suffix = 0;

	while (dSize >= 1024.0 && suffix < 4) {
		dSize /= 1024.0;
		suffix++;
	}

	if (suffix == 0) {
		sprintf(ret, "%u bytes", (unsigned int)size);
	}
	else {
		if (fabs(dSize - floor(dSize)) < 1e-3)
			sprintf(ret, "%.0f%s", dSize, suffixStr[suffix - 1]);
		else
			sprintf(ret, "%.2f%s", dSize, suffixStr[suffix - 1]);
	}
	return ret;

}

// Return a in [-PI,PI]
double RoundAngle(double a) {

	double r = a;
	while (r<-PI) r += 2.0*PI;
	while (r> PI) r -= 2.0*PI;
	return r;

}

#ifdef NOTUSED

// -------------------------------------------------------
// Assembly optimisation for FPU
// -------------------------------------------------------
int  SolveIASM(double *u ,double *v,double *w,
	VERTEX3D *nuv,VERTEX3D *U,VERTEX3D *V,VERTEX3D *W,VERTEX3D *Z) {

	double det,iDet,detU,detV,detW,tu,tv,tw,dummy;

	__asm {

		mov ebx , dword ptr [nuv]
			mov edx , dword ptr [W]

			//  double det = DET33( U->x , V->x , W->x ,
			//                      U->y , V->y , W->y ,
			//                      U->z , V->z , W->z );

			fld         qword ptr [ebx]
			fmul        qword ptr [edx]
			fld         qword ptr [ebx+8]
			fmul        qword ptr [edx+8]
			faddp       st(1),st
			fld         qword ptr [ebx+10h]
			fmul        qword ptr [edx+10h]
			faddp       st(1),st
			fstp        qword ptr [det]

			//  if(det==0.0) return 0;   // U,V and W are coplanar

			fldz
			fld         qword ptr [det]
			fcomip      st,st(1)
			je          not_found
			fstp        qword ptr [dummy]

			//  double iDet = 1.0 / det;

			fld1
			fdiv        qword ptr [det]
			fstp        qword ptr [iDet]

			//  double detU = DET33( Z->x , V->x , W->x ,
			//                       Z->y , V->y , W->y ,
			//                       Z->z , V->z , W->z );

			mov ebx , dword ptr [U]
			mov ecx , dword ptr [V]
			mov edi , dword ptr [Z]

			fld         qword ptr [ecx+8]
			fmul        qword ptr [edx+10h]
			fld         qword ptr [ecx+10h]
			fmul        qword ptr [edx+8]
			fsubp       st(1),st
			fmul        qword ptr [edi]
			fld         qword ptr [edx+8]
			fmul        qword ptr [edi+10h]
			fld         qword ptr [edx+10h]
			fmul        qword ptr [edi+8]
			fsubp       st(1),st
			fmul        qword ptr [ecx]
			faddp       st(1),st
			fld         qword ptr [edi+8]
			fmul        qword ptr [ecx+10h]
			fld         qword ptr [edi+10h]
			fmul        qword ptr [ecx+8]
			fsubp       st(1),st
			fmul        qword ptr [edx]
			faddp       st(1),st
			fstp        qword ptr [detU]

			//  double tu = detU*iDet;

			fld         qword ptr [detU]
			fmul        qword ptr [iDet]
			fstp        qword ptr [tu]

			//  if( tu<0.0 || tu>1.0 ) return 0; // Bounding rectange exclusion

			fldz
			fld         qword ptr [tu]
			fcomip      st,st(1)
			jb          not_found
			fstp        qword ptr [dummy]
			fld1
			fld         qword ptr [tu]
			fcomip      st,st(1)
			ja          not_found
			fstp        qword ptr [dummy]

			//  double detV = DET33( U->x , Z->x , W->x ,
			//                       U->y , Z->y , W->y ,
			//                       U->z , Z->z , W->z );

			fld         qword ptr [edi+8]
			fmul        qword ptr [edx+10h]
			fld         qword ptr [edi+10h]
			fmul        qword ptr [edx+8]
			fsubp       st(1),st
			fmul        qword ptr [ebx]
			fld         qword ptr [edx+8]
			fmul        qword ptr [ebx+10h] 
			fld         qword ptr [edx+10h] 
			fmul        qword ptr [ebx+8] 
			fsubp       st(1),st 
			fmul        qword ptr [edi]
			faddp       st(1),st
			fld         qword ptr [ebx+8] 
			fmul        qword ptr [edi+10h]
			fld         qword ptr [ebx+10h] 
			fmul        qword ptr [edi+8] 
			fsubp       st(1),st
			fmul        qword ptr [edx] 
			faddp       st(1),st
			fstp        qword ptr [detV] 

			//  double tv = detV*iDet;

			fld         qword ptr [detV]
			fmul        qword ptr [iDet]
			fstp        qword ptr [tv]

			// if( tv<0.0 || tv>1.0 ) return 0; // Bounding rectange exclusion

			fldz
			fld         qword ptr [tv]
			fcomip      st,st(1)
			jb          not_found
			fstp        qword ptr [dummy]
			fld1
			fld         qword ptr [tv]
			fcomip      st,st(1)
			ja          not_found
			fstp        qword ptr [dummy]

			//  double detW = DET33( U->x , V->x , Z->x ,
			//                       U->y , V->y , Z->y ,
			//                       U->z , V->z , Z->z );

			fld         qword ptr [ecx+8]
			fmul        qword ptr [edi+10h]
			fld         qword ptr [ecx+10h]
			fmul        qword ptr [edi+8]
			fsubp       st(1),st
			fmul        qword ptr [ebx]
			fld         qword ptr [edi+8]
			fmul        qword ptr [ebx+10h]
			fld         qword ptr [edi+10h]
			fmul        qword ptr [ebx+8]
			fsubp       st(1),st
			fmul        qword ptr [ecx]
			faddp       st(1),st
			fld         qword ptr [ebx+8]
			fmul        qword ptr [ecx+10h]
			fld         qword ptr [ebx+10h]
			fmul        qword ptr [ecx+8]
			fsubp       st(1),st
			fmul        qword ptr [edi]
			faddp       st(1),st
			fstp        qword ptr [detW]

			//  double tw = detW*iDet;

			fld         qword ptr [detW] 
			fmul        qword ptr [iDet] 
			fstp        qword ptr [tw] 

			//  if(tw<0.0) return 0; // Plane is in the back

			fldz
			fld         qword ptr [tw]
			fcomip      st,st(1)
			jb          not_found
			fstp        qword ptr [dummy]

			// Intersection found
			mov         eax,dword ptr [u]
			fld         qword ptr [tu]
			fstp        qword ptr [eax]
			mov         eax,dword ptr [v]
			fld         qword ptr [tv]
			fstp        qword ptr [eax]
			mov         eax,dword ptr [w]
			fld         qword ptr [tw]
			fstp        qword ptr [eax]

	}

	return 1;

not_found:
	__asm fstp qword ptr [dummy]  // Pop FPU stack
		return 0;
}

// SSE2 SolveI implementation: Performance are similar to FPU, not used.

// Compute determinant and store it in xmm0 low quadword (REM: all xmm registers are affected)
#define DET33SSE(_11,_12,_13,_21,_22,_23,_31,_32,_33) \
	__asm {                 \
	__asm movq    xmm0,_23  \
	__asm movq    xmm1,_31  \
	__asm movq    xmm2,_33  \
	__asm movq    xmm3,_21  \
	__asm movhpd  xmm0,_22  \
	__asm movhpd  xmm1,_33  \
	__asm movhpd  xmm2,_32  \
	__asm movhpd  xmm3,_23  \
	__asm movq    xmm4,_21  \
	__asm movq    xmm5,_32  \
	__asm mulpd   xmm0,xmm1 \
	__asm mulpd   xmm2,xmm3 \
	__asm mulsd   xmm4,xmm5 \
	__asm movq    xmm6,_31  \
	__asm movq    xmm7,_22  \
	__asm mulsd   xmm6,xmm7 \
	__asm subsd   xmm4,xmm6 \
	__asm subpd   xmm0,xmm2 \
	__asm movq    xmm1,_12  \
	__asm movhpd  xmm1,_11  \
	__asm movq    xmm2,_13  \
	__asm mulpd   xmm0,xmm1 \
	__asm mulsd   xmm4,xmm2 \
	__asm movapd  xmm2,xmm0 \
	__asm shufpd  xmm0,xmm0,1 \
	__asm addsd   xmm0,xmm2 \
	__asm addsd   xmm0,xmm4 \
}

int    SolveISSE2(double *u ,double *v,double *w,         // Result
	VERTEX3D *O,VERTEX3D *U,VERTEX3D *V,    // Rectangle (O,u,v) u in [0,1], v in [0,1]
	VERTEX3D *rayDir,VERTEX3D *rayPos)      // Ray
{
	double ux = U->x;
	double uy = U->y;
	double uz = U->z;

	double vx = V->x;
	double vy = V->y;
	double vz = V->z;

	double rx = -rayDir->x;
	double ry = -rayDir->y;
	double rz = -rayDir->z;

	double det;
	double one = 1.0;

	DET33SSE( ux , vx , rx ,
		uy , vy , ry ,
		uz , vz , rz );

	__asm {
		pxor    xmm1,xmm1
			comisd  xmm0,xmm1
			je _notfound
			movsd det,xmm0
	}

	double odx = rayPos->x - O->x;
	double ody = rayPos->y - O->y;
	double odz = rayPos->z - O->z;

	DET33SSE( odx , vx , rx ,
		ody , vy , ry ,
		odz , vz , rz );

	__asm {
		divsd   xmm0,det
			pxor    xmm1,xmm1
			comisd  xmm0,xmm1
			jb _notfound
			movq    xmm1,one
			comisd  xmm0,xmm1
			ja _notfound
			mov   edi,u
			movsd [edi],xmm0
	}

	DET33SSE( ux , odx , rx ,
		uy , ody , ry ,
		uz , odz , rz );

	__asm {
		divsd   xmm0,det
			pxor    xmm1,xmm1
			comisd  xmm0,xmm1
			jb _notfound
			movq    xmm1,one
			comisd  xmm0,xmm1
			ja _notfound
			mov   edi,v
			movsd [edi],xmm0
	}

	DET33SSE( ux , vx , odx ,
		uy , vy , ody ,
		uz , vz , odz );

	__asm {
		divsd   xmm0,det
			pxor    xmm1,xmm1
			comisd  xmm0,xmm1
			jb _notfound
			mov   edi,w
			movsd [edi],xmm0
	}

	__asm emms;
	return 1;

_notfound:
	__asm emms;
	return 0;

}

#endif
