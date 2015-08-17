//Common geometry handling/editing features, shared between Molflow and Synrad

#include "Geometry.h"
#include "GLApp\GLMessageBox.h"
#include "Molflow.h"

extern MolFlow *mApp;

void Geometry::CheckCollinear() {
	char tmp[256];
	// Check collinear polygon
	int nbCollinear = 0;
	for (int i = 0; i < GetNbFacet(); i++) {
		if (GetFacet(i)->collinear) nbCollinear++;
	}
	BOOL ok = FALSE;
	if (nbCollinear) {
		sprintf(tmp, "%d null polygon(s) found !\nThese polygons have all vertices on a single line, thus they do nothing.\nDelete them?", nbCollinear);
		ok = GLMessageBox::Display(tmp, "Info", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK;
	}
	if (ok) RemoveCollinear();
}

void Geometry::CheckNonSimple() {
	char tmp[256];
	// Check non simple polygon
	int *nonSimpleList = (int *)malloc(GetNbFacet()*sizeof(int));
	int nbNonSimple = 0;
	for (int i = 0; i < GetNbFacet(); i++) {
		if (GetFacet(i)->sh.sign == 0.0)
			nonSimpleList[nbNonSimple++] = i;
	}
	BOOL ok = FALSE;
	if (nbNonSimple) {
		sprintf(tmp, "%d non simple (or null) polygon(s) found !\nSome tasks may not work properly\n"
			"Should I try to correct them (vertex shifting)?", nbNonSimple);
		ok = GLMessageBox::Display(tmp, "Warning", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONWARNING) == GLDLG_OK;
	}

	if (ok) CorrectNonSimple(nonSimpleList, nbNonSimple);
	SAFE_FREE(nonSimpleList);
}

void Geometry::CheckIsolatedVertex() {
	int nbI = HasIsolatedVertices();
	if (nbI) {
		char tmp[256];
		sprintf(tmp, "Remove %d isolated vertices ?", nbI);
		if (GLMessageBox::Display(tmp, "Question", GLDLG_OK | GLDLG_CANCEL, GLDLG_ICONINFO) == GLDLG_OK) {
			DeleteIsolatedVertices(FALSE);
		}
	}
}

void Geometry::CorrectNonSimple(int *nonSimpleList, int nbNonSimple) {
	mApp->changedSinceSave = TRUE;
	Facet *f;
	for (int i = 0; i < nbNonSimple; i++) {
		f = GetFacet(nonSimpleList[i]);
		if (f->sh.sign == 0.0) {
			int j = 0;
			while ((j < f->sh.nbIndex) && (f->sh.sign == 0.0)) {
				f->ShiftVertex();
				InitializeGeometry(nonSimpleList[i]);
				//f->DetectOrientation();
				j++;
			}
		}
	}
	//BuildGLList();
}

int Geometry::GetNbVertex() {
	return sh.nbVertex;
}

VERTEX3D Geometry::GetFacetCenter(int facet) {

	return facets[facet]->sh.center;

}

int Geometry::GetNbStructure() {
	return sh.nbSuper;
}

char *Geometry::GetStructureName(int idx) {
	return strName[idx];
}

void Geometry::CreatePolyFromVertices_Convex() {
	//creates facet from selected vertices

	mApp->changedSinceSave = TRUE;
	nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx, 0xFF, sh.nbVertex*sizeof(int));
	for (int i = 0; i < sh.nbVertex; i++) {
		//VERTEX3D *v = GetVertex(i);
		if (vertices3[i].selected) {
			vIdx[nbSelectedVertex] = i;
			nbSelectedVertex++;
		}
	}

	if (nbSelectedVertex < 3) {
		char errMsg[512];
		sprintf(errMsg, "Select at least 3 vertices.");
		throw Error(errMsg);
		return;
	}//at least three vertices

	VERTEX3D U, V, N;
	U.x = vertices3[vIdx[0]].x - vertices3[vIdx[1]].x;
	U.y = vertices3[vIdx[0]].y - vertices3[vIdx[1]].y;
	U.z = vertices3[vIdx[0]].z - vertices3[vIdx[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U, 1.0 / nU); // Normalize U

	int i2 = 2;
	double nV;
	do {
		V.x = vertices3[vIdx[0]].x - vertices3[vIdx[i2]].x;
		V.y = vertices3[vIdx[0]].y - vertices3[vIdx[i2]].y;
		V.z = vertices3[vIdx[0]].z - vertices3[vIdx[i2]].z;
		nV = Norme(&V);
		ScalarMult(&V, 1.0 / nV); // Normalize V
		i2++;
	} while (Dot(&U, &V) > 0.99 && i2 < nbSelectedVertex); //if U and V are almost the same, the projection would be inaccurate


	//Now we have the U,V plane, let's define it by computing the normal vector:
	Cross(&N, &V, &U); //We have a normal vector
	double nN = Norme(&N);
	ScalarMult(&N, 1.0 / nN); // Normalize N

	Cross(&V, &N, &U); //Make V perpendicular to U and N (and still in the U,V plane)
	nV = Norme(&V);
	ScalarMult(&V, 1.0 / nV); // Normalize V

	VERTEX2D *projected = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));
	VERTEX2D *debug = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));

	//Get coordinates in the U,V system
	for (int i = 0; i < nbSelectedVertex; i++) {
		ProjectVertex(&(vertices3[vIdx[i]]), &(projected[i]), U, V, vertices3[vIdx[0]]);
	}

	//Graham scan here on the projected[] array
	int *returnList = (int *)malloc(nbSelectedVertex*sizeof(int));
	grahamMain(projected, nbSelectedVertex, returnList);
	int ii, loopLength;
	for (ii = 0; ii<nbSelectedVertex; ii++) {
		if (returnList[ii] == returnList[0] && ii>0) break;
	}
	loopLength = ii;
	//End graham scan


	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet - 1] = new Facet(loopLength);
	facets[sh.nbFacet - 1]->sh.sticking = 0.0;
	facets[sh.nbFacet - 1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet - 1]->selected = TRUE;
	nbSelected = 1;
	for (int i = 0; i < loopLength; i++) {
		facets[sh.nbFacet - 1]->indices[i] = vIdx[returnList[i]];
	}
	SAFE_FREE(vIdx);
	SAFE_FREE(projected);
	SAFE_FREE(debug);
	SAFE_FREE(returnList);

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet - 1);
	mApp->facetList->ScrollToVisible(sh.nbFacet - 1, 1, FALSE);
}

void Geometry::CreatePolyFromVertices_Order() {
	//creates facet from selected vertices

	mApp->changedSinceSave = TRUE;
	/*nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx,0xFF,sh.nbVertex*sizeof(int));
	for(int i=0;i<sh.nbVertex;i++ ) {
	//VERTEX3D *v = GetVertex(i);
	if( vertices3[i].selected ) {
	vIdx[nbSelectedVertex] = i;
	nbSelectedVertex++;
	}
	}

	if (selectedVertexList.size()<3) {
	char errMsg[512];
	sprintf(errMsg,"Select at least 3 vertices.");
	throw Error(errMsg);
	return;
	}//at least three vertices

	VERTEX3D U,V,N;
	U.x = vertices3[selectedVertexList[0]].x - vertices3[selectedVertexList[1]].x;
	U.y = vertices3[selectedVertexList[0]].y - vertices3[selectedVertexList[1]].y;
	U.z = vertices3[selectedVertexList[0]].z - vertices3[selectedVertexList[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U,1.0/nU); // Normalize U

	int i2=2;
	double nV;
	do {
	V.x = vertices3[selectedVertexList[0]].x - vertices3[selectedVertexList[i2]].x;
	V.y = vertices3[selectedVertexList[0]].y - vertices3[selectedVertexList[i2]].y;
	V.z = vertices3[selectedVertexList[0]].z - vertices3[selectedVertexList[i2]].z;
	nV = Norme(&V);
	ScalarMult(&V,1.0/nV); // Normalize V
	i2++;
	} while (Dot(&U,&V)>0.99 && i2<nbSelectedVertex); //if U and V are almost the same, the projection would be inaccurate


	//Now we have the U,V plane, let's define it by computing the normal vector:
	Cross(&N,&V,&U); //We have a normal vector
	double nN = Norme(&N);
	ScalarMult(&N,1.0/nN); // Normalize N

	Cross(&V,&N,&U); //Make V perpendicular to U and N (and still in the U,V plane)
	nV = Norme(&V);
	ScalarMult(&V,1.0/nV); // Normalize V

	VERTEX2D *projected = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));
	VERTEX2D *debug = (VERTEX2D *)malloc(nbSelectedVertex * sizeof(VERTEX2D));

	//Get coordinates in the U,V system
	for (int i=0;i<nbSelectedVertex;i++) {
	ProjectVertex(&(vertices3[vIdx[i]]),&(projected[i]),U,V,vertices3[vIdx[0]]);
	}

	//Graham scan here on the projected[] array
	int *returnList = (int *)malloc(nbSelectedVertex*sizeof(int));
	grahamMain(projected,nbSelectedVertex,returnList);
	int ii,loopLength;
	for (ii=0;ii<nbSelectedVertex;ii++) {
	if (returnList[ii]==returnList[0] && ii>0) break;
	}
	loopLength=ii;
	//End graham scan
	*/

	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet - 1] = new Facet((int)selectedVertexList.size());
	facets[sh.nbFacet - 1]->sh.sticking = 0.0;
	facets[sh.nbFacet - 1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet - 1]->selected = TRUE;
	nbSelected = 1;
	for (size_t i = 0; i < selectedVertexList.size(); i++) {
		facets[sh.nbFacet - 1]->indices[i] = selectedVertexList[i];
	}
	/*SAFE_FREE(vIdx);
	SAFE_FREE(projected);
	SAFE_FREE(debug);
	SAFE_FREE(returnList);*/

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet - 1);
	mApp->facetList->ScrollToVisible(sh.nbFacet - 1, 1, FALSE);
}

void Geometry::CreateDifference() {
	//creates facet from selected vertices

	mApp->changedSinceSave = TRUE;
	nbSelectedVertex = 0;

	if (nbSelected != 2) {
		char errMsg[512];
		sprintf(errMsg, "Select exactly 2 facets.");
		throw Error(errMsg);
		return;
	}//at least three vertices

	int firstFacet = -1;
	int secondFacet = -1;;
	for (int i = 0; i < sh.nbFacet && secondFacet < 0; i++)
	{
		if (facets[i]->selected) {
			if (firstFacet < 0) firstFacet = i;
			else (secondFacet = i);
		}
	}

	//TO DO:
	//swap if normals not collinear
	//shift vertex to make nice cut

	//a new facet
	sh.nbFacet++;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	facets[sh.nbFacet - 1] = new Facet(facets[firstFacet]->sh.nbIndex + facets[secondFacet]->sh.nbIndex + 2);
	facets[sh.nbFacet - 1]->sh.sticking = 0.0;
	facets[sh.nbFacet - 1]->sh.sticking = DES_NONE;
	//set selection
	UnSelectAll();
	facets[sh.nbFacet - 1]->selected = TRUE;
	nbSelected = 1;
	//one circle on first facet
	int counter = 0;
	for (int i = 0; i < facets[firstFacet]->sh.nbIndex; i++)
		facets[sh.nbFacet - 1]->indices[counter++] = facets[firstFacet]->indices[i];
	//close circle by adding the first vertex again
	facets[sh.nbFacet - 1]->indices[counter++] = facets[firstFacet]->indices[0];
	//reverse circle on second facet
	for (int i = facets[secondFacet]->sh.nbIndex - 1; i >= 0; i--)
		facets[sh.nbFacet - 1]->indices[counter++] = facets[secondFacet]->GetIndex(i + 1);
	//close circle by adding the first vertex again
	facets[sh.nbFacet - 1]->indices[counter++] = facets[secondFacet]->indices[0];

	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet - 1);
	mApp->facetList->ScrollToVisible(sh.nbFacet - 1, 1, FALSE);
}

void Geometry::SelectCoplanar(int width, int height, double tolerance) {


	nbSelectedVertex = 0;

	int *vIdx = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(vIdx, 0xFF, sh.nbVertex*sizeof(int));
	for (int i = 0; i < sh.nbVertex; i++) {
		//VERTEX3D *v = GetVertex(i);
		if (vertices3[i].selected) {
			vIdx[nbSelectedVertex] = i;
			nbSelectedVertex++;
		}
	}

	VERTEX3D U, V, N;
	U.x = vertices3[vIdx[0]].x - vertices3[vIdx[1]].x;
	U.y = vertices3[vIdx[0]].y - vertices3[vIdx[1]].y;
	U.z = vertices3[vIdx[0]].z - vertices3[vIdx[1]].z;
	double nU = Norme(&U);
	ScalarMult(&U, 1.0 / nU); // Normalize U

	V.x = vertices3[vIdx[0]].x - vertices3[vIdx[2]].x;
	V.y = vertices3[vIdx[0]].y - vertices3[vIdx[2]].y;
	V.z = vertices3[vIdx[0]].z - vertices3[vIdx[2]].z;
	double nV = Norme(&V);
	ScalarMult(&V, 1.0 / nV); // Normalize V

	Cross(&N, &V, &U); //We have a normal vector
	double nN = Norme(&N);
	if (nN < 1e-8) {
		GLMessageBox::Display("Sorry, the 3 selected vertices are on a line.", "Can't define plane", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	ScalarMult(&N, 1.0 / nN); // Normalize N

	// Plane equation
	double A = N.x;
	double B = N.y;
	double C = N.z;
	VERTEX3D p0 = vertices3[vIdx[0]];
	double D = -Dot(&N, &p0);

	//double denominator=sqrt(pow(A,2)+pow(B,2)+pow(C,2));
	double distance;

	int outX, outY;

	for (int i = 0; i < sh.nbVertex; i++) {
		VERTEX3D *v = GetVertex(i);
		BOOL onScreen = GLToolkit::Get2DScreenCoord((float)v->x, (float)v->y, (float)v->z, &outX, &outY); //To improve
		onScreen = (onScreen && outX >= 0 && outY >= 0 && outX <= (width) && (outY <= height));
		if (onScreen) {
			distance = abs(A*v->x + B*v->y + C*v->z + D);
			if (distance < tolerance) { //vertex is on the plane

				vertices3[i].selected = TRUE;

			}
			else {
				vertices3[i].selected = FALSE;
			}
		}
		else {
			vertices3[i].selected = FALSE;
		}
	}
	SAFE_FREE(vIdx);
}

VERTEX3D *Geometry::GetVertex(int idx) {
	return vertices3 + idx;
}

Facet *Geometry::GetFacet(int facet) {
	if (facet >= sh.nbFacet || facet < 0) {
		char errMsg[512];
		sprintf(errMsg, "Geometry::GetFacet()\nA process tried to access facet #%d that doesn't exist.\nAutoSaving and probably crashing...", facet + 1);
		GLMessageBox::Display(errMsg, "Error", GLDLG_OK, GLDLG_ICONERROR);
		mApp->AutoSave();
		throw Error(errMsg);
	}
	return facets[facet];
}

int Geometry::GetNbFacet() {
	return sh.nbFacet;
}

AABB Geometry::GetBB() {

	if (viewStruct < 0) {

		return bb;

	}
	else {

		// BB of selected struture
		AABB sbb;

		sbb.min.x = 1e100;
		sbb.min.y = 1e100;
		sbb.min.z = 1e100;
		sbb.max.x = -1e100;
		sbb.max.y = -1e100;
		sbb.max.z = -1e100;

		// Axis Aligned Bounding Box
		for (int i = 0; i < sh.nbFacet; i++) {
			Facet *f = facets[i];
			if (f->sh.superIdx == viewStruct) {
				for (int j = 0; j < f->sh.nbIndex; j++) {
					VERTEX3D p = vertices3[f->indices[j]];
					if (p.x < sbb.min.x) sbb.min.x = p.x;
					if (p.y < sbb.min.y) sbb.min.y = p.y;
					if (p.z < sbb.min.z) sbb.min.z = p.z;
					if (p.x > sbb.max.x) sbb.max.x = p.x;
					if (p.y > sbb.max.y) sbb.max.y = p.y;
					if (p.z > sbb.max.z) sbb.max.z = p.z;
				}
			}
		}



		return sbb;
	}

}

VERTEX3D Geometry::GetCenter() {

	if (viewStruct < 0) {

		return center;


	}
	else {

		VERTEX3D r;
		AABB sbb = GetBB();
		r.x = (sbb.max.x + sbb.min.x) / 2.0;
		r.y = (sbb.max.y + sbb.min.y) / 2.0;
		r.z = (sbb.max.z + sbb.min.z) / 2.0;

		return r;

	}
}

int Geometry::AddRefVertex(VERTEX3D *p, VERTEX3D *refs, int *nbRef) {

	BOOL found = FALSE;
	int i = 0;
	//VERTEX3D n;
	double v2 = vThreshold*vThreshold;

	while (i < *nbRef && !found) {
		//Sub(&n,p,refs + i);
		double dx = abs(p->x - (refs + i)->x);
		if (dx < vThreshold) {
			double dy = abs(p->y - (refs + i)->y);
			if (dy < vThreshold) {
				double dz = abs(p->z - (refs + i)->z);
				if (dz < vThreshold) {
					found = (dx*dx + dy*dy + dz*dz < v2);
				}
			}
		}
		if (!found) i++;
	}

	if (!found) {
		// Add a new reference vertex
		refs[*nbRef] = *p;
		*nbRef = *nbRef + 1;
	}

	return i;

}

void Geometry::CollapseVertex(GLProgress *prg, double totalWork) {
	mApp->changedSinceSave = TRUE;
	if (!isLoaded) return;
	// Collapse neighbor vertices
	VERTEX3D *refs = (VERTEX3D *)malloc(sh.nbVertex*sizeof(VERTEX3D));
	if (!refs) throw Error("Out of memory: CollapseVertex");
	int      *idx  = (int *)malloc(sh.nbVertex*sizeof(int));
	if (!idx) throw Error("Out of memory: CollapseVertex");
	int       nbRef=0;

	// Collapse
	prg->SetMessage("Collapsing vertices...");
	for (int i = 0; i < sh.nbVertex; i++) {
		prg->SetProgress(((double)i / (double)sh.nbVertex) / totalWork);
		idx[i] = AddRefVertex(vertices3 + i, refs, &nbRef);
	}

	// Create the new vertex array
	SAFE_FREE(vertices3);
	vertices3 = (VERTEX3D *)malloc(nbRef * sizeof(VERTEX3D));
	if (!vertices3) throw Error("Out of memory: CollapseVertex");
	//UnselectAllVertex();

	memcpy(vertices3, refs, nbRef * sizeof(VERTEX3D));
	sh.nbVertex = nbRef;

	// Update facets indices
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		prg->SetProgress(((double)i / (double)sh.nbFacet) * 0.05 + 0.45);
		for (int j = 0; j < f->sh.nbIndex; j++)
			f->indices[j] = idx[f->indices[j]];
	}

	free(idx);
	free(refs);

}

BOOL Geometry::GetCommonEdges(Facet *f1, Facet *f2, int *c1, int *c2, int *chainLength) {

	// Detect common edge between facet
	int p11, p12, p21, p22, lgth, si, sj;
	int maxLength = 0;

	for (int i = 0; i < f1->sh.nbIndex; i++) {

		p11 = f1->GetIndex(i);
		p12 = f1->GetIndex(i + 1);

		for (int j = 0; j < f2->sh.nbIndex; j++) {

			p21 = f2->GetIndex(j);
			p22 = f2->GetIndex(j + 1);

			if (p11 == p22 && p12 == p21) {

				// Common edge found
				si = i;
				sj = j;
				lgth = 1;

				// Loop until the end of the common edge chain
				i += 2;
				j -= 1;
				p12 = f1->GetIndex(i);
				p21 = f2->GetIndex(j);
				BOOL ok = (p12 == p21);
				while (lgth < f1->sh.nbIndex && lgth<f2->sh.nbIndex && ok) {
					p12 = f1->GetIndex(i);
					p21 = f2->GetIndex(j);
					ok = (p12 == p21);
					if (ok) {
						i++; j--;
						lgth++;
					}
				}

				if (lgth>maxLength) {
					*c1 = si;
					*c2 = sj;
					maxLength = lgth;
				}

			}

		}

	}

	if (maxLength > 0) {
		*chainLength = maxLength;
		return TRUE;
	}

	return FALSE;

}

void Geometry::MoveVertexTo(int idx, double x, double y, double z) {
	vertices3[idx].x = x;
	vertices3[idx].y = y;
	vertices3[idx].z = z;

}

void Geometry::SwapNormal() {

	if (!IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (GetNbSelected() <= 0) return;
	mApp->changedSinceSave = TRUE;
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->selected) {
			f->SwapNormal();
			InitializeGeometry(i);// Reinitialise geom
			try {
				SetFacetTexture(i, f->tRatio, f->hasMesh);
			} catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
			}
		}
	}

	DeleteGLLists(TRUE, TRUE);
	BuildGLList();

}

void Geometry::Extrude(VERTEX3D offset,double distance) {

	//creates facet from selected vertices

	mApp->changedSinceSave = TRUE;

	int oldNbFacet = sh.nbFacet;
	for (int i = 0; i < oldNbFacet; i++)
	{
		if (facets[i]->selected) {
			int firstFacet = i;
			facets[firstFacet]->selected = FALSE;
			//Update direction if necessary
			VERTEX3D dir2;

			if (IS_ZERO(offset.x) && IS_ZERO(offset.y) && IS_ZERO(offset.z)) { //Use facet normal to determine offset
				dir2.x = facets[firstFacet]->sh.N.x*distance;
				dir2.y = facets[firstFacet]->sh.N.y*distance;
				dir2.z = facets[firstFacet]->sh.N.z*distance;
			}
			else { //Use provided offset
				dir2.x = offset.x;
				dir2.y = offset.y;
				dir2.z = offset.z;
			}

			//Copy facet
			VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((sh.nbVertex + facets[firstFacet]->sh.nbIndex) * sizeof(VERTEX3D)); //create new, extended vertex array
			memmove(tmp_vertices3, vertices3, (sh.nbVertex)*sizeof(VERTEX3D)); //copy old vertices
			SAFE_FREE(vertices3); //delete old array
			vertices3 = tmp_vertices3; //make new array the official vertex holder


			for (int j = 0; j < facets[firstFacet]->sh.nbIndex; j++) { //copy vertex coordinates
				vertices3[sh.nbVertex + j].x = vertices3[facets[firstFacet]->indices[j]].x + dir2.x;
				vertices3[sh.nbVertex + j].y = vertices3[facets[firstFacet]->indices[j]].y + dir2.y;
				vertices3[sh.nbVertex + j].z = vertices3[facets[firstFacet]->indices[j]].z + dir2.z;
				vertices3[sh.nbVertex + j].selected = FALSE;
			}

			//Copy facet
			int secondFacet = sh.nbFacet; //last facet
			facets = (Facet **)realloc(facets, (sh.nbFacet+1+facets[firstFacet]->sh.nbIndex) * sizeof(Facet *));
			facets[secondFacet] = new Facet(facets[firstFacet]->sh.nbIndex);
			facets[secondFacet]->selected = TRUE;
			for (int j = 0; j < facets[firstFacet]->sh.nbIndex; j++)
				facets[secondFacet]->indices[facets[firstFacet]->sh.nbIndex-1-j] = sh.nbVertex + j; //assign new vertices to new facet in inverse order

			int direction=1;
			if (Dot(&dir2, &facets[firstFacet]->sh.N) *distance< 0.0) direction *= -1; //extrusion towards normal or opposite?
			for (int j = 0; j < facets[firstFacet]->sh.nbIndex; j++) {
				facets[secondFacet + 1 + j] = new Facet(4);
				facets[secondFacet + 1 + j]->indices[0] = facets[firstFacet]->indices[j%facets[firstFacet]->sh.nbIndex];
				facets[secondFacet + 1 + j]->indices[1] = facets[secondFacet]->indices[(2*facets[firstFacet]->sh.nbIndex - 1 - j) %facets[firstFacet]->sh.nbIndex];
				facets[secondFacet + 1 + j]->indices[2] = facets[secondFacet]->indices[(2*facets[firstFacet]->sh.nbIndex - 2 - j) % facets[firstFacet]->sh.nbIndex];
				facets[secondFacet + 1 + j]->indices[3] = facets[firstFacet]->indices[(j + 1)%facets[firstFacet]->sh.nbIndex];
				facets[secondFacet + 1 + j]->selected = TRUE;
			}
			sh.nbVertex += facets[firstFacet]->sh.nbIndex; //update number of vertices
			sh.nbFacet += facets[firstFacet]->sh.nbIndex + 1;
		}
	}

	
	InitializeGeometry();
	mApp->UpdateFacetParams(TRUE);
	UpdateSelection();
	mApp->facetList->SetSelectedRow(sh.nbFacet - 1);
	mApp->facetList->ScrollToVisible(sh.nbFacet - 1, 1, FALSE);

}

void Geometry::ShiftVertex() {

	if (!IsLoaded()) {
		GLMessageBox::Display("No geometry loaded.", "No geometry", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}
	if (GetNbSelected() <= 0) return;
	mApp->changedSinceSave = TRUE;
	DeleteGLLists(TRUE, TRUE);
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		if (f->selected) {
			f->ShiftVertex();
			InitializeGeometry(i);// Reinitialise geom
			try {
				SetFacetTexture(i, f->tRatio, f->hasMesh);
			}
			catch (Error &e) {
				GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
			}

		}
	}
	
	BuildGLList();
}

void Geometry::Merge(int nbV, int nbF, VERTEX3D *nV, Facet **nF) {
	mApp->changedSinceSave = TRUE;
	// Merge the current geometry with the specified one
	if (!nbV || !nbF) return;

	// Reallocate mem
	Facet   **nFacets = (Facet **)malloc((sh.nbFacet + nbF) * sizeof(Facet *));
	VERTEX3D *nVertices3 = (VERTEX3D *)malloc((sh.nbVertex + nbV) * sizeof(VERTEX3D));


	if (sh.nbFacet) memcpy(nFacets, facets, sizeof(Facet *)* sh.nbFacet);
	memcpy(nFacets + sh.nbFacet, nF, sizeof(Facet *)* nbF);

	if (sh.nbVertex) memcpy(nVertices3, vertices3, sizeof(VERTEX3D)* sh.nbVertex);
	memcpy(nVertices3 + sh.nbVertex, nV, sizeof(VERTEX3D)* nbV);

	SAFE_FREE(facets);
	SAFE_FREE(vertices3);
	facets = nFacets;
	vertices3 = nVertices3;
	//UnselectAllVertex();

	// Shift indices
	for (int i = 0; i < nbF; i++) {
		Facet *f = facets[sh.nbFacet + i];
		for (int j = 0; j < f->sh.nbIndex; j++)
			f->indices[j] += sh.nbVertex;
	}

	sh.nbVertex += nbV;
	sh.nbFacet += nbF;

}

void Geometry::RemoveLinkFacet() { //unused

	// Remove facet used as link for superstructure (not needed)
	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->IsLinkFacet()) nb++;

	if (nb == 0) return;

	Facet **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->IsLinkFacet()) {
			delete facets[i];

		}
		else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

}

int Geometry::HasIsolatedVertices() {

	// Check if there are unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check, 0, sh.nbVertex*sizeof(int));

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		for (int j = 0; j < f->sh.nbIndex; j++) {
			check[f->indices[j]]++;
		}
	}

	int nbUnused = 0;
	for (int i = 0; i < sh.nbVertex; i++) {
		if (!check[i]) nbUnused++;
	}

	SAFE_FREE(check);
	return nbUnused;

}

void  Geometry::DeleteIsolatedVertices(BOOL selectedOnly) {
	mApp->changedSinceSave = TRUE;
	// Remove unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check, 0, sh.nbVertex*sizeof(int));

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		for (int j = 0; j < f->sh.nbIndex; j++) {
			check[f->indices[j]]++;
		}
	}

	int nbUnused = 0;
	for (int i = 0; i < sh.nbVertex; i++) {
		if (!check[i] && !(selectedOnly&&!vertices3[i].selected)) nbUnused++;
	}

	int nbVert = sh.nbVertex - nbUnused;


	if (nbVert == 0) {
		// Remove all
		SAFE_FREE(check);
		Clear();
		return;
	}


	// Update facet indices
	int *newId = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(newId, 0, sh.nbVertex*sizeof(int));
	for (int i = 0, n = 0; i < sh.nbVertex; i++) {
		if (check[i] || (selectedOnly && !vertices3[i].selected)) {
			newId[i] = n;
			n++;
		}
	}
	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		for (int j = 0; j < f->sh.nbIndex; j++) {
			f->indices[j] = newId[f->indices[j]];
		}
	}



	VERTEX3D *nVert = (VERTEX3D *)malloc(nbVert*sizeof(VERTEX3D));

	_ASSERTE(nVert);
	for (int i = 0, n = 0; i < sh.nbVertex; i++) {
		if (check[i] || (selectedOnly && !vertices3[i].selected)) {
			nVert[n] = vertices3[i];
			n++;
		}
	}

	SAFE_FREE(vertices3);
	vertices3 = nVert;
	sh.nbVertex = nbVert;

	SAFE_FREE(check);
	SAFE_FREE(newId);

	// Delete old resources
	//DeleteGLLists(TRUE,TRUE);

	//InitializeGeometry();

}

void  Geometry::SelectIsolatedVertices() {

	UnselectAllVertex();
	// Select unused vertices
	int *check = (int *)malloc(sh.nbVertex*sizeof(int));
	memset(check, 0, sh.nbVertex*sizeof(int));

	for (int i = 0; i < sh.nbFacet; i++) {
		Facet *f = facets[i];
		for (int j = 0; j < f->sh.nbIndex; j++) {
			check[f->indices[j]]++;
		}
	}

	for (int i = 0; i < sh.nbVertex; i++) {
		if (!check[i]) vertices3[i].selected = TRUE;
	}

	SAFE_FREE(check);
}

void Geometry::RemoveCollinear() {

	mApp->changedSinceSave = TRUE;
	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->collinear) nb++;

	if (nb == 0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->collinear) {
			delete facets[i];
			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);

		}
		else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	BuildGLList();
	mApp->UpdateModelParams();

}

void Geometry::RemoveSelectedVertex() {

	mApp->changedSinceSave = TRUE;

	//Analyze facets
	std::vector<int> facetsToRemove, facetsToChange;
	for (int f = 0; f < sh.nbFacet; f++) {
		int nbSelVertex = 0;
		for (int i = 0; i < facets[f]->sh.nbIndex; i++)
			if (vertices3[facets[f]->indices[i]].selected)
				nbSelVertex++;
		if (nbSelVertex) {
			facetsToChange.push_back(f);
			if ((facets[f]->sh.nbIndex - nbSelVertex) <= 2)
				facetsToRemove.push_back(f);
		}
	}

	for (size_t c = 0; c < facetsToChange.size(); c++) {
		Facet* f = facets[facetsToChange[c]];
		int nbRemove = 0;
		for (size_t i = 0; (int)i < f->sh.nbIndex; i++) //count how many to remove			
			if (vertices3[f->indices[i]].selected)
				nbRemove++;
		int *newIndices = (int *)malloc((f->sh.nbIndex - nbRemove)*sizeof(int));
		int nb = 0;
		for (size_t i = 0; (int)i < f->sh.nbIndex; i++)
			if (!vertices3[f->indices[i]].selected) newIndices[nb++] = f->indices[i];

		SAFE_FREE(f->indices); f->indices = newIndices;
		SAFE_FREE(f->vertices2);
		SAFE_FREE(f->visible);
		f->sh.nbIndex -= nbRemove;
		f->vertices2 = (VERTEX2D *)malloc(f->sh.nbIndex*sizeof(VERTEX2D));
		memset(f->vertices2, 0, f->sh.nbIndex * sizeof(VERTEX2D));
		f->visible = (BOOL *)malloc(f->sh.nbIndex*sizeof(BOOL));
		_ASSERTE(f->visible);
		memset(f->visible, 0xFF, f->sh.nbIndex*sizeof(BOOL));
	}

	if (facetsToRemove.size()) {
		Facet   **newFacets = (Facet **)malloc((sh.nbFacet - facetsToRemove.size()) * sizeof(Facet *));
		size_t nextToRemove = 0;
		size_t nextToAdd = 0;
		for (size_t f = 0; (int)f < sh.nbFacet; f++) {
			if (nextToRemove < facetsToRemove.size() && f == facetsToRemove[nextToRemove]) {
				delete facets[f];
				mApp->RenumberSelections(nextToAdd);
				mApp->RenumberFormulas(nextToAdd);
				nextToRemove++;

			}
			else {
				newFacets[nextToAdd++] = facets[f];
			}
		}
		SAFE_DELETE(facets);
		facets = newFacets;
		sh.nbFacet -= facetsToRemove.size();
	}

	DeleteIsolatedVertices(TRUE);

	DeleteGLLists(TRUE, TRUE);

	BuildGLList();
	CalcTotalOutGassing();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::RemoveSelected() {

	mApp->changedSinceSave = TRUE;
	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->selected) nb++;

	if (nb == 0) return;
	/*
	if(sh.nbFacet-nb==0) {
	// Remove all
	Clear();
	return;
	}
	*/

	Facet   **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			delete facets[i];
			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);

		}
		else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	BuildGLList();
	CalcTotalOutGassing();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

int Geometry::ExplodeSelected(BOOL toMap, int desType, double exponent, double *values) {


	mApp->changedSinceSave = TRUE;
	if (nbSelected == 0) return -1;

	// Check that all facet has a mesh
	BOOL ok = TRUE;
	int idx = 0;
	while (ok && idx < sh.nbFacet) {
		if (facets[idx]->selected)
			ok = facets[idx]->hasMesh;
		idx++;
	}
	if (!ok) return -2;

	int nb = 0;
	int FtoAdd = 0;
	int VtoAdd = 0;
	Facet::FACETGROUP *blocks = (Facet::FACETGROUP *)malloc(nbSelected * sizeof(Facet::FACETGROUP));

	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			facets[i]->Explode(blocks + nb);
			FtoAdd += blocks[nb].nbF;
			VtoAdd += blocks[nb].nbV;
			nb++;
		}
	}

	// Update vertex array
	VERTEX3D *ptrVert;
	int       vIdx;
	VERTEX3D *nVert = (VERTEX3D *)malloc((sh.nbVertex + VtoAdd)*sizeof(VERTEX3D));
	memcpy(nVert, vertices3, sh.nbVertex*sizeof(VERTEX3D));

	ptrVert = nVert + sh.nbVertex;
	vIdx = sh.nbVertex;
	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			facets[i]->FillVertexArray(ptrVert);
			for (int j = 0; j < blocks[nb].nbF; j++) {
				for (int k = 0; k < blocks[nb].facets[j]->sh.nbIndex; k++) {
					blocks[nb].facets[j]->indices[k] = vIdx + k;
				}
				vIdx += blocks[nb].facets[j]->sh.nbIndex;
			}
			ptrVert += blocks[nb].nbV;
			nb++;
		}
	}
	SAFE_FREE(vertices3);
	vertices3 = nVert;
	for (int i = sh.nbVertex; i < sh.nbVertex + VtoAdd; i++)
		vertices3[i].selected = FALSE;
	sh.nbVertex += VtoAdd;

	// Update facet
	Facet   **f = (Facet **)malloc((sh.nbFacet + FtoAdd - nbSelected) * sizeof(Facet *));

	// Delete selected
	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			delete facets[i];
			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);

		}
		else {
			f[nb++] = facets[i];
		}
	}

	// Add new facets
	int count = 0;
	for (int i = 0; i < nbSelected; i++) {
		for (int j = 0; j<blocks[i].nbF; j++) {
			f[nb++] = blocks[i].facets[j];
			if (toMap) { //set outgassing values
				f[nb - 1]->sh.flow = *(values + count++) *0.100; //0.1: mbar*l/s->Pa*m3/s
				if (f[nb - 1]->sh.flow>0.0) {
					f[nb - 1]->sh.desorbType = desType + 1;
					f[nb - 1]->selected = TRUE;
					if (f[nb - 1]->sh.desorbType == DES_COSINE_N) f[nb - 1]->sh.desorbTypeN = exponent;
				}
				else {
					f[nb - 1]->sh.desorbType = DES_NONE;
					f[nb - 1]->selected = FALSE;
				}
			}
		}
	}

	// Free allocated memory
	for (int i = 0; i < nbSelected; i++) {
		SAFE_FREE(blocks[i].facets);
	}
	SAFE_FREE(blocks);

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	InitializeGeometry();

	return 0;

}

void Geometry::RemoveNullFacet() {

	// Remove degenerated facet (area~0.0)
	int nb = 0;
	double areaMin = vThreshold*vThreshold;
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->sh.area < areaMin) nb++;

	if (nb == 0) return;

	Facet   **f = (Facet **)malloc((sh.nbFacet - nb) * sizeof(Facet *));

	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->sh.area < areaMin) {
			delete facets[i];

			mApp->RenumberSelections(nb);
			mApp->RenumberFormulas(nb);

		}
		else {
			f[nb++] = facets[i];
		}
	}

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	BuildGLList();

}

void Geometry::AlignFacets(int* selection, int nbSelected, int Facet_source, int Facet_dest, int Anchor_source, int Anchor_dest,
	int Aligner_source, int Aligner_dest, BOOL invertNormal, BOOL invertDir1, BOOL invertDir2, BOOL copy, Worker *worker) {

	double counter = 0.0;
	double selected = (double)GetNbSelected();
	if (selected<1E-30) return;
	GLProgress *prgAlign = new GLProgress("Aligning facets...", "Please wait");
	prgAlign->SetProgress(0.0);
	prgAlign->SetVisible(TRUE);
	if (!mApp->AskToReset(worker)) return;
	if (copy) CloneSelectedFacets(); //move
	BOOL *alreadyMoved = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyMoved, FALSE, sh.nbVertex*sizeof(BOOL*));

	//Translating facets to align anchors
	int temp;
	if (invertDir1) { //change anchor and direction on source
		temp = Aligner_source;
		Aligner_source = Anchor_source;
		Anchor_source = temp;
	}
	if (invertDir2) { //change anchor and direction on destination
		temp = Aligner_dest;
		Aligner_dest = Anchor_dest;
		Anchor_dest = temp;
	}
	VERTEX3D Translation;
	Sub(&Translation, GetVertex(Anchor_dest), GetVertex(Anchor_source));

	int nb = 0;
	for (int i = 0; i < nbSelected; i++) {
		counter += 0.333;
		prgAlign->SetProgress(counter / selected);
		for (int j = 0; j < facets[selection[i]]->sh.nbIndex; j++) {
			if (!alreadyMoved[facets[selection[i]]->indices[j]]) {
				vertices3[facets[selection[i]]->indices[j]].x += Translation.x;
				vertices3[facets[selection[i]]->indices[j]].y += Translation.y;
				vertices3[facets[selection[i]]->indices[j]].z += Translation.z;
				alreadyMoved[facets[selection[i]]->indices[j]] = TRUE;
			}
		}
	}

	SAFE_FREE(alreadyMoved);

	//Rotating to match normal vectors
	VERTEX3D Axis;
	VERTEX3D Normal;
	double angle;
	Normal = facets[Facet_dest]->sh.N;
	if (invertNormal) ScalarMult(&Normal, -1.0);
	Cross(&Axis, &(facets[Facet_source]->sh.N), &Normal);
	if (Norme(&Axis)<1e-5) { //The two normals are either collinear or the opposite
		if ((Dot(&(facets[Facet_dest]->sh.N), &(facets[Facet_source]->sh.N))>0.99999 && (!invertNormal)) ||
			(Dot(&(facets[Facet_dest]->sh.N), &(facets[Facet_source]->sh.N)) < 0.00001 && invertNormal)) { //no rotation needed
			Axis.x = 1.0;
			Axis.y = 0.0;
			Axis.z = 0.0;
			angle = 0.0;

		}
		else { //180deg rotation needed
			Cross(&Axis, &(facets[Facet_source]->sh.U), &(facets[Facet_source]->sh.N));
			Normalize(&Axis);
			angle = 180.0;
		}

	}
	else {
		Normalize(&Axis);
		angle = Dot(&(facets[Facet_source]->sh.N), &(facets[Facet_dest]->sh.N)) / Norme(&(facets[Facet_source]->sh.N)) / Norme(&(facets[Facet_dest]->sh.N));
		//BOOL opposite=(angle<0.0);
		angle = acos(angle);
		angle = angle / PI * 180;
		if (invertNormal) angle = 180.0 - angle;
	}


	BOOL *alreadyRotated = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyRotated, FALSE, sh.nbVertex*sizeof(BOOL*));


	nb = 0;
	for (int i = 0; i < nbSelected; i++) {
		counter += 0.333;
		prgAlign->SetProgress(counter / selected);
		for (int j = 0; j < facets[selection[i]]->sh.nbIndex; j++) {
			if (!alreadyRotated[facets[selection[i]]->indices[j]]) {
				//rotation comes here
				Rotate(&(vertices3[facets[selection[i]]->indices[j]]), *(GetVertex(Anchor_dest)), Axis, angle);
				alreadyRotated[facets[selection[i]]->indices[j]] = TRUE;
			}
		}
	}

	SAFE_FREE(alreadyRotated);

	//Rotating to match direction points

	VERTEX3D Dir1, Dir2;
	Sub(&Dir1, GetVertex(Aligner_dest), GetVertex(Anchor_dest));
	Sub(&Dir2, GetVertex(Aligner_source), GetVertex(Anchor_source));
	Cross(&Axis, &Dir2, &Dir1);
	if (Norme(&Axis)<1e-5) { //The two directions are either collinear or the opposite
		if (Dot(&Dir1, &Dir2)>0.99999) { //no rotation needed
			Axis.x = 1.0;
			Axis.y = 0.0;
			Axis.z = 0.0;
			angle = 0.0;

		}
		else { //180deg rotation needed
			//construct a vector perpendicular to the normal
			Axis = facets[Facet_source]->sh.N;
			Normalize(&Axis);
			angle = 180.0;
		}

	}
	else {
		Normalize(&Axis);
		angle = Dot(&Dir1, &Dir2) / Norme(&Dir1) / Norme(&Dir2);
		//BOOL opposite=(angle<0.0);
		angle = acos(angle);
		angle = angle / PI * 180;
		//if (invertNormal) angle=180.0-angle;
	}

	BOOL *alreadyRotated2 = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyRotated2, FALSE, sh.nbVertex*sizeof(BOOL*));


	nb = 0;
	for (int i = 0; i < nbSelected; i++) {
		counter += 0.333;
		prgAlign->SetProgress(counter / selected);
		for (int j = 0; j < facets[selection[i]]->sh.nbIndex; j++) {
			if (!alreadyRotated2[facets[selection[i]]->indices[j]]) {
				//rotation comes here
				Rotate(&(vertices3[facets[selection[i]]->indices[j]]), *(GetVertex(Anchor_dest)), Axis, angle);
				alreadyRotated2[facets[selection[i]]->indices[j]] = TRUE;
			}
		}
	}

	SAFE_FREE(alreadyRotated2);

	InitializeGeometry();
	//update textures
	/*try {
		for (int i = 0; i < nbSelected; i++)
			SetFacetTexture(selection[i], facets[selection[i]]->tRatio, facets[selection[i]]->hasMesh);
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}*/
	prgAlign->SetVisible(FALSE);
	SAFE_DELETE(prgAlign);
}

void Geometry::MoveSelectedFacets(double dX, double dY, double dZ, BOOL copy, Worker *worker) {

	GLProgress *prgMove = new GLProgress("Moving selected facets...", "Please wait");
	prgMove->SetProgress(0.0);
	prgMove->SetVisible(TRUE);
	if (!(dX == 0.0&&dY == 0.0&&dZ == 0.0)) {
		if (!mApp->AskToReset(worker)) return;
		int nbSelFacet = 0;
		if (copy) CloneSelectedFacets(); //move
		double counter = 1.0;
		double selected = (double)GetNbSelected();
		if (selected == 0.0) return;

		BOOL *alreadyMoved = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
		memset(alreadyMoved, FALSE, sh.nbVertex*sizeof(BOOL*));


		int nb = 0;
		for (int i = 0; i < sh.nbFacet; i++) {
			if (facets[i]->selected) {
				counter += 1.0;
				prgMove->SetProgress(counter / selected);
				for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
					if (!alreadyMoved[facets[i]->indices[j]]) {
						vertices3[facets[i]->indices[j]].x += dX;
						vertices3[facets[i]->indices[j]].y += dY;
						vertices3[facets[i]->indices[j]].z += dZ;
						alreadyMoved[facets[i]->indices[j]] = TRUE;
					}
				}
			}
		}

		SAFE_FREE(alreadyMoved);

		InitializeGeometry();
		//update textures
		/*try {
			for (int i = 0; i < sh.nbFacet; i++) if (facets[i]->selected) SetFacetTexture(i, facets[i]->tRatio, facets[i]->hasMesh);
		}
		catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}*/
	}
	prgMove->SetVisible(FALSE);
	SAFE_DELETE(prgMove);
}

void Geometry::MirrorSelectedFacets(VERTEX3D P0, VERTEX3D N, BOOL copy, Worker *worker) {

	double selected = (double)GetNbSelected();
	double counter = 0.0;
	if (selected == 0.0) return;
	GLProgress *prgMirror = new GLProgress("Mirroring selected facets...", "Please wait");
	prgMirror->SetProgress(0.0);
	prgMirror->SetVisible(TRUE);

	if (!mApp->AskToReset(worker)) return;
	int nbSelFacet = 0;
	if (copy) CloneSelectedFacets(); //move
	BOOL *alreadyMirrored = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyMirrored, FALSE, sh.nbVertex*sizeof(BOOL*));


	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			counter += 1.0;
			prgMirror->SetProgress(counter / selected);
			nbSelFacet++;
			for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
				if (!alreadyMirrored[facets[i]->indices[j]]) {
					//mirroring comes here
					Mirror(&(vertices3[facets[i]->indices[j]]), P0, N);
					alreadyMirrored[facets[i]->indices[j]] = TRUE;
				}
			}
		}
	}

	SAFE_FREE(alreadyMirrored);
	if (nbSelFacet == 0) return;
	SwapNormal();
	InitializeGeometry();
	//update textures
	/*try {
		for (int i = 0; i < sh.nbFacet; i++) if (facets[i]->selected) SetFacetTexture(i, facets[i]->tRatio, facets[i]->hasMesh);
	}
	catch (Error &e) {
		GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
		return;
	}*/

	prgMirror->SetVisible(FALSE);
	SAFE_DELETE(prgMirror);
}

void Geometry::RotateSelectedFacets(VERTEX3D AXIS_P0, VERTEX3D AXIS_DIR, double theta, BOOL copy, Worker *worker) {

	double selected = (double)GetNbSelected();
	double counter = 0.0;
	if (selected == 0.0) return;
	GLProgress *prgRotate = new GLProgress("Rotating selected facets...", "Please wait");
	prgRotate->SetProgress(0.0);
	prgRotate->SetVisible(TRUE);

	if (theta != 0.0) {
		if (!mApp->AskToReset(worker)) return;
		if (copy) CloneSelectedFacets(); //move
		BOOL *alreadyRotated = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
		memset(alreadyRotated, FALSE, sh.nbVertex*sizeof(BOOL*));


		int nb = 0;
		for (int i = 0; i < sh.nbFacet; i++) {
			if (facets[i]->selected) {
				counter += 1.0;
				prgRotate->SetProgress(counter / selected);
				for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
					if (!alreadyRotated[facets[i]->indices[j]]) {
						//rotation comes here
						Rotate(&(vertices3[facets[i]->indices[j]]), AXIS_P0, AXIS_DIR, theta);
						alreadyRotated[facets[i]->indices[j]] = TRUE;
					}
				}
			}
		}

		SAFE_FREE(alreadyRotated);
		InitializeGeometry();
		//update textures
		/*try {
			for (int i = 0; i < sh.nbFacet; i++) if (facets[i]->selected) SetFacetTexture(i, facets[i]->tRatio, facets[i]->hasMesh);

		}
		catch (Error &e) {
			GLMessageBox::Display((char *)e.GetMsg(), "Error", GLDLG_OK, GLDLG_ICONERROR);
			return;
		}*/
	}
	prgRotate->SetVisible(FALSE);
	SAFE_DELETE(prgRotate);
}

void Geometry::CloneSelectedFacets() { //create clone of selected facets
	double counter = 0.0;
	double selected = (double)GetNbSelected();
	if (selected == 0.0) return;
	int *copyId = (int*)malloc(sh.nbVertex*sizeof(int*));
	memset(copyId, -1, sh.nbVertex*sizeof(int*));

	//count how many new vertices to create

	int nb = sh.nbVertex - 1;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			counter += 1.0;
			for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
				if (copyId[facets[i]->indices[j]] == -1) {
					nb++;
					copyId[facets[i]->indices[j]] = nb;
				}
			}
		}
	}
	/*
	VERTEX3D *newVertices=(VERTEX3D*)malloc((nb+1)*sizeof(VERTEX3D*));
	memcpy(newVertices,vertices3,sh.nbVertex*sizeof(VERTEX3D));
	SAFE_FREE(vertices3);
	*/
	//vertices3=(VERTEX3D *)realloc(vertices3,(nb+1)*sizeof(VERTEX3D));
	VERTEX3D *tmp_vertices3 = (VERTEX3D *)malloc((nb + 1) * sizeof(VERTEX3D)); //create new, extended vertex array
	memmove(tmp_vertices3, vertices3, (sh.nbVertex)*sizeof(VERTEX3D)); //copy old vertices
	memset(tmp_vertices3 + sh.nbVertex, 0, (nb + 1 - sh.nbVertex) * sizeof(VERTEX3D));  //zero out remaining bits (not necessary, will be overwritten anyway)

	SAFE_FREE(vertices3); //delete old array
	vertices3 = tmp_vertices3; //make new array the official vertex holder
	for (int i = 0; i < sh.nbVertex; i++) {
		if (copyId[i] != -1) {
			vertices3[copyId[i]].x = vertices3[i].x;
			vertices3[copyId[i]].y = vertices3[i].y;
			vertices3[copyId[i]].z = vertices3[i].z;
			vertices3[copyId[i]].selected = vertices3[i].selected;
		}
	}
	sh.nbVertex = nb + 1; //update number of vertices
	sh.nbFacet += (int)selected;
	facets = (Facet **)realloc(facets, sh.nbFacet * sizeof(Facet *));
	int nb2 = sh.nbFacet - (int)selected - 1; //copy new facets
	for (int i = 0; i < sh.nbFacet - (int)selected; i++) {
		if (facets[i]->selected) {
			nb2++;
			facets[nb2] = new Facet(facets[i]->sh.nbIndex);
			facets[nb2]->Copy(facets[i], FALSE);
			//copy indices
			for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
				facets[nb2]->indices[j] = facets[i]->indices[j];
			}
			facets[i]->selected = FALSE;
		}
	}
	for (int i = (sh.nbFacet - (int)selected); i < sh.nbFacet; i++) {
		for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
			facets[i]->indices[j] = copyId[facets[i]->indices[j]];
		}
	}
	CalcTotalOutGassing();

	//Debug memory check
	//_ASSERTE (!_CrtDumpMemoryLeaks());;
	_ASSERTE(_CrtCheckMemory());
}

void Geometry::MoveSelectedVertex(double dX, double dY, double dZ, BOOL copy, Worker *worker) {


	if (!(dX == 0.0&&dY == 0.0&&dZ == 0.0)) {
		if (!mApp->AskToReset(worker)) return;
		mApp->changedSinceSave = TRUE;
		if (!copy) { //move
			for (int i = 0; i < sh.nbVertex; i++) {
				if (vertices3[i].selected) {
					vertices3[i].x += dX;
					vertices3[i].y += dY;
					vertices3[i].z += dZ;
				}
			}
			InitializeGeometry();



		}
		else { //copy
			int nbVertexOri = sh.nbVertex;
			for (int i = 0; i < nbVertexOri; i++) {
				if (vertices3[i].selected) {
					AddVertex(vertices3[i].x + dX, vertices3[i].y + dY, vertices3[i].z + dZ);
				}
			}
		}

	}
}

void Geometry::AddVertex(double X, double Y, double Z) {

	mApp->changedSinceSave = TRUE;

	//a new vertex
	sh.nbVertex++;
	VERTEX3D *verticesNew = (VERTEX3D *)malloc(sh.nbVertex * sizeof(VERTEX3D));
	memcpy(verticesNew, vertices3, (sh.nbVertex - 1) * sizeof(VERTEX3D)); //copy old vertices
	SAFE_FREE(vertices3);
	verticesNew[sh.nbVertex - 1].x = X;
	verticesNew[sh.nbVertex - 1].y = Y;
	verticesNew[sh.nbVertex - 1].z = Z;
	verticesNew[sh.nbVertex - 1].selected = TRUE;
	vertices3 = verticesNew;

	//InitializeGeometry();

}

void Geometry::GetSelection(int **selection, int *nbSel){
	int sel = 0;
	*nbSel = GetNbSelected();
	int *selected = (int *)malloc((*nbSel)*sizeof(int));
	for (int i = 0; i < sh.nbFacet; i++)
		if (facets[i]->selected) selected[sel++] = i;
	*selection = selected;
}

void Geometry::SetSelection(int **selection, int *nbSel){

	UnSelectAll();
	for (int i = 0; i < *nbSel; i++) {
		int toSelect = (*selection)[i];
		if (toSelect < sh.nbFacet) facets[toSelect]->selected = TRUE;
	}
	UpdateSelection();
	if (*nbSel>0) mApp->facetList->ScrollToVisible((*selection)[*nbSel - 1], 0, TRUE); //in facet list, select the last facet of selection group
	mApp->UpdateFacetParams(TRUE);
}

void Geometry::AddStruct(char *name) {
	strName[sh.nbSuper++] = _strdup(name);
	BuildGLList();
}

void Geometry::DelStruct(int numToDel) {

	RemoveFromStruct(numToDel);
	CheckIsolatedVertex();
	mApp->UpdateModelParams();

	for (int i = 0; i<sh.nbFacet; i++) {
		if (facets[i]->sh.superIdx>numToDel) facets[i]->sh.superIdx--;
		if (facets[i]->sh.superDest > numToDel) facets[i]->sh.superDest--;
	}
	for (int j = numToDel; j < (sh.nbSuper - 1); j++)
	{
		strName[j] = _strdup(strName[j + 1]);
	}
	sh.nbSuper--;
	BuildGLList();
}

void Geometry::ScaleSelectedVertices(VERTEX3D invariant, double factor, BOOL copy, Worker *worker) {



	if (!mApp->AskToReset(worker)) return;
	mApp->changedSinceSave = TRUE;
	if (!copy) { //scale
		for (int i = 0; i < sh.nbVertex; i++) {
			if (vertices3[i].selected) {
				vertices3[i].x = invariant.x + factor*(vertices3[i].x - invariant.x);
				vertices3[i].y = invariant.y + factor*(vertices3[i].y - invariant.y);
				vertices3[i].z = invariant.z + factor*(vertices3[i].z - invariant.z);
			}
		}

	}
	else { //scale and copy
		int nbVertexOri = sh.nbVertex;
		for (int i = 0; i < nbVertexOri; i++) {
			if (vertices3[i].selected) {
				AddVertex(invariant.x + factor*(vertices3[i].x - invariant.x),
					invariant.y + factor*(vertices3[i].y - invariant.y),
					invariant.z + factor*(vertices3[i].z - invariant.z));
			}
		}
	}
	InitializeGeometry();
}

void Geometry::ScaleSelectedFacets(VERTEX3D invariant, double factorX, double factorY, double factorZ, BOOL copy, Worker *worker) {


	GLProgress *prgMove = new GLProgress("Scaling selected facets...", "Please wait");
	prgMove->SetProgress(0.0);
	prgMove->SetVisible(TRUE);

	if (!mApp->AskToReset(worker)) return;
	int nbSelFacet = 0;
	if (copy) CloneSelectedFacets(); //move
	double counter = 1.0;
	double selected = (double)GetNbSelected();
	if (selected == 0.0) return;

	BOOL *alreadyMoved = (BOOL*)malloc(sh.nbVertex*sizeof(BOOL*));
	memset(alreadyMoved, FALSE, sh.nbVertex*sizeof(BOOL*));


	int nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			counter += 1.0;
			prgMove->SetProgress(counter / selected);
			for (int j = 0; j < facets[i]->sh.nbIndex; j++) {
				if (!alreadyMoved[facets[i]->indices[j]]) {
					vertices3[facets[i]->indices[j]].x = invariant.x + factorX*(vertices3[facets[i]->indices[j]].x - invariant.x);
					vertices3[facets[i]->indices[j]].y = invariant.y + factorY*(vertices3[facets[i]->indices[j]].y - invariant.y);
					vertices3[facets[i]->indices[j]].z = invariant.z + factorZ*(vertices3[facets[i]->indices[j]].z - invariant.z);
					alreadyMoved[facets[i]->indices[j]] = TRUE;
				}
			}
		}
	}

	SAFE_FREE(alreadyMoved);

	InitializeGeometry();
	//update textures
	//for(int i=0;i<sh.nbFacet;i++) if(facets[i]->selected) SetFacetTexture(i,facets[i]->tRatio,facets[i]->hasMesh);	   

	prgMove->SetVisible(FALSE);
	SAFE_DELETE(prgMove);
}


void Geometry::SplitSelectedFacets(VERTEX3D base, VERTEX3D normal, Worker *worker) {
	/*
	mApp->changedSinceSave = TRUE;
	if (nbSelected == 0) return -1;

	// Check that all facet has a mesh
	BOOL ok = TRUE;
	int idx = 0;
	while (ok && idx < sh.nbFacet) {
		if (facets[idx]->selected)
			ok = facets[idx]->hasMesh;
		idx++;
	}
	if (!ok) return -2;

	int nb = 0;
	int FtoAdd = 0;
	int VtoAdd = 0;
	Facet::FACETGROUP *blocks = (Facet::FACETGROUP *)malloc(nbSelected * sizeof(Facet::FACETGROUP));

	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			facets[i]->Explode(blocks + nb);
			FtoAdd += blocks[nb].nbF;
			VtoAdd += blocks[nb].nbV;
			nb++;
		}
	}

	// Update vertex array
	VERTEX3D *ptrVert;
	int       vIdx;
	VERTEX3D *nVert = (VERTEX3D *)malloc((sh.nbVertex + VtoAdd)*sizeof(VERTEX3D));
	memcpy(nVert, vertices3, sh.nbVertex*sizeof(VERTEX3D));

	ptrVert = nVert + sh.nbVertex;
	vIdx = sh.nbVertex;
	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			facets[i]->FillVertexArray(ptrVert);
			for (int j = 0; j < blocks[nb].nbF; j++) {
				for (int k = 0; k < blocks[nb].facets[j]->sh.nbIndex; k++) {
					blocks[nb].facets[j]->indices[k] = vIdx + k;
				}
				vIdx += blocks[nb].facets[j]->sh.nbIndex;
			}
			ptrVert += blocks[nb].nbV;
			nb++;
		}
	}
	SAFE_FREE(vertices3);
	vertices3 = nVert;
	for (int i = sh.nbVertex; i < sh.nbVertex + VtoAdd; i++)
		vertices3[i].selected = FALSE;
	sh.nbVertex += VtoAdd;

	// Update facet
	Facet   **f = (Facet **)malloc((sh.nbFacet + FtoAdd - nbSelected) * sizeof(Facet *));

	// Delete selected
	nb = 0;
	for (int i = 0; i < sh.nbFacet; i++) {
		if (facets[i]->selected) {
			delete facets[i];
			mApp->RenumberSelections(i);
			mApp->RenumberFormulas(i);
		}
		else {
			f[nb++] = facets[i];
		}
	}

	// Add new facets
	int count = 0;
	for (int i = 0; i < nbSelected; i++) {
		for (int j = 0; j<blocks[i].nbF; j++) {
			f[nb++] = blocks[i].facets[j];
			if (toMap) { //set outgassing values
				f[nb - 1]->sh.flow = *(values + count++) *0.100; //0.1: mbar*l/s->Pa*m3/s
				if (f[nb - 1]->sh.flow>0.0) {
					f[nb - 1]->sh.desorbType = desType + 1;
					f[nb - 1]->selected = TRUE;
					if (f[nb - 1]->sh.desorbType == DES_COSINE_N) f[nb - 1]->sh.desorbTypeN = exponent;
				}
				else {
					f[nb - 1]->sh.desorbType = DES_NONE;
					f[nb - 1]->selected = FALSE;
				}
			}
		}
	}

	// Free allocated memory
	for (int i = 0; i < nbSelected; i++) {
		SAFE_FREE(blocks[i].facets);
	}
	SAFE_FREE(blocks);

	SAFE_FREE(facets);
	facets = f;
	sh.nbFacet = nb;

	// Delete old resources
	DeleteGLLists(TRUE, TRUE);

	InitializeGeometry();

	return 0;
	*/
}