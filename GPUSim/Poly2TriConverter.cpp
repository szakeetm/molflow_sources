//
// Created by pbahr on 13/02/2020.
//

#include "Poly2TriConverter.h"

#define DET22(_11,_12,_21,_22) ( (_11)*(_22) - (_21)*(_12) )
#define DET33(_11,_12,_13,_21,_22,_23,_31,_32,_33)  \
  ((_11)*( (_22)*(_33) - (_32)*(_23) ) +            \
   (_12)*( (_23)*(_31) - (_33)*(_21) ) +            \
   (_13)*( (_21)*(_32) - (_31)*(_22) ))


//TODO:: Make use of MathTools.h/.cpp
size_t  IDX(int i, size_t nb) {
    return (i < 0) ? (nb + i) : (i%nb);
}

size_t IDX(size_t i, size_t nb) {
    return i%nb;
}

size_t Next(int i, size_t nb) {
    return (i+1)%nb;
}

size_t Next(size_t i, size_t nb) {
    return (i+1)%nb;
}

size_t Previous(int i, size_t nb) {
    return IDX(i - 1, nb);
}

size_t Previous(size_t i, size_t nb) {
    return IDX((int)i - 1, nb);
}

bool IsConvex(const std::vector<float2>& points,size_t idx) {

    // Check if p.pts[idx] is a convex vertex (calculate the sign of the oriented angle)

    size_t i1 = Previous(idx,points.size());
    size_t i2 = IDX(idx, points.size());
    size_t i3 = Next(idx, points.size());

    double d = DET22(points[i1].x - points[i2].x,points[i3].x - points[i2].x,
                     points[i1].y - points[i2].y,points[i3].y - points[i2].y);

    //return (d*p.sign)>=0.0;
    return d <= 0.0;
}

bool IsInPoly(const float2 &p, const std::vector<float2>& polyPoints) {

    // Fast method to check if a point is inside a polygon or not.
    // Works with convex and concave polys, orientation independent
    int n_updown = 0;
    int n_found = 0;

    size_t nbSizeMinusOne = polyPoints.size() - 1;
    for (size_t j = 0; j < nbSizeMinusOne; j++) {
        const float2& p1 = polyPoints[j];
        const float2& p2 = polyPoints[j+1];

        /*
        if (p1.x < p2.x) {
             minx = p1.x;
             maxx = p2.x;
        }
        else {
             minx = p2.x;
             maxx = p1.x;
        }
        if (p.x > minx && p.x <= maxx) {
        */

        if (p.x<p1.x != p.x<p2.x) {
            double slope = (p2.y - p1.y) / (p2.x - p1.x);
            if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
                n_updown++;
            }
            else {
                n_updown--;
            }
            n_found++;
        }
    }
    //Last point. Repeating code because it's the fastest and this function is heavily used
    const float2& p1 = polyPoints[nbSizeMinusOne];
    const float2& p2 = polyPoints[0];

    if (p.x<p1.x != p.x<p2.x) {
        double slope = (p2.y - p1.y) / (p2.x - p1.x);
        if ((slope * p.x - p.y) < (slope * p1.x - p1.y)) {
            n_updown++;
        }
        else {
            n_updown--;
        }
        n_found++;
    }

    return ((n_found / 2) & 1) ^ ((n_updown / 2) & 1);
}

bool ContainsConcave(const std::vector<float2>& points,int i1,int i2,int i3)
{
    // Determine if the specified triangle contains or not a concave point
    size_t _i1 = IDX(i1, points.size());
    size_t _i2 = IDX(i2, points.size());
    size_t _i3 = IDX(i3, points.size());

    const float2& p1 = points[_i1];
    const float2& p2 = points[_i2];
    const float2& p3 = points[_i3];

    int found = 0;
    int i = 0;
    while(!found && i<points.size()) {
        if( i!=_i1 && i!=_i2 && i!=_i3 ) {
            if (IsInPoly(points[i], { p1,p2,p3 }))
                found = !IsConvex(points,i);
        }
        i++;
    }

    return found;
}

int  FindEar(const std::vector<float2>& points){

    int i = 0;
    bool earFound = false;
    while (i < points.size() && !earFound) {
        if (IsConvex(points, i))
            earFound = !ContainsConcave(points, i - 1, i, i + 1);
        if (!earFound) i++;
    }

    // REM: Theoritically, it should always find an ear (2-Ears theorem).
    // However on degenerated geometry (flat poly) it may not find one.
    // Returns first point in case of failure.
    if (earFound)
        return i;
    else
        return 0;
}

// Return indices to create a new triangle
int3 GetTriangleFromEar(const flowgeom::TempFacet& facet, const std::vector<float2>& points, const int ear) {

    int3 indices;
    indices.x = facet.indices[Previous(ear, points.size())];
    indices.y = facet.indices[IDX(ear, points.size())];
    indices.z = facet.indices[Next(ear, points.size())];

    return indices;
}

std::vector<int3> Poly2TriConverter::Triangulate(flowgeom::TempFacet& facet) {

    // Triangulate a facet (rendering purpose)
    // The facet must have at least 3 points
    // Use the very simple "Two-Ears" theorem. It computes in O(n^2).
    std::vector<int3> triangles;
    std::vector<float2>& vertices = facet.vertices2;
    std::vector<uint32_t>& indices = facet.indices;

    // Perform triangulation
    while (vertices.size() > 3) {
        int e = FindEar(vertices);
        //DrawEar(f, p, e, addTextureCoord);

        // Create new triangle facet and copy polygon parameters, but change indices
        int3 triangleIndices = GetTriangleFromEar(facet, vertices, e);
        triangles.push_back(triangleIndices);

        // Remove the ear
        vertices.erase(vertices.begin() + e);
        indices.erase(indices.begin() + e);

    }

    // Draw the last ear
    int3 triangleIndices = GetTriangleFromEar(facet, vertices, 0);
    triangles.push_back(triangleIndices);

    return triangles;
}

// Update facet list of geometry by removing polygon facets and replacing them with triangular facets with the same properties
int Poly2TriConverter::PolygonsToTriangles(flowgpu::PolygonMesh *polygonMesh, flowgpu::TriangleMesh *triangleMesh) {

    std::vector<flowgeom::Polygon> convertedTris;
    std::vector<flowgeom::Polygon>& polygons = polygonMesh->poly;
    for (size_t facetIndex = 0; facetIndex < polygons.size(); facetIndex++) {
        size_t nbVert = polygons[facetIndex].nbVertices;

        if (nbVert == 3){
            convertedTris.push_back(polygons[facetIndex]);
            triangleMesh->indices.push_back(make_int3(polygonMesh->indices[polygons[facetIndex].indexOffset],
                                             polygonMesh->indices[polygons[facetIndex].indexOffset+1],
                                             polygonMesh->indices[polygons[facetIndex].indexOffset+2]));

        }
        else if (nbVert > 3) {

            // Prepare new temp facet to pass to the triangulation algorithm
            flowgeom::TempFacet temp;
            for(int i = 0; i < nbVert; ++i){
                temp.indices.push_back(polygonMesh->indices[polygons[facetIndex].indexOffset + i]);
                temp.vertices2.push_back(polygonMesh->vertices2d[polygons[facetIndex].indexOffset + i]);
            }

            if(nbVert != temp.vertices2.size())
                std::cout << "[WARNING] Polygon with "<<temp.vertices2.size()<<" vertices should have "<< nbVert <<std::endl;

            std::vector<int3> triangleIndices = Triangulate(temp);
            triangleMesh->indices.insert(std::end(triangleMesh->indices),std::begin(triangleIndices),std::end(triangleIndices));
            std::vector<flowgeom::Polygon> newTris;
            for(auto& tri : triangleIndices){

                // TODO: This should be checked before triangulation on the polygon itself
                if(tri.x == tri.y || tri.x == tri.z || tri.y == tri.z){
                    std::cout << "[WARNING] Triangle with same vertices was created! PolyIndex: "<< facetIndex <<std::endl;
                    std::cout << "[WARNING] Vertices: "<< tri.x << " , " << tri.y << " , " << tri.z << std::endl;
                    continue;
                }
                flowgeom::Polygon newPoly(3);
                newPoly.parentIndex = facetIndex;
                newPoly.indexOffset = polygons[facetIndex].indexOffset;
                newTris.push_back(newPoly);
            }
            convertedTris.insert(std::end(convertedTris),std::begin(newTris),std::end(newTris));

            if(newTris.size() > nbVert - 2)
                std::cout << "[WARNING] Polygon with "<<nbVert<<" vertices was split into "<< newTris.size() << " triangles!" <<std::endl;
            //delete geometry->facets[facetIndex];
        }
        else{
            //triangleFacets.push_back(geometry->facets[facetIndex]);
        }
    }

    // Lookup parent IDs in the original polygon mesh
    for(auto& tri : convertedTris){
        for(std::vector<flowgeom::Polygon>::iterator polyIter = polygons.begin(); polyIter != polygons.end(); ++polyIter){
            if(tri.parentIndex == (*polyIter).parentIndex){
                tri.copyParametersFrom(*polyIter);

                break;
            }
        }
    }

    // Second lookup to delete
    for(std::vector<flowgeom::Polygon>::iterator polyIter = polygons.begin(); polyIter != polygons.end(); ++polyIter){
        if((*polyIter).nbVertices == 3){
            //std::cout << "Deleting Tri# "<<(*polyIter).parentIndex<< " from PolyList"<<std::endl;
            polygons.erase(polyIter);
            break;
        }
    }
    for(auto& tri : convertedTris){
        for(std::vector<flowgeom::Polygon>::iterator polyIter = polygons.begin(); polyIter != polygons.end(); ++polyIter){
            if(tri.parentIndex == (*polyIter).parentIndex){
                //std::cout << "Deleting Poly# "<<tri.parentIndex<<std::endl;
                polygons.erase(polyIter);
                break;
            }
        }
    }

    std::cout << "Amount of n>3 Polygons after triangulation: "<<polygons.size()<<std::endl;
    std::cout << "Amount of Triangles after triangulation: "<<convertedTris.size()<<std::endl;


    /*geometry->sh.nbFacet = triangleFacets.size();
    geometry->facets = (Facet **)realloc(geometry->facets, geometry->sh.nbFacet * sizeof(Facet *));

    // Update facet list
    int i = 0;
    for(Facet* facet : triangleFacets){
        geometry->facets[i++] = facet;
    }

    // to recalculate various facet properties
    geometry->InitializeGeometry();*/

    triangleMesh->poly.insert(std::end(triangleMesh->poly),std::begin(convertedTris),std::end(convertedTris));
    return convertedTris.size();
}

// Update facet list of geometry by removing polygon facets and replacing them with triangular facets with the same properties
//TODO: Parameter should be the Model (with vertices etc.)
std::vector<flowgeom::Polygon> Poly2TriConverter::PolygonsToTriangles(std::vector<flowgeom::TempFacet>& facets){

    std::vector<flowgeom::Polygon> convertedTris;
    for (size_t facetIndex = 0; facetIndex < facets.size(); facetIndex++) {
        size_t nb = facets[facetIndex].indices.size();
        if (nb > 3) {
            // Create new triangle facets and invalidate old polygon facet
            std::vector<int3> triangleIndices = Triangulate(facets[facetIndex]);
            std::vector<flowgeom::Polygon> newTris;
            for(auto& tri : triangleIndices){
                flowgeom::Polygon newPoly(3);
                newPoly.parentIndex = facetIndex;
                newTris.push_back(newPoly);
            }
            convertedTris.insert(std::end(convertedTris),std::begin(newTris),std::end(newTris));
            //delete geometry->facets[facetIndex];
        }
        else{
            //triangleFacets.push_back(geometry->facets[facetIndex]);
        }
    }

    /*geometry->sh.nbFacet = triangleFacets.size();
    geometry->facets = (Facet **)realloc(geometry->facets, geometry->sh.nbFacet * sizeof(Facet *));

    // Update facet list
    int i = 0;
    for(Facet* facet : triangleFacets){
        geometry->facets[i++] = facet;
    }

    // to recalculate various facet properties
    geometry->InitializeGeometry();*/

    return convertedTris;
}