//
// Created by pbahr on 13/02/2020.
//

#ifndef MOLFLOW_PROJ_POLY2TRICONVERTER_H
#define MOLFLOW_PROJ_POLY2TRICONVERTER_H

#include "ModelReader.h"

class Poly2TriConverter {
    static std::vector<int3> Triangulate(flowgeom::TempFacet& facet);
public:
    static std::vector<flowgeom::Polygon> PolygonsToTriangles(std::vector<flowgeom::TempFacet>& facets);
    static int
    PolygonsToTriangles(flowgpu::PolygonMesh *polygonMesh, flowgpu::TriangleMesh *triangleMesh);
};


#endif //MOLFLOW_PROJ_POLY2TRICONVERTER_H
