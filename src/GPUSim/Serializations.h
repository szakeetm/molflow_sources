//
// Created by pbahr on 1/12/22.
//

#ifndef MOLFLOW_PROJ_SERIALIZATIONS_H
#define MOLFLOW_PROJ_SERIALIZATIONS_H

#include "cereal/cereal.hpp"
#include "OptixPolygon.h"

template<class Archive>
void serialize(Archive & archive,
               flowgpu::TempFacetProperties & prop) {
    archive(
            CEREAL_NVP(
                    prop.sticking),       // Sticking (0=>reflection  , 1=>absorption)   - can be overridden by time-dependent parameter
            CEREAL_NVP(prop.opacity),        // opacity  (0=>transparent , 1=>opaque)
            CEREAL_NVP(prop.area),          // Facet area (m^2)

            CEREAL_NVP(prop.profileType),    // Profile type
            CEREAL_NVP(prop.superIdx),       // Super structure index (Indexed from 0)
            CEREAL_NVP(prop.superDest),      // Super structure destination index (Indexed from 1, 0=>current)
            CEREAL_NVP(
                    prop.teleportDest),   // Teleport destination facet id (for periodic boundary condition) (Indexed from 1, 0=>none, -1=>teleport to where it came from)

            CEREAL_NVP(prop.countAbs),       // Count absoprtion (MC texture)
            CEREAL_NVP(prop.countRefl),      // Count reflection (MC texture)
            CEREAL_NVP(prop.countTrans),     // Count transparent (MC texture)
            CEREAL_NVP(prop.countDirection),

            // Flags
            CEREAL_NVP(prop.is2sided),     // 2 sided
            CEREAL_NVP(prop.isProfile),    // Profile facet
            CEREAL_NVP(prop.isTextured),   // texture
            CEREAL_NVP(
                    prop.isVolatile),   // Volatile facet (absorbtion facet which does not affect particule trajectory)

            // Geometry
            CEREAL_NVP(prop.nbIndex),   // Number of index/vertex
            //CEREAL_NVP(prop.sign),      // Facet vertex rotation (see Facet::DetectOrientation())

            // Plane basis (O,U,V) (See Geometry::InitializeGeometry() for info)
            CEREAL_NVP(prop.O),  // Origin
            CEREAL_NVP(prop.U),  // U vector
            CEREAL_NVP(prop.V),  // V vector
            CEREAL_NVP(prop.nU), // Normalized U
            CEREAL_NVP(prop.nV), // Normalized V

            // Normal vector
            CEREAL_NVP(prop.N),    // normalized
            CEREAL_NVP(prop.Nuv),  // normal to (u,v) not normlized

            // Hit/Abs/Des/Density recording on 2D texture map
            CEREAL_NVP(prop.texWidth),    // Rounded texture resolution (U)
            CEREAL_NVP(prop.texHeight),   // Rounded texture resolution (V)
            CEREAL_NVP(prop.texWidthD),   // Actual texture resolution (U)
            CEREAL_NVP(prop.texHeightD),  // Actual texture resolution (V)


            // Molflow-specific facet parameters
            CEREAL_NVP(
                    prop.temperature),    // Facet temperature (Kelvin)                  - can be overridden by time-dependent parameter
            CEREAL_NVP(
                    prop.outgassing),           // (in unit *m^3/s)                      - can be overridden by time-dependent parameter

            CEREAL_NVP(prop.desorbType),     // Desorption type
            CEREAL_NVP(prop.desorbTypeN),    // Exponent in Cos^N desorption type

            CEREAL_NVP(prop.countDes),       // Count desoprtion (MC texture)

            CEREAL_NVP(prop.totalOutgassing) //total outgassing for the given facet
    );
}

template<class Archive>
void serialize(Archive & archive,
               flowgpu::Polygon & poly)
{
    archive(
            poly.nbVertices,
            poly.indexOffset,
            poly.O,
            poly.U,
            poly.V,
            poly.Nuv,
            poly.nU,
            poly.nV,
            poly.N
    );
}

#endif //MOLFLOW_PROJ_SERIALIZATIONS_H
