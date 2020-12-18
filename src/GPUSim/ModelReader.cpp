//
// Created by pbahr on 12/02/2020.
//

#include "ModelReader.h"
#include "Poly2TriConverter.h"
//#include "Facet_shared.h"

// debug output
#include <fstream>

#include "helper_math.h"
#include "../MolflowTypes.h"

#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/tuple.hpp>

/*#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>*/

template<class Archive>
void serialize(Archive &archive,
               float2 &m) {
    archive(m.x, m.y);
}

template<class Archive>
void serialize(Archive &archive,
               float3 &m) {
    archive(m.x, m.y, m.z);
}

template<class Archive>
void serialize(Archive &archive,
               int3 &m) {
    archive(m.x, m.y, m.z);
}

inline void vector3d_to_float3(float3 &t, const Vector3d &o) {
    t.x = o.x;
    t.y = o.y;
    t.z = o.z;

    return;
}

inline float3 vector3d_to_float3(const Vector3d &o) {
    float3 t;
    t.x = o.x;
    t.y = o.y;
    t.z = o.z;

    return t;
}

inline void vector3d_to_double3(double3 &t, const Vector3d &o) {
    t.x = o.x;
    t.y = o.y;
    t.z = o.z;

    return;
}

inline double3 vector3d_to_double3(const Vector3d &o) {
    double3 t;
    t.x = o.x;
    t.y = o.y;
    t.z = o.z;

    return t;
}
/*! \namespace flowgeom - Molflow Geometry code */
namespace flowgpu {

/*! --- Initialise model with a Molflow-exported geometry --- */
    flowgpu::Model *initializeModel(std::string fileName) {

        std::cout << "#GPUTestsuite: Loading input file: " << fileName << std::endl;

        flowgpu::Model *model = new flowgpu::Model();
        std::ifstream file(fileName);
        cereal::XMLInputArchive archive(file);

#ifdef WITHTRIANGLES
        model->triangle_meshes.push_back(new flowgpu::TriangleMesh());
        archive(
                cereal::make_nvp("poly", model->triangle_meshes[0]->poly),
                cereal::make_nvp("facetProbabilities", model->triangle_meshes[0]->facetProbabilities),
                cereal::make_nvp("cdfs", model->triangle_meshes[0]->cdfs),
                //cereal::make_nvp("vertices2d", nullptr) ,
                cereal::make_nvp("vertices3d", model->triangle_meshes[0]->vertices3d),
                cereal::make_nvp("indices", model->triangle_meshes[0]->indices),
                cereal::make_nvp("nbFacets", model->triangle_meshes[0]->nbFacets),
                cereal::make_nvp("nbVertices", model->triangle_meshes[0]->nbVertices),
                cereal::make_nvp("nbFacetsTotal", model->nbFacets_total),
                cereal::make_nvp("nbVerticesTotal", model->nbVertices_total),
                cereal::make_nvp("useMaxwellDistribution", model->parametersGlobal.useMaxwellDistribution)
        );
#else
        model->poly_meshes.push_back(new flowgpu::PolygonMesh());
        archive(
                cereal::make_nvp("poly", model->poly_meshes[0]->poly),
                cereal::make_nvp("facetProbabilities", model->poly_meshes[0]->facetProbabilities),
                cereal::make_nvp("cdfs", model->poly_meshes[0]->cdfs),
                cereal::make_nvp("vertices2d", model->poly_meshes[0]->vertices2d),
                cereal::make_nvp("vertices3d", model->poly_meshes[0]->vertices3d),
                cereal::make_nvp("indices", model->poly_meshes[0]->indices),
                cereal::make_nvp("nbFacets", model->poly_meshes[0]->nbFacets),
                cereal::make_nvp("nbVertices", model->poly_meshes[0]->nbVertices),
                cereal::make_nvp("nbFacetsTotal", model->nbFacets_total),
                cereal::make_nvp("nbVerticesTotal", model->nbVertices_total),
                cereal::make_nvp("useMaxwellDistribution", model->parametersGlobal.useMaxwellDistribution)
        );
#endif

        std::cout << "#GPUTestsuite: Loading completed!" << std::endl;

        return model;
    }

    void convertFacet2Poly(const std::vector<TempFacet> &facets, std::vector<flowgpu::Polygon> &convertedPolygons) {

        int32_t vertCount = 0;
        std::vector<size_t> deleteFacets; //mark some facets for deleting e.g. outgassing on linked

        for (int i = 0; i < facets.size(); ++i) {
            auto &temp = facets[i];

            flowgpu::Polygon polygon(temp.vertices2.size());
            polygon.facProps.stickingFactor = temp.facetProperties.sticking;
            polygon.facProps.temperature = temp.facetProperties.temperature;

            //TODO: For now, treat 'link facets' like transparent facets
            if (temp.facetProperties.superDest) {
                polygon.facProps.opacity = 0.0f;
            } else {
                polygon.facProps.opacity = temp.facetProperties.opacity;
            }
            polygon.facProps.is2sided = temp.facetProperties.is2sided;

            // desorption
            polygon.desProps.desorbType = temp.facetProperties.desorbType;
            polygon.desProps.cosineExponent = temp.facetProperties.desorbTypeN;

            if (polygon.nbVertices != temp.facetProperties.nbIndex) {
                polygon.nbVertices = temp.facetProperties.nbIndex;
                std::cout << "Parsing error! Vert size != nbIndex" << std::endl;
                exit(0);
            }

            polygon.indexOffset = vertCount;
            vector3d_to_float3(polygon.O, temp.facetProperties.O);
            vector3d_to_float3(polygon.U, temp.facetProperties.U);
            vector3d_to_float3(polygon.V, temp.facetProperties.V);
            vector3d_to_float3(polygon.Nuv, temp.facetProperties.Nuv);
            vector3d_to_float3(polygon.nU, temp.facetProperties.nU);
            vector3d_to_float3(polygon.nV, temp.facetProperties.nV);
            vector3d_to_float3(polygon.N, temp.facetProperties.N);

            vector3d_to_double3(polygon.Ox64, temp.facetProperties.O);
            vector3d_to_double3(polygon.Ux64, temp.facetProperties.U);
            vector3d_to_double3(polygon.Vx64, temp.facetProperties.V);
            vector3d_to_double3(polygon.Nuvx64, temp.facetProperties.Nuv);
            vector3d_to_double3(polygon.nUx64, temp.facetProperties.U);
            vector3d_to_double3(polygon.nVx64, temp.facetProperties.V);
            vector3d_to_double3(polygon.Nx64, temp.facetProperties.N);

            /*polygon.O = temp.facetProperties.O;
            polygon.U = temp.facetProperties.U;
            polygon.V = temp.facetProperties.V;
            polygon.Nuv = temp.facetProperties.Nuv;
            polygon.nU = temp.facetProperties.nU;
            polygon.nV = temp.facetProperties.nV;
            polygon.N = temp.facetProperties.N;*/

            polygon.parentIndex = i;
            vertCount += polygon.nbVertices;

            convertedPolygons.push_back(std::move(polygon));

        }

    }

    //! Calculate outgassing values in relation to (tri_area / poly_area)
    void CalculateRelativePolygonOutgassing(const std::vector<TempFacet> &facets, flowgpu::PolygonMesh *polyMesh) {
        float fullOutgassing = 0;
        int facetIndex = 0;
        for (auto &facet : polyMesh->poly) {

            // outgassing of a triangle is only a percentage of the original polygon's

            float fullOutgassing_inc = fullOutgassing +
                                       (facets[facet.parentIndex].facetProperties.outgassing) /
                                       (1.38E-23 * facets[facet.parentIndex].facetProperties.temperature);
            polyMesh->facetProbabilities.push_back(make_float2(fullOutgassing, fullOutgassing_inc));
            fullOutgassing = fullOutgassing_inc;

            ++facetIndex;
        }
        for (auto &facetProb : polyMesh->facetProbabilities) {
            facetProb.x /= fullOutgassing; // normalize to [0,1]
            facetProb.y /= fullOutgassing; // normalize to [0,1]
        }
    };

    //! Calculate outgassing values in relation to (tri_area / poly_area)
    void CalculateRelativeTriangleOutgassing(const std::vector<TempFacet> &facets, flowgpu::TriangleMesh *triMesh) {
        float fullOutgassing = 0;
        int facetIndex = 0;
        for (auto &facet : triMesh->poly) {

            // Calculate triangle area
            auto &triIndices = triMesh->indices[facetIndex];
            auto &a = triMesh->vertices3d[triIndices.x];
            auto &b = triMesh->vertices3d[triIndices.y];
            auto &c = triMesh->vertices3d[triIndices.z];


            float3 ab = make_float3((b.x - a.x), (b.y - a.y), (b.z - a.z));
            float3 ac = make_float3((c.x - a.x), (c.y - a.y), (c.z - a.z));
            float area = 0.5 * length(cross(ab, ac));

            // outgassing of a triangle is only a percentage of the original polygon's
            float areaPercentageOfPoly = (area / facets[facet.parentIndex].facetProperties.area);

            float fullOutgassing_inc = fullOutgassing +
                                       (facets[facet.parentIndex].facetProperties.outgassing * areaPercentageOfPoly) /
                                       (1.38E-23 * facets[facet.parentIndex].facetProperties.temperature);
            triMesh->facetProbabilities.push_back(make_float2(fullOutgassing, fullOutgassing_inc));
            fullOutgassing = fullOutgassing_inc;

            ++facetIndex;
        }
        for (auto &facetProb : triMesh->facetProbabilities) {
            facetProb.x /= fullOutgassing; // normalize to [0,1]
            facetProb.y /= fullOutgassing; // normalize to [0,1]
        }
    };

    int InitializeTexture(TempFacet &facet,
                          std::vector<Texel> &texture,
                          FacetTexture &facetTex,
                          std::vector<float> &texInc) {

        //Textures
        if (facet.facetProperties.isTextured) {
            size_t nbE = facet.facetProperties.texWidth * facet.facetProperties.texHeight;
            facetTex.texHeight = facet.facetProperties.texHeight;
            facetTex.texWidth = facet.facetProperties.texWidth;
            facetTex.texHeightD = facet.facetProperties.texHeightD;
            facetTex.texWidthD = facet.facetProperties.texWidthD;

            // Add a cutoff to allow for texture positions [0.0,1.0], opposed to [0.0,1.0[
            const float cutOff = 0.9999999f;
            facetTex.texHeightD *= cutOff;
            facetTex.texWidthD *= cutOff;

            if (facet.texelInc.size() != nbE) {
                printf("texture inc vector has weird size: %d should be %d\n", facet.texelInc.size(), nbE);
                return 0;
            }

            try {
                texture = std::vector<Texel>(nbE);
                for (auto it = facet.texelInc.begin(); it != facet.texelInc.end(); it++)
                    texInc.emplace_back(*it);
                //texInc = facet.texelInc;
            }
            catch (...) {
                printf("Not enough memory to load textures\n");
                return 0;
            }
        }
        return 1;
    }

    int InitializeProfile(TempFacet &facet, std::vector<Texel> &profile) {

        //Profile
        if (facet.facetProperties.isProfile) {
            //const size_t PROFILE_SIZE = 100;

            try {
                // PROFILE_SIZE bins per profile
                profile.resize(PROFILE_SIZE);
            }
            catch (...) {
                printf("Not enough memory to load profiles\n");
                return 0;
            }
        }
        return 1;
    }

    int parseGeomFromSerialization(flowgpu::Model *model, std::vector<flowgpu::TempFacet> &facets,
                                   std::vector<float3> &vertices3d) {
        // First create a regular polygonmesh
        // transform Molflow facet data to simulation polygons
        // transform polygons to triangles
        flowgpu::PolygonMesh *polyMesh = new flowgpu::PolygonMesh();
        convertFacet2Poly(facets, polyMesh->poly);

        int i = 0;
        int indexOffset = 0;

        for (auto &facet : facets) {
            polyMesh->poly[i++].indexOffset = indexOffset;
            indexOffset += facet.indices.size();

            for (auto ind : facet.indices) {
                polyMesh->indices.push_back(ind);
            }
            for (auto vert : facet.vertices2) {
                polyMesh->vertices2d.push_back(make_float2(vert.u, vert.v));
                polyMesh->vertices2d64.push_back(make_double2(vert.u, vert.v));

            }
        }
        polyMesh->vertices3d = vertices3d;

        polyMesh->nbFacets = model->geomProperties.nbFacet;
        polyMesh->nbVertices = model->geomProperties.nbVertex;

#ifdef WITHTRIANGLES
        // Now create Triangle Mesh
        flowgpu::TriangleMesh *triMesh = new flowgpu::TriangleMesh();
        try {
            Poly2TriConverter::PolygonsToTriangles(polyMesh, triMesh);
        }
        catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            return -1;
        }
        triMesh->vertices3d = polyMesh->vertices3d;
        triMesh->nbVertices = triMesh->poly.size() * 3;
        triMesh->nbFacets = triMesh->poly.size();

        triMesh->cdfs.push_back(0);

        model->triangle_meshes.push_back(triMesh);
#endif

        if (!polyMesh->poly.empty()) {
            polyMesh->cdfs.push_back(0);
            model->poly_meshes.push_back(polyMesh);
        }
        for (auto &polyMesh : model->poly_meshes) {
            model->nbFacets_total += polyMesh->nbFacets;
            model->nbVertices_total += polyMesh->nbVertices;
        }
        for (auto &triMesh : model->triangle_meshes) {
            model->nbFacets_total += triMesh->nbFacets;
            model->nbVertices_total += triMesh->nbVertices;
        }
        model->geomProperties.nbFacet = model->nbFacets_total;
        model->geomProperties.nbVertex = model->nbVertices_total;

        for (auto &facet : facets) {
            model->tri_facetOffset.emplace_back(facet.facetProperties.hitOffset);
        }

#ifdef WITHTRIANGLES
        //--- Calculate outgassing values in relation to (tri_area / poly_area)
        CalculateRelativeTriangleOutgassing(facets, triMesh);
#else
        CalculateRelativePolygonOutgassing(facets, polyMesh);
#endif
        for (auto &mesh : model->poly_meshes) {
            std::vector<flowgpu::FacetType> sbtIndices; // Facet Type
            for (auto &polygon : mesh->poly) {
                //if(polygon.facProps.is2sided && polygon.facProps.opacity == 0.0f){
                if (polygon.facProps.opacity == 0.0f) {
                    //std::cout << sbtIndices.size() << " > is transparent " << std::endl;
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_TRANS);
                } else {
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_SOLID);
                }
            }
            mesh->sbtIndices = sbtIndices;
        }
        for (auto &mesh : model->triangle_meshes) {
            std::vector<flowgpu::FacetType> sbtIndices; // Facet Type
            for (auto &triangle : mesh->poly) {
                //if(triangle.facProps.is2sided && triangle.facProps.opacity == 0.0f){
                if (triangle.facProps.opacity == 0.0f) {
                    std::cout << sbtIndices.size() << " > is transparent " << std::endl;
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_TRANS);
                } else {
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_SOLID);
                }
            }
            mesh->sbtIndices = sbtIndices;
        }

        // Textures
        if (!model->textures.empty()) {
            std::cout << "[WARNING] Textures would get added to non-empty vector!" << std::endl;
            return 1;
        } else {
            int textureOffset = 0;
            int texelOffset = 0;
            for (int facetInd = 0; facetInd < facets.size(); ++facetInd) {
                auto &facet = facets[facetInd];
                if (facet.facetProperties.isTextured) {

                    std::vector<Texel> texture;
                    std::vector<float> texInc;
                    FacetTexture facetTex; // TODO: Map Offset to Triangle
                    if (!InitializeTexture(facet, texture, facetTex, texInc)) {
                        printf("[ERROR] Initializing facetTex #%d\n", facetInd);
                        exit(0);
                    }

                    facetTex.bbMin = make_float3(std::numeric_limits<float>::max());
                    facetTex.bbMax = make_float3(std::numeric_limits<float>::lowest());

                    for (auto &ind : facet.indices) {
                        facetTex.bbMin = fminf(facetTex.bbMin, vertices3d[ind]);
                        facetTex.bbMax = fmaxf(facetTex.bbMax, vertices3d[ind]);
                    }

                    std::cout << "Texture BBox: " << facetTex.bbMin.x << " , " << facetTex.bbMin.y << " , "
                              << facetTex.bbMin.z << " | " <<
                              facetTex.bbMax.x << " , " << facetTex.bbMax.y << " , " << facetTex.bbMax.z << std::endl;


                    facetTex.texelOffset = texelOffset;
                    texelOffset += texture.size();

                    model->textures.insert(std::end(model->textures), std::begin(texture), std::end(texture));
                    model->texInc.insert(std::end(model->texInc), std::begin(texInc), std::end(texInc));
                    model->facetTex.push_back(facetTex);

                    for (auto &polyMesh : model->poly_meshes) {
                        for (auto &polygon : polyMesh->poly) {
                            if (polygon.parentIndex == facetInd) {
                                polygon.texProps.textureOffset = textureOffset;
                                polygon.texProps.textureSize = texture.size();
                                if (facet.facetProperties.countAbs)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countAbs;
                                if (facet.facetProperties.countRefl)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countRefl;
                                if (facet.facetProperties.countTrans)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countTrans;
                                if (facet.facetProperties.countDirection)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countDirection;
                                if (facet.facetProperties.countDes)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countDes;
                            }
                        }
                    }
                    for (auto &triMesh : model->triangle_meshes) {
                        uint32_t triInd = 0;
                        for (auto &triangle : triMesh->poly) {
                            if (triangle.parentIndex == facetInd) {
                                triangle.texProps.textureOffset = textureOffset;
                                triangle.texProps.textureSize = texture.size();
                                if (facet.facetProperties.countAbs)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countAbs;
                                if (facet.facetProperties.countRefl)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countRefl;
                                if (facet.facetProperties.countTrans)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countTrans;
                                if (facet.facetProperties.countDirection)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countDirection;
                                if (facet.facetProperties.countDes)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countDes;

                                if (triangle.texProps.textureFlags != noTexture) {
                                    // TODO: Map Texcoords

                                    /*std::cout << triangle.parentIndex << " | " << triInd << " => ";
                                    std::cout << " Origin " << triangle.O.x << " , " << triangle.O.y << " , " << triangle.O.z << " => ";
                                    std::cout << " U " << triangle.U.x << " , " << triangle.U.y << " , " << triangle.U.z << " => ";
                                    std::cout << " V " << triangle.V.x << " , " << triangle.V.y << " , " << triangle.V.z << " => ";
                                    //facetTex. = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
                                    float3 c1 = (triMesh->vertices3d[triMesh->indices[triInd].x] - facetTex.bbMin) / (facetTex.bbMax-facetTex.bbMin);
                                    float3 c2 = (triMesh->vertices3d[triMesh->indices[triInd].y] - facetTex.bbMin) / (facetTex.bbMax-facetTex.bbMin);
                                    float3 c3 = (triMesh->vertices3d[triMesh->indices[triInd].z] - facetTex.bbMin) / (facetTex.bbMax-facetTex.bbMin);

                                    std::cout << c1.x << " , " << c1.y << " , " << c1.z << " | ";
                                    std::cout << c2.x << " , " << c2.y << " , " << c2.z << " | ";
                                    std::cout << c3.x << " , " << c3.y << " , " << c3.z << " | ";
                                    std::cout << "---" << std::endl;*/
                                }

                            }

                            triInd++;
                        }
                    }
                    textureOffset += 1;
                } // one more texture

                // tex coords for all triangles
                {

                    for (auto &triMesh : model->triangle_meshes) {
                        uint32_t triInd = 0;
                        for (auto &triangle : triMesh->poly) {
                            if (triangle.parentIndex == facetInd) {
                                // TODO: Map Texcoords

                                /*std::cout << triangle.parentIndex << " | " << triInd << " => ";
                                std::cout << " Origin " << triangle.O.x << " , " << triangle.O.y << " , "
                                          << triangle.O.z << " => ";
                                std::cout << " U " << triangle.U.x << " , " << triangle.U.y << " , " << triangle.U.z
                                          << " => ";
                                std::cout << " V " << triangle.V.x << " , " << triangle.V.y << " , " << triangle.V.z
                                          << " => ";*/
                                //facetTex. = destination->sh.O + u * destination->sh.U + v * destination->sh.V;
                                /*float3 c1 = (triMesh->vertices3d[triMesh->indices[triInd].x] - facetTex.bbMin) /
                                            (facetTex.bbMax - facetTex.bbMin);
                                float3 c2 = (triMesh->vertices3d[triMesh->indices[triInd].y] - facetTex.bbMin) /
                                            (facetTex.bbMax - facetTex.bbMin);
                                float3 c3 = (triMesh->vertices3d[triMesh->indices[triInd].z] - facetTex.bbMin) /
                                            (facetTex.bbMax - facetTex.bbMin);*/
                                //float3 c4 = (triMesh->vertices3d[triMesh->indices[triInd].x] - triangle.O) / (facetTex.bbMax-facetTex.bbMin);
                                /*double u =  0.0;
                                double v = 0.0;
                                if(facet.facetProperties.V.x * facet.facetProperties.U.y != 0.0) {
                                    double det = facet.facetProperties.U.x * facet.facetProperties.V.y - facet.facetProperties.U.y * facet.facetProperties.V.x; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) - vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    double detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));
                                }
                                else if(facet.facetProperties.V.y * facet.facetProperties.U.z != 0.0) {
                                    double det = facet.facetProperties.U.y * facet.facetProperties.V.z - facet.facetProperties.U.z * facet.facetProperties.V.y; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) - vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    double detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));
                                }
                                else if(facet.facetProperties.V.z * facet.facetProperties.U.x != 0.0) {
                                    double det = facet.facetProperties.U.z * facet.facetProperties.V.x - facet.facetProperties.U.x * facet.facetProperties.V.z; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) - vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    double detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) - vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(make_float2(static_cast<float>(u),static_cast<float>(v)));
                                }*/

                                /*
                                std::cout << "\n";
                                std::cout << triMesh->vertices3d[triMesh->indices[triInd].x].x << " , " << triMesh->vertices3d[triMesh->indices[triInd].x].y << " , " << triMesh->vertices3d[triMesh->indices[triInd].x].z << " | ";
                                std::cout << triMesh->vertices3d[triMesh->indices[triInd].y].x << " , " << triMesh->vertices3d[triMesh->indices[triInd].y].y << " , " << triMesh->vertices3d[triMesh->indices[triInd].y].z << " | ";
                                std::cout << triMesh->vertices3d[triMesh->indices[triInd].z].x << " , " << triMesh->vertices3d[triMesh->indices[triInd].z].y << " , " << triMesh->vertices3d[triMesh->indices[triInd].z].z << " --- ";
        */
                                /*std::cout << c1.x << " , " << c1.y << " , " << c1.z << " | ";
                                std::cout << c2.x << " , " << c2.y << " , " << c2.z << " | ";
                                std::cout << c3.x << " , " << c3.y << " , " << c3.z << " | ";
                                std::cout << "---" << std::endl;*/



                                // texture coords for all triangles even without textures
                                double u = 0.0;
                                double v = 0.0;
                                if (facet.facetProperties.V.x * facet.facetProperties.U.y != 0.0) {
                                    double det = facet.facetProperties.U.x * facet.facetProperties.V.y -
                                                 facet.facetProperties.U.y *
                                                 facet.facetProperties.V.x; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) -
                                                vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    double detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;

                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.x * facet.facetProperties.V.y - b.y * facet.facetProperties.V.x;
                                    detV = facet.facetProperties.U.x * b.y - facet.facetProperties.U.y * b.x;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));
                                } else if (facet.facetProperties.V.y * facet.facetProperties.U.z != 0.0) {
                                    double det = facet.facetProperties.U.y * facet.facetProperties.V.z -
                                                 facet.facetProperties.U.z *
                                                 facet.facetProperties.V.y; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) -
                                                vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    double detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.y * facet.facetProperties.V.z - b.z * facet.facetProperties.V.y;
                                    detV = facet.facetProperties.U.y * b.z - facet.facetProperties.U.z * b.y;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));
                                } else if (facet.facetProperties.V.z * facet.facetProperties.U.x != 0.0) {
                                    double det = facet.facetProperties.U.z * facet.facetProperties.V.x -
                                                 facet.facetProperties.U.x *
                                                 facet.facetProperties.V.z; // TODO: Pre calculate

                                    // Vert 1
                                    double3 b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].x]) -
                                                vector3d_to_double3(facet.facetProperties.O);
                                    double detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    double detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 2
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].y]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));

                                    // Vert 3
                                    b = make_double3(triMesh->vertices3d[triMesh->indices[triInd].z]) -
                                        vector3d_to_double3(facet.facetProperties.O);
                                    detU = b.z * facet.facetProperties.V.x - b.x * facet.facetProperties.V.z;
                                    detV = facet.facetProperties.U.z * b.x - facet.facetProperties.U.x * b.z;
                                    v = (detV) / (det);
                                    u = (detU) / det;
                                    triMesh->texCoords.emplace_back(
                                            make_float2(static_cast<float>(u), static_cast<float>(v)));
                                }

                                /*float2 tex = *(triMesh->texCoords.rbegin() + 2);
                                std::cout << "["<<triMesh->texCoords.size()<<"] " << tex.x << " , " << tex.y << " | ";
                                tex = *(triMesh->texCoords.rbegin() + 1);
                                std::cout << tex.x << " , " << tex.y << " | ";
                                tex = *(triMesh->texCoords.rbegin());
                                std::cout << tex.x << " , " << tex.y << " | ";
                                std::cout << "---" << std::endl;*/
                            }
                            triInd++;
                        }
                    }
                }
            }
#ifdef BOUND_CHECK
            model->nbTexel_total = texelOffset;
#endif
        }

        // Profiles
        if (!model->profiles.empty()) {
            std::cout << "[WARNING] Profiles would get added to non-empty vector!" << std::endl;
            return 1;
        } else {
            int currentProfOffset = 0;
            for (int facetInd = 0; facetInd < facets.size(); ++facetInd) {
                auto &facet = facets[facetInd];
                if (facet.facetProperties.isProfile) {

                    std::vector<Texel> profile;
                    if (!InitializeProfile(facet, profile)) {
                        printf("[ERROR] Initializing facetProf #%d\n", facetInd);
                        exit(0);
                    }

                    // append to a global continuous structure
                    model->profiles.insert(std::end(model->profiles), std::begin(profile), std::end(profile));

                    for (auto &polygonMesh : model->poly_meshes) {
                        for (auto &polygon : polygonMesh->poly) {
                            if (polygon.parentIndex == facetInd) {
                                polygon.profProps.profileOffset = currentProfOffset;
                                polygon.profProps.profileType = facet.facetProperties.profileType;
                            }
                        }
                    }
                    for (auto &triangleMesh : model->triangle_meshes) {
                        for (auto &triangle : triangleMesh->poly) {
                            if (triangle.parentIndex == facetInd) {
                                triangle.profProps.profileOffset = currentProfOffset;
                                triangle.profProps.profileType = facet.facetProperties.profileType;
                            }
                        }
                    }
                    currentProfOffset += PROFILE_SIZE;
                } // one more profile
            }
#ifdef BOUND_CHECK
            model->nbProfSlices_total = currentProfOffset;
#endif
        }

        /*for(auto& mesh : model->triangle_meshes){
            uint32_t triInd = 0;
            for(auto& triangle : mesh->poly){
                *//*if(triangle.facProps.stickingFactor > 0){
                    std::cout << triangle.parentIndex << " | " << triInd << " => ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].x].x << " , " << mesh->vertices3d[mesh->indices[triInd].x].y << " , " << mesh->vertices3d[mesh->indices[triInd].x].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].y].x << " , " << mesh->vertices3d[mesh->indices[triInd].y].y << " , " << mesh->vertices3d[mesh->indices[triInd].y].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].z].x << " , " << mesh->vertices3d[mesh->indices[triInd].z].y << " , " << mesh->vertices3d[mesh->indices[triInd].z].z << " | ";
                    std::cout << std::endl;
                }*//*
                if(triangle.texProps.textureFlags!= noTexture){
                    std::cout << triangle.parentIndex << " | " << triInd << " => ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].x].x << " , " << mesh->vertices3d[mesh->indices[triInd].x].y << " , " << mesh->vertices3d[mesh->indices[triInd].x].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].y].x << " , " << mesh->vertices3d[mesh->indices[triInd].y].y << " , " << mesh->vertices3d[mesh->indices[triInd].y].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].z].x << " , " << mesh->vertices3d[mesh->indices[triInd].z].y << " , " << mesh->vertices3d[mesh->indices[triInd].z].z << " | ";
                    std::cout << std::endl;
                }
                triInd++;
            }
        }*/

        return 0;
    }

    //! Load simulation data (geometry etc.) from Molflow's serialization output
    flowgpu::Model *loadFromSerialization(const std::string &inputString) {
        std::stringstream inputStream;
        inputStream << inputString;
        cereal::BinaryInputArchive inputArchive(inputStream);

        flowgpu::Model *model = new flowgpu::Model();
        std::vector<Vector3d> vertices3d;
        std::vector<TempFacet> facets;

        //Worker params
        inputArchive(cereal::make_nvp("wp", model->wp));
        inputArchive(cereal::make_nvp("ontheflyParams", model->ontheflyParams));
        inputArchive(cereal::make_nvp("CDFs", model->CDFs));
        inputArchive(cereal::make_nvp("IDs", model->IDs));
        inputArchive(cereal::make_nvp("parameters", model->parameters));
        //inputArchive(cereal::make_nvp("temperatures",model->temperatures));
        inputArchive(cereal::make_nvp("moments", model->moments));
        //inputArchive(cereal::make_nvp("desorptionParameterIDs",model->desorptionParameterIDs));

        //Geometry
        inputArchive(cereal::make_nvp("GeomProperties", model->geomProperties));
        inputArchive(vertices3d);

        model->structures.resize(model->geomProperties.nbSuper); //Create structures

        //Facets
        facets.resize(model->geomProperties.nbFacet);

        for (int facInd = 0; facInd < model->geomProperties.nbFacet; ++facInd) {
            inputArchive(cereal::make_nvp("facet" + std::to_string(facInd), facets[facInd]));
        }

        std::vector<float3> vertices3f;
        for (std::vector<Vector3d>::iterator it = vertices3d.begin(); it != vertices3d.end(); it++)
            vertices3f.emplace_back(make_float3(it->x, it->y, it->z));
        vertices3d.clear();
        if (parseGeomFromSerialization(model, facets, vertices3f)) {
            delete model;
            model = nullptr;
            return model;
        }
        std::cout << "#ModelReader: Gas mass: " << model->wp.gasMass << std::endl;
        std::cout << "#ModelReader: Maxwell: " << model->wp.useMaxwellDistribution << std::endl;
        std::cout << "#ModelReader: Name: " << model->geomProperties.name << std::endl;
        std::cout << "#ModelReader: #Vertex: " << vertices3d.size() << std::endl;
        std::cout << "#ModelReader: #Facets: " << model->geomProperties.nbFacet << std::endl;
        size_t nbTexelCount = 0;
        for (int facInd = 0; facInd < facets.size(); ++facInd) {
            if ((facets[facInd].texelInc.empty() && facets[facInd].facetProperties.isTextured) ||
                (!facets[facInd].texelInc.empty() && !facets[facInd].facetProperties.isTextured)) {
                std::cerr << "#ModelReader: [ERROR] Texture flag and texel increments don't align! " << std::endl;
                _sleep(10000);
                exit(0);
            }
            if (facets[facInd].facetProperties.isTextured) {
                std::cout << "#ModelReader: Facet#" << facInd << " #Texels: " << facets[facInd].texelInc.size()
                          << std::endl;
                nbTexelCount += facets[facInd].texelInc.size();
            }
        }
        for (auto &mesh : model->triangle_meshes) {
            std::cout << "#ModelReader: #Tri " << mesh->poly.size() << std::endl;
            size_t nbProf = 0;
            size_t nbTex = 0;
            size_t triCount = 0;
            for (auto &tri : mesh->poly) {
                if (tri.profProps.profileType != PROFILE_FLAGS::noProfile)
                    ++nbProf;
                if (tri.texProps.textureFlags != TEXTURE_FLAGS::noTexture) {
                    ++nbTex;
                    std::cout << "#ModelReader: Tri#" << triCount++ << "[" << tri.parentIndex << "] #TexOffset: "
                              << tri.texProps.textureOffset << std::endl;
                }
            }
            std::cout << "#ModelReader: #ProfiledTris " << nbProf << std::endl;
            std::cout << "#ModelReader: #TexturedTris " << nbTex << std::endl;
        }

        std::cout << "#ModelReader: #TextureCells: " << model->textures.size() << std::endl;

        if (nbTexelCount != model->textures.size()) {
            std::cerr << "#ModelReader: [ERROR] Texture count out of sync: " << nbTexelCount << " / "
                      << model->textures.size() << std::endl;
            _sleep(10000);
            exit(0);
        }
        return model;
    }

    //! Load simulation data (geometry etc.) from Molflow's serialization file output
    flowgpu::Model *loadFromExternalSerialization(const std::string &fileName) {
        std::ifstream file(fileName);
        cereal::XMLInputArchive inputarchive(file);

        flowgpu::Model *model = new flowgpu::Model();
        std::vector<float3> vertices3d;
        std::vector<TempFacet> facets;

        //Worker params
        //inputarchive(cereal::make_nvp("wp",model->parametersGlobal));
        inputarchive(cereal::make_nvp("gasMass", model->parametersGlobal.gasMass));
        inputarchive(cereal::make_nvp("useMaxwellDistribution", model->parametersGlobal.useMaxwellDistribution));
        inputarchive(cereal::make_nvp("GeomProperties", model->geomProperties));
        inputarchive(vertices3d);
        facets.resize(model->geomProperties.nbFacet);
        for (int facInd = 0; facInd < model->geomProperties.nbFacet; ++facInd) {
            inputarchive(cereal::make_nvp("facet" + std::to_string(facInd), facets[facInd]));
        }

        std::cout << "#ModelReader: Gas mass: " << model->parametersGlobal.gasMass << std::endl;
        std::cout << "#ModelReader: Maxwell: " << model->parametersGlobal.useMaxwellDistribution << std::endl;
        std::cout << "#ModelReader: Name: " << model->geomProperties.name << std::endl;
        std::cout << "#ModelReader: #Vertex: " << vertices3d.size() << std::endl;
        std::cout << "#ModelReader: #Facets: " << model->geomProperties.nbFacet << std::endl;
        for (int facInd = 0; facInd < facets.size(); ++facInd) {
            if (!facets[facInd].texelInc.empty())
                std::cout << "#ModelReader: Facet#" << facInd << " #Texels: " << facets[facInd].texelInc.size()
                          << std::endl;
        }

        parseGeomFromSerialization(model, facets, vertices3d);
        std::cout << "#ModelReader: #TextureCells: " << model->textures.size() << std::endl;

        return model;
    }

}// namespace