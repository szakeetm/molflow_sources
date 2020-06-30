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
void serialize(Archive & archive,
               float2 & m)
{
    archive( m.x, m.y);
}

template<class Archive>
void serialize(Archive & archive,
               float3 & m)
{
    archive( m.x, m.y, m.z);
}

template<class Archive>
void serialize(Archive & archive,
               int3 & m)
{
    archive( m.x, m.y, m.z);
}

inline void vector3d_to_float3(float3& t, const Vector3d& o){
    t.x = o.x;
    t.y = o.y;
    t.z = o.z;

    return;
}


/*! \namespace flowgeom - Molflow Geometry code */
namespace flowgpu {

/*! --- Initialise model with a Molflow-exported geometry --- */
    flowgpu::Model* initializeModel(std::string fileName){

    std::cout << "#GPUTestsuite: Loading input file: " << fileName << std::endl;

    flowgpu::Model* model = new flowgpu::Model();
    std::ifstream file( fileName );
    cereal::XMLInputArchive archive( file );

#ifdef WITHTRIANGLES
    model->triangle_meshes.push_back(new flowgpu::TriangleMesh());
    archive(
            cereal::make_nvp("poly", model->triangle_meshes[0]->poly) ,
            cereal::make_nvp("facetProbabilities", model->triangle_meshes[0]->facetProbabilities) ,
            cereal::make_nvp("cdfs", model->triangle_meshes[0]->cdfs) ,
            //cereal::make_nvp("vertices2d", nullptr) ,
            cereal::make_nvp("vertices3d", model->triangle_meshes[0]->vertices3d) ,
            cereal::make_nvp("indices", model->triangle_meshes[0]->indices) ,
            cereal::make_nvp("nbFacets", model->triangle_meshes[0]->nbFacets) ,
            cereal::make_nvp("nbVertices", model->triangle_meshes[0]->nbVertices) ,
            cereal::make_nvp("nbFacetsTotal", model->nbFacets_total) ,
            cereal::make_nvp("nbVerticesTotal", model->nbVertices_total) ,
            cereal::make_nvp("useMaxwellDistribution", model->parametersGlobal.useMaxwellDistribution)
    );
#else
        model->poly_meshes.push_back(new flowgpu::PolygonMesh());
        archive(
                cereal::make_nvp("poly", model->poly_meshes[0]->poly) ,
                cereal::make_nvp("facetProbabilities", model->poly_meshes[0]->facetProbabilities) ,
                cereal::make_nvp("cdfs", model->poly_meshes[0]->cdfs) ,
                cereal::make_nvp("vertices2d", model->poly_meshes[0]->vertices2d) ,
                cereal::make_nvp("vertices3d", model->poly_meshes[0]->vertices3d) ,
                cereal::make_nvp("indices", model->poly_meshes[0]->indices) ,
                cereal::make_nvp("nbFacets", model->poly_meshes[0]->nbFacets) ,
                cereal::make_nvp("nbVertices", model->poly_meshes[0]->nbVertices) ,
                cereal::make_nvp("nbFacetsTotal", model->nbFacets_total) ,
                cereal::make_nvp("nbVerticesTotal", model->nbVertices_total) ,
                cereal::make_nvp("useMaxwellDistribution", model->parametersGlobal.useMaxwellDistribution)
        );
#endif

        std::cout << "#GPUTestsuite: Loading completed!" << std::endl;

        return model;
    }

    void convertFacet2Poly(const std::vector<TempFacet>& facets, std::vector<flowgpu::Polygon>& convertedPolygons){

        int32_t vertCount = 0;
        for(int i = 0; i < facets.size(); ++i){
            auto& temp = facets[i];

            flowgpu::Polygon polygon(temp.vertices2.size());
            polygon.facProps.stickingFactor = temp.facetProperties.sticking;
            polygon.facProps.temperature = temp.facetProperties.temperature;
            polygon.facProps.opacity = temp.facetProperties.opacity;

            polygon.facProps.is2sided = temp.facetProperties.is2sided;

            if(polygon.nbVertices != temp.facetProperties.nbIndex){
                polygon.nbVertices = temp.facetProperties.nbIndex;
                std::cout << "Parsing error! Vert size != nbIndex"<<std::endl;
                exit(0);
            }

            polygon.indexOffset = vertCount;
            vector3d_to_float3(polygon.O,temp.facetProperties.O);
            vector3d_to_float3(polygon.U,temp.facetProperties.U);
            vector3d_to_float3(polygon.V,temp.facetProperties.V);
            vector3d_to_float3(polygon.Nuv,temp.facetProperties.Nuv);
            vector3d_to_float3(polygon.nU,temp.facetProperties.nU);
            vector3d_to_float3(polygon.nV,temp.facetProperties.nV);
            vector3d_to_float3(polygon.N,temp.facetProperties.N);
            /*polygon.O = temp.facetProperties.O;
            polygon.U = temp.facetProperties.U;
            polygon.V = temp.facetProperties.V;
            polygon.Nuv = temp.facetProperties.Nuv;
            polygon.nU = temp.facetProperties.nU;
            polygon.nV = temp.facetProperties.nV;
            polygon.N = temp.facetProperties.N;*/

            polygon.parentIndex = i;

            convertedPolygons.push_back(std::move(polygon));

            vertCount += polygon.nbVertices;
        }

    }

    //! Calculate outgassing values in relation to (tri_area / poly_area)
    void CalculateRelativeTriangleOutgassing(const std::vector<TempFacet>& facets, flowgpu::TriangleMesh* triMesh){
        float fullOutgassing = 0;
        int facetIndex = 0;
        for(auto& facet : triMesh->poly){

            // Calculate triangle area
            auto& triIndices = triMesh->indices[facetIndex];
            auto& a = triMesh->vertices3d[triIndices.x];
            auto& b = triMesh->vertices3d[triIndices.y];
            auto& c = triMesh->vertices3d[triIndices.z];


            float3 ab = make_float3((b.x - a.x) , (b.y - a.y) , (b.z - a.z));
            float3 ac = make_float3((c.x - a.x) , (c.y - a.y) , (c.z - a.z));
            float area = 0.5 * length(cross(ab,ac));

            // outgassing of a triangle is only a percentage of the original polygon's
            float areaPercentageOfPoly = (area/facets[facet.parentIndex].facetProperties.area);

            float fullOutgassing_inc = fullOutgassing + (facets[facet.parentIndex].facetProperties.outgassing * areaPercentageOfPoly)/ (1.38E-23*facets[facet.parentIndex].facetProperties.temperature);
            triMesh->facetProbabilities.push_back(make_float2(fullOutgassing, fullOutgassing_inc));
            fullOutgassing = fullOutgassing_inc;

            ++facetIndex;
        }
        for(auto& facetProb : triMesh->facetProbabilities){
            facetProb.x /= fullOutgassing; // normalize to [0,1]
            facetProb.y /= fullOutgassing; // normalize to [0,1]
        }
    };

    int InitializeTexture(TempFacet& facet,
            std::vector<Texel>& texture,
            FacetTexture& facetTex,
            std::vector<float>& texInc)
    {

        //Textures
        if (facet.facetProperties.isTextured) {
            size_t nbE = facet.facetProperties.texWidth*facet.facetProperties.texHeight;
            facetTex.texHeight = facet.facetProperties.texHeight;
            facetTex.texWidth = facet.facetProperties.texWidth;
            facetTex.texHeightD = facet.facetProperties.texHeightD;
            facetTex.texWidthD = facet.facetProperties.texWidthD;

            if(facet.texelInc.size() != nbE){
                printf("texture inc vector has weird size: %d should be %d\n", facet.texelInc.size() , nbE);
                return 0;
            }

            try {
                texture = std::vector<Texel>(nbE);
                for(std::vector<double>::iterator it = facet.texelInc.begin();it != facet.texelInc.end();it++)
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

    int InitializeProfile(TempFacet &facet, std::vector<Texel> &profile)
    {

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

    int parseGeomFromSerialization(flowgpu::Model* model, std::vector<flowgpu::TempFacet>& facets, std::vector<float3>& vertices3d){
        // First create a regular polygonmesh
        // transform Molflow facet data to simulation polygons
        // transform polygons to triangles
        flowgpu::PolygonMesh* polyMesh = new flowgpu::PolygonMesh();
        convertFacet2Poly(facets,polyMesh->poly);

        int i= 0;
        int indexOffset = 0;

        for(auto& facet : facets){
            polyMesh->poly[i++].indexOffset = indexOffset;
            indexOffset += facet.indices.size();

            for(auto ind : facet.indices){
                polyMesh->indices.push_back(ind);
            }
            for(auto vert : facet.vertices2){
                polyMesh->vertices2d.push_back(make_float2(vert.u,vert.v));
            }
        }
        polyMesh->vertices3d = vertices3d;

        polyMesh->nbFacets = model->geomProperties.nbFacet;
        polyMesh->nbVertices = model->geomProperties.nbVertex;

        // Now create Triangle Mesh
        flowgpu::TriangleMesh* triMesh = new flowgpu::TriangleMesh();
        try {
            Poly2TriConverter::PolygonsToTriangles(polyMesh, triMesh);
        }
        catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
            return -1;
        }
        triMesh->vertices3d = polyMesh->vertices3d;
        triMesh->nbVertices = triMesh->poly.size() * 3;
        triMesh->nbFacets = triMesh->poly.size();

        triMesh->cdfs.push_back(0);

        if(!polyMesh->poly.empty())
            model->poly_meshes.push_back(polyMesh);
        model->triangle_meshes.push_back(triMesh);

        for(auto& polyMesh : model->poly_meshes){
            model->nbFacets_total += polyMesh->nbFacets;
            model->nbVertices_total += polyMesh->nbVertices;
        }
        for(auto& triMesh : model->triangle_meshes){
            model->nbFacets_total += triMesh->nbFacets;
            model->nbVertices_total += triMesh->nbVertices;
        }
        model->geomProperties.nbFacet = model->nbFacets_total;
        model->geomProperties.nbVertex = model->nbVertices_total;

        for(auto& facet : facets){
            model->tri_facetOffset.emplace_back(facet.facetProperties.hitOffset);
        }

        //--- Calculate outgassing values in relation to (tri_area / poly_area)
        CalculateRelativeTriangleOutgassing(facets,triMesh);

        for(auto& mesh : model->poly_meshes){
            std::vector<flowgpu::FacetType> sbtIndices; // Facet Type
            for(auto& polygon : mesh->poly){
                if(polygon.facProps.is2sided && polygon.facProps.opacity == 0.0f){
                    //std::cout << sbtIndices.size() << " > is transparent " << std::endl;
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_TRANS);
                }
                else{
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_SOLID);
                }
            }
            mesh->sbtIndices = sbtIndices;
        }
        for(auto& mesh : model->triangle_meshes){
            std::vector<flowgpu::FacetType> sbtIndices; // Facet Type
            for(auto& triangle : mesh->poly){
                if(triangle.facProps.is2sided && triangle.facProps.opacity == 0.0f){
                    std::cout << sbtIndices.size() << " > is transparent " << std::endl;
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_TRANS);
                }
                else{
                    sbtIndices.emplace_back(FacetType::FACET_TYPE_SOLID);
                }
            }
            mesh->sbtIndices = sbtIndices;
        }

        // Textures
        if(!model->textures.empty()){
            std::cout << "[WARNING] Textures would get added to non-empty vector!"<< std::endl;
            return 1;
        }
        else{
            int textureOffset = 0;
            int texelOffset = 0;
            for(int facetInd = 0; facetInd < facets.size(); ++facetInd){
                auto& facet = facets[facetInd];
                if(facet.facetProperties.isTextured){

                    std::vector<Texel> texture;
                    std::vector<float> texInc;
                    FacetTexture facetTex; // TODO: Map Offset to Triangle
                    if(!InitializeTexture(facet, texture, facetTex, texInc)){
                        printf("[ERROR] Initializing facetTex #%d\n", facetInd);
                        exit(0);
                    }

                    facetTex.bbMin = make_float3(std::numeric_limits<float>::max());
                    facetTex.bbMax = make_float3(std::numeric_limits<float>::min());

                    for(auto& ind : facet.indices){
                        facetTex.bbMin = fminf(facetTex.bbMin, vertices3d[ind]);
                        facetTex.bbMax = fmaxf(facetTex.bbMax, vertices3d[ind]);
                    }

                    facetTex.texelOffset = texelOffset;
                    texelOffset += texture.size();

                    model->textures.insert(std::end(model->textures),std::begin(texture),std::end(texture));
                    model->texInc.insert(std::end(model->texInc),std::begin(texInc),std::end(texInc));
                    model->facetTex.push_back(facetTex);

                    for(auto& polyMesh : model->poly_meshes){
                        for(auto& polygon : polyMesh->poly){
                            if(polygon.parentIndex == facetInd){
                                polygon.texProps.textureOffset = textureOffset;
                                polygon.texProps.textureSize = texture.size();
                                if(facet.facetProperties.countAbs)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countAbs;
                                if(facet.facetProperties.countRefl)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countRefl;
                                if(facet.facetProperties.countTrans)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countTrans;
                                if(facet.facetProperties.countDirection)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countDirection;
                                if(facet.facetProperties.countDes)
                                    polygon.texProps.textureFlags |= TEXTURE_FLAGS::countDes;
                            }
                        }
                    }
                    for(auto& triMesh : model->triangle_meshes){
                        for(auto& triangle : triMesh->poly){
                            if(triangle.parentIndex == facetInd){
                                triangle.texProps.textureOffset = textureOffset;
                                triangle.texProps.textureSize = texture.size();
                                if(facet.facetProperties.countAbs)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countAbs;
                                if(facet.facetProperties.countRefl)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countRefl;
                                if(facet.facetProperties.countTrans)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countTrans;
                                if(facet.facetProperties.countDirection)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countDirection;
                                if(facet.facetProperties.countDes)
                                    triangle.texProps.textureFlags |= TEXTURE_FLAGS::countDes;
                            }
                        }
                    }
                    textureOffset += 1;
                } // one more texture
            }
        }

        // Profiles
        if(!model->profiles.empty()){
            std::cout << "[WARNING] Profiles would get added to non-empty vector!"<< std::endl;
            return 1;
        }
        else{
            int currentProfOffset = 0;
            for(int facetInd = 0; facetInd < facets.size(); ++facetInd){
                auto& facet = facets[facetInd];
                if(facet.facetProperties.isProfile){

                    std::vector<Texel> profile;
                    if(!InitializeProfile(facet, profile)){
                        printf("[ERROR] Initializing facetProf #%d\n", facetInd);
                        exit(0);
                    }

                    // append to a global continuous structure
                    model->profiles.insert(std::end(model->profiles), std::begin(profile), std::end(profile));

                    for(auto& polygonMesh : model->poly_meshes){
                        for(auto& polygon : polygonMesh->poly){
                            if(polygon.parentIndex == facetInd){
                                polygon.profProps.profileOffset = currentProfOffset;
                                polygon.profProps.profileType = facet.facetProperties.profileType;
                            }
                        }
                    }
                    for(auto& triangleMesh : model->triangle_meshes){
                        for(auto& triangle : triangleMesh->poly){
                            if(triangle.parentIndex == facetInd){
                                triangle.profProps.profileOffset = currentProfOffset;
                                triangle.profProps.profileType = facet.facetProperties.profileType;
                            }
                        }
                    }
                    currentProfOffset += PROFILE_SIZE;
                } // one more profile
            }
        }

        for(auto& mesh : model->triangle_meshes){
            uint32_t triInd = 0;
            for(auto& triangle : mesh->poly){
                if(triangle.facProps.stickingFactor > 0){
                    std::cout << triangle.parentIndex << " | " << triInd << " => ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].x].x << " , " << mesh->vertices3d[mesh->indices[triInd].x].y << " , " << mesh->vertices3d[mesh->indices[triInd].x].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].y].x << " , " << mesh->vertices3d[mesh->indices[triInd].y].y << " , " << mesh->vertices3d[mesh->indices[triInd].y].z << " | ";
                    std::cout << mesh->vertices3d[mesh->indices[triInd].z].x << " , " << mesh->vertices3d[mesh->indices[triInd].z].y << " , " << mesh->vertices3d[mesh->indices[triInd].z].z << " | ";
                    std::cout << std::endl;
                }
                triInd++;
            }
        }

        return 0;
    }

    //! Load simulation data (geometry etc.) from Molflow's serialization output
    flowgpu::Model* loadFromSerialization(std::string inputString){
        std::stringstream inputStream;
        inputStream << inputString;
        cereal::BinaryInputArchive inputArchive(inputStream);

        flowgpu::Model* model = new flowgpu::Model();
        std::vector<Vector3d> vertices3d;
        std::vector<TempFacet> facets;

        //Worker params
        inputArchive(cereal::make_nvp("wp",model->wp));
        inputArchive(cereal::make_nvp("ontheflyParams",model->ontheflyParams));
        inputArchive(cereal::make_nvp("CDFs",model->CDFs));
        inputArchive(cereal::make_nvp("IDs",model->IDs));
        inputArchive(cereal::make_nvp("parameters",model->parameters));
        //inputArchive(cereal::make_nvp("temperatures",model->temperatures));
        inputArchive(cereal::make_nvp("moments",model->moments));
        //inputArchive(cereal::make_nvp("desorptionParameterIDs",model->desorptionParameterIDs));

        //Geometry
        inputArchive(cereal::make_nvp("GeomProperties",model->geomProperties));
        inputArchive(vertices3d);

        model->structures.resize(model->geomProperties.nbSuper); //Create structures

        //Facets
        facets.resize(model->geomProperties.nbFacet);

        for(int facInd = 0; facInd < model->geomProperties.nbFacet;++facInd){
            inputArchive(cereal::make_nvp("facet"+std::to_string(facInd),facets[facInd]));
        }

        std::vector<float3> vertices3f;
        for(std::vector<Vector3d>::iterator it=vertices3d.begin(); it!=vertices3d.end();it++)
            vertices3f.emplace_back(make_float3(it->x,it->y,it->z));
        vertices3d.clear();
        if(parseGeomFromSerialization(model, facets, vertices3f)){
            delete model;
            model = nullptr;
            return model;
        }
        std::cout << "#ModelReader: Gas mass: " << model->wp.gasMass << std::endl;
        std::cout << "#ModelReader: Maxwell: " << model->wp.useMaxwellDistribution << std::endl;
        std::cout << "#ModelReader: Name: " << model->geomProperties.name << std::endl;
        std::cout << "#ModelReader: #Vertex: " << vertices3d.size() << std::endl;
        std::cout << "#ModelReader: #Facets: " << model->geomProperties.nbFacet << std::endl;
        for(int facInd = 0; facInd < facets.size(); ++facInd){
            if(!facets[facInd].texelInc.empty())
                std::cout << "#ModelReader: Facet#" << facInd << " #Texels: " << facets[facInd].texelInc.size() << std::endl;
            if(facets[facInd].facetProperties.isTextured)
                std::cout << "#ModelReader: Facet#" << facInd << " #Texels: " << facets[facInd].texelInc.size() << std::endl;
        }
        for(auto& mesh : model->triangle_meshes){
            std::cout << "#ModelReader: #Tri " << mesh->poly.size()<< std:: endl;
            size_t nbProf = 0;
            for(auto& tri : mesh->poly){
                if(tri.profProps.profileType != PROFILE_FLAGS::noProfile)
                    ++nbProf;
            }
            std::cout << "#ModelReader: #ProfiledTris " << nbProf << std::endl;
        }

        std::cout << "#ModelReader: #TextureCells: " << model->textures.size() << std::endl;

        return model;
    }

    //! Load simulation data (geometry etc.) from Molflow's serialization file output
    flowgpu::Model* loadFromExternalSerialization(std::string fileName){
        std::ifstream file( fileName );
        cereal::XMLInputArchive inputarchive(file);

        flowgpu::Model* model = new flowgpu::Model();
        std::vector<float3> vertices3d;
        std::vector<TempFacet> facets;

        //Worker params
        //inputarchive(cereal::make_nvp("wp",model->parametersGlobal));
        inputarchive(cereal::make_nvp("gasMass",model->parametersGlobal.gasMass));
        inputarchive(cereal::make_nvp("useMaxwellDistribution",model->parametersGlobal.useMaxwellDistribution));
        inputarchive(cereal::make_nvp("GeomProperties",model->geomProperties));
        inputarchive(vertices3d);
        facets.resize(model->geomProperties.nbFacet);
        for(int facInd = 0; facInd < model->geomProperties.nbFacet;++facInd){
            inputarchive(cereal::make_nvp("facet"+std::to_string(facInd),facets[facInd]));
        }

        std::cout << "#ModelReader: Gas mass: " << model->parametersGlobal.gasMass << std::endl;
        std::cout << "#ModelReader: Maxwell: " << model->parametersGlobal.useMaxwellDistribution << std::endl;
        std::cout << "#ModelReader: Name: " << model->geomProperties.name << std::endl;
        std::cout << "#ModelReader: #Vertex: " << vertices3d.size() << std::endl;
        std::cout << "#ModelReader: #Facets: " << model->geomProperties.nbFacet << std::endl;
        for(int facInd = 0; facInd < facets.size(); ++facInd){
            if(!facets[facInd].texelInc.empty())
                std::cout << "#ModelReader: Facet#" << facInd << " #Texels: " << facets[facInd].texelInc.size() << std::endl;
        }

        parseGeomFromSerialization(model, facets, vertices3d);
        std::cout << "#ModelReader: #TextureCells: " << model->textures.size() << std::endl;

        return model;
    }

}// namespace