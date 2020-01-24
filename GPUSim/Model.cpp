//
// Created by pbahr on 08/10/2019.
//

#include "Model.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "3rdParty/tiny_obj_loader.h"
//std
#include <set>
//#include <Facet_shared.h>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>

namespace std {
    inline bool operator<(const tinyobj::index_t &a,
                          const tinyobj::index_t &b)
    {
        if (a.vertex_index < b.vertex_index) return true;
        if (a.vertex_index > b.vertex_index) return false;

        if (a.normal_index < b.normal_index) return true;
        if (a.normal_index > b.normal_index) return false;

        if (a.texcoord_index < b.texcoord_index) return true;
        if (a.texcoord_index > b.texcoord_index) return false;

        return false;
    }
}

/*! \namespace flowgpu - Molflow GPU code */
namespace flowgpu {

    /*inline void vector3d_to_vec3f(gdt::vec3f& t, const Vector3d& o){
        t.x = o.x;
        t.y = o.y;
        t.z = o.z;

        return;
    }*/

    /*Model *loadFromMolflow(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs){
        Model *model = new Model;
        // Loop over Vertices3
        //thrust::host_vector<int> X(10);

        int polyNb = 0;
        int indNb = 0;
        int vertNb = 0;
        for(int i = 0; i < structures.size(); i++){
            polyNb += structures[i].facets.size();
            for(int f = 0; f < structures[i].facets.size(); f++){
                indNb += structures[i].facets[f].indices.size();
                vertNb += structures[i].facets[f].vertices2.size();
            }
        }

        for ( int s = 0; s < structures.size(); s++){
            //TODO: one mesh per structure

            PolygonMesh *mesh = new PolygonMesh;
            mesh->nbFacets = polyNb;
            mesh->nbIndices = indNb;
            mesh->nbVertices = indNb;

            std::vector<vec3f>(geomVertices.size()).swap(mesh->vertices3d);
            mesh->poly.resize(mesh->nbFacets);
            mesh->indices.resize(mesh->nbIndices);
            mesh->vertices2d.resize(mesh->nbVertices);

            int size_cdf = 0;
            for(auto& cdf : CDFs){
                size_cdf += 2*cdf.size(); // size for temp bins, v bins and both pair values
            }
            mesh->cdfs.resize(size_cdf);
            mesh->facetProbabilities.resize(mesh->nbFacets);

            int32_t vertCount = 0;
            double fullOutgassing = 0;

            for( int f = 0; f < structures[s].facets.size(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                const SubprocessFacet& fac = structures[s].facets[f];
                Polygon newPoly(fac.indices.size());
                newPoly.vertOffset = vertCount;
                newPoly.nbVertices = fac.indices.size();
                newPoly.stickingFactor = fac.sh.sticking;

                vector3d_to_vec3f(newPoly.O, fac.sh.O);
                vector3d_to_vec3f(newPoly.U, fac.sh.U);
                vector3d_to_vec3f(newPoly.V, fac.sh.V);
                vector3d_to_vec3f(newPoly.Nuv, fac.sh.Nuv);
                vector3d_to_vec3f(newPoly.nU, fac.sh.nU);
                vector3d_to_vec3f(newPoly.nV, fac.sh.nV);
                vector3d_to_vec3f(newPoly.N, fac.sh.N);

                int counter = 0;
                for(size_t vertInd : fac.indices){

                    // load vertex corresponding to index
                    vec3f newVert;
                    vector3d_to_vec3f(newVert, geomVertices[vertInd]);

                    //newVert *= 1.0f / dot(newVert,newVert);

                    mesh->vertices3d[vertInd] = newVert;
                    mesh->indices[vertCount] = vertInd;
                    mesh->vertices2d[vertCount] = vec2f(fac.vertices2[counter].u, fac.vertices2[counter].v);
                    counter++;
                    vertCount++;
                }
                mesh->poly[f] = std::move(newPoly);

                double facOutgassing = wp.latestMoment*fac.sh.outgassing / (1.38E-23*fac.sh.temperature);
                mesh->facetProbabilities[f].s = fullOutgassing; // lower border
                fullOutgassing += facOutgassing;
                mesh->facetProbabilities[f].t = fullOutgassing; // upper border
            }

            for( int f = 0; f < structures[s].facets.size(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                const SubprocessFacet& fac = structures[s].facets[f];
                mesh->facetProbabilities[f].s /= fullOutgassing; // normalize from [0,1]
                mesh->facetProbabilities[f].t /= fullOutgassing; // normalize from [0,1]
                //std::cout << "2. "<< fullOutgassing << " - ["<<mesh->facetProbabilities[f].s<<","<<mesh->facetProbabilities[f].t<<"]"<<std::endl;

            }

            int index = 0;
            for(auto& cdf : CDFs){
                for(auto& pair : cdf){
                    mesh->cdfs[index++] = pair.first; // %2==0
                    mesh->cdfs[index++] = pair.second; // %2==1
                }
            }

            // final sanity check
            if (mesh->vertices3d.empty() || mesh->indices.empty() || mesh->vertices2d.empty())
                delete mesh;
            else
                model->poly_meshes.push_back(mesh);
        }

        // calculate total amount of facets
        model->nbFacets_total = 0;
        for (PolygonMesh* mesh : model->poly_meshes)
            model->nbFacets_total += mesh->nbFacets;

        // calculate global bounds for the whole model
        for (PolygonMesh* mesh : model->poly_meshes)
            for (const vec3f& vtx : mesh->vertices3d)
                model->bounds.extend(vtx);


            // load some global settings (for now into model)
            model->useMaxwell = wp.useMaxwellDistribution;


        return model;
    }*/
}
