//
// Created by pbahr on 08/10/2019.
//

#include "Model.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "3rdParty/tiny_obj_loader.h"
//std
#include <set>
#include <Facet_shared.h>

/*struct Polygon
{
    Polygon() = default;
    Polygon( std::vector<flowgpu::vec3f> vertex):
            vertices( vertex )
    {

    }
    std::vector<flowgpu::vec3f> vertices;
};*/

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



    /*! find vertex with given position, normal, texcoord, and return
        its vertex ID, or, if it doesn't exit, add it to the mesh, and
        its just-created index */
    int addVertex(TriangleMesh *mesh,
                  tinyobj::attrib_t &attributes,
                  const tinyobj::index_t &idx,
                  std::map<tinyobj::index_t,int> &knownVertices){

        if (knownVertices.find(idx) != knownVertices.end())
            return knownVertices[idx];

        const vec3f *vertex_array   = (const vec3f*)attributes.vertices.data();
        const vec3f *normal_array   = (const vec3f*)attributes.normals.data();
        const vec2f *texcoord_array = (const vec2f*)attributes.texcoords.data();

        int newID = mesh->vertex.size();
        knownVertices[idx] = newID;

        mesh->vertex.push_back(vertex_array[idx.vertex_index]);
        if (idx.normal_index >= 0) {
            while (mesh->normal.size() < mesh->vertex.size())
                mesh->normal.push_back(normal_array[idx.normal_index]);
        }
        if (idx.texcoord_index >= 0) {
            while (mesh->texcoord.size() < mesh->vertex.size())
                mesh->texcoord.push_back(texcoord_array[idx.texcoord_index]);
        }

        // just for sanity's sake:
        if (mesh->texcoord.size() > 0)
            mesh->texcoord.resize(mesh->vertex.size());
        // just for sanity's sake:
        if (mesh->normal.size() > 0)
            mesh->normal.resize(mesh->vertex.size());

        return newID;
    }

    Model *loadFromMolflow(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures){
        Model *model = new Model;
        // Loop over Vertices3

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
            //std::cout << "---- "<<mesh->nbFacets<<" -- "<<mesh->nbIndices<<" -- "<<mesh->nbVertices<< std::endl;

            int32_t vertCount = 0;
            for( int f = 0; f < structures[s].facets.size(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                const SubprocessFacet& fac = structures[s].facets[f];
                Polygon newPoly(fac.indices.size());
                newPoly.vertOffset = vertCount;
                newPoly.nbVertices = fac.indices.size();

                newPoly.O.x = fac.sh.O.x;
                newPoly.O.y = fac.sh.O.y;
                newPoly.O.z = fac.sh.O.z;
                newPoly.U.x = fac.sh.U.x;
                newPoly.U.y = fac.sh.U.y;
                newPoly.U.z = fac.sh.U.z;
                newPoly.V.x = fac.sh.V.x;
                newPoly.V.y = fac.sh.V.y;
                newPoly.V.z = fac.sh.V.z;
                newPoly.Nuv.x = fac.sh.Nuv.x;
                newPoly.Nuv.y = fac.sh.Nuv.y;
                newPoly.Nuv.z = fac.sh.Nuv.z;

                int counter = 0;
                for(size_t vertInd : fac.indices){

                    // load vertex corresponding to index
                    vec3f newVert;
                    newVert.x = geomVertices[vertInd].x;
                    newVert.y = geomVertices[vertInd].y;
                    newVert.z = geomVertices[vertInd].z;

                    //newVert *= 1.0f / dot(newVert,newVert);

                   // std::cout << newVert << " -> " << vertInd<<"/"<<mesh->vertices3d.size()<<" -- Vert push"<<std::endl;
                    mesh->vertices3d[vertInd] = newVert;
                    mesh->indices[vertCount] = vertInd;
                    mesh->vertices2d[vertCount] = vec2f(fac.vertices2[counter].u, fac.vertices2[counter].v);
                    counter++;
                    vertCount++;
                }

                //std::cout << f<< " -- "<<counter<<" / "<<vertCount<<" --> "<<newPoly.vertOffset<< std::endl;

                mesh->poly[f] = std::move(newPoly);

            }

            // final sanity check
            if (mesh->vertices3d.empty() || mesh->indices.empty() || mesh->vertices2d.empty())
                delete mesh;
            else
                model->poly_meshes.push_back(mesh);
        }

        // calculate global bounds for the whole model
        for (PolygonMesh* mesh : model->poly_meshes)
            for (const vec3f& vtx : mesh->vertices3d)
                model->bounds.extend(vtx);

        return model;
    }
}
