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
    Polygon( std::vector<osc::vec3f> vertex):
            vertices( vertex )
    {

    }
    std::vector<osc::vec3f> vertices;
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

/*! \namespace osc - Optix Siggraph Course */
namespace osc {



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

    Model *loadOBJ(const std::string &objFile){
        Model *model = new Model;

        const std::string mtlDir
                = objFile.substr(0,objFile.rfind('/')+1);
        PRINT(mtlDir);

        tinyobj::attrib_t attributes;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        std::string err = "";

        bool readOK
                = tinyobj::LoadObj(&attributes,
                                   &shapes,
                                   &materials,
                                   &err,
                                   &err,
                                   objFile.c_str(),
                                   mtlDir.c_str(),
                                    true); //triangulate
        if (!readOK) {
            throw std::runtime_error("Could not read OBJ model from "+objFile+":"+mtlDir+" : "+err);
        }

        if (materials.empty())
            throw std::runtime_error("could not parse materials ...");

        std::cout << "Done loading obj file - found " << shapes.size() << " shapes with " << materials.size() << " materials" << std::endl;
        for (int shapeID=0;shapeID<(int)shapes.size();shapeID++) {
            tinyobj::shape_t &shape = shapes[shapeID];

            std::set<int> materialIDs;
            for (auto faceMatID : shape.mesh.material_ids)
                materialIDs.insert(faceMatID);

            std::map<tinyobj::index_t,int> knownVertices;

            for (int materialID : materialIDs) {
                TriangleMesh *mesh = new TriangleMesh;

                for (int faceID=0;faceID<shape.mesh.material_ids.size();faceID++) {
                    if (shape.mesh.material_ids[faceID] != materialID) continue;
                    tinyobj::index_t idx0 = shape.mesh.indices[3*faceID+0];
                    tinyobj::index_t idx1 = shape.mesh.indices[3*faceID+1];
                    tinyobj::index_t idx2 = shape.mesh.indices[3*faceID+2];

                    vec3i idx(addVertex(mesh, attributes, idx0, knownVertices),
                              addVertex(mesh, attributes, idx1, knownVertices),
                              addVertex(mesh, attributes, idx2, knownVertices));
                    mesh->index.push_back(idx);
                    mesh->diffuse = (const vec3f&)materials[materialID].diffuse;
                    mesh->diffuse = gdt::randomColor(materialID);
                }

                if (mesh->vertex.empty())
                    delete mesh;
                else
                    model->meshes.push_back(mesh);
            }
        }

        // of course, you should be using tbb::parallel_for for stuff
        // like this:
        for (auto mesh : model->meshes)
            for (auto vtx : mesh->vertex)
                model->bounds.extend(vtx);

        std::cout << "created a total of " << model->meshes.size() << " meshes" << std::endl;
        return model;
    }

    Model *loadFromMolflow(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures){
        Model *model = new Model;
        // Loop over Vertices3

        PolygonMesh *mesh = new PolygonMesh;
        mesh->nbFacets = 0;
        for (const SuperStructure& structure : structures){
            mesh->nbFacets += structure.facets.size();
        }

        mesh->vertices3d.resize(geomVertices.size());
        mesh->poly.reserve(mesh->nbFacets);

        int polyNb = 0;
        for (const SuperStructure& structure : structures){
            //mesh->nbFacets += structure.facets.size();
            for(const SubprocessFacet& fac : structure.facets){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values

                Polygon newPoly;
                newPoly.nbIndices = fac.indices.size();
                newPoly.nbVertices = fac.vertices2.size();
                newPoly.indices = new int32_t[fac.indices.size()];
                newPoly.vertices2d = new vec2f[fac.vertices2.size()];
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


                //
                int counter = 0;

                for(auto vertInd : fac.indices){
                    newPoly.indices[counter++] = vertInd;

                    // load vertex corresponding to index
                    vec3f newVert;
                    newVert.x = geomVertices[vertInd].x;
                    newVert.y = geomVertices[vertInd].y;
                    newVert.z = geomVertices[vertInd].z;

                    mesh->vertices3d[vertInd] = std::move(newVert);
                }

                counter = 0;
                for(auto vert : fac.vertices2){
                    newPoly.vertices2d[counter++].u = vert.u;
                    newPoly.vertices2d[counter++].v = vert.v;
                }
                mesh->poly.push_back(std::move(newPoly));

                polyNb++;
            }
        }

        /*for(const Vector3d& loadedVertex : geomVertices){
            // Load Molflow style vertex into optix mesh
            vec3f newVert;
            newVert.x = loadedVertex.x;
            newVert.y = loadedVertex.y;
            newVert.z = loadedVertex.z;
            mesh->vertices3d.push_back(std::move(newVert));
        }*/

        if (mesh->vertices3d.empty())
            delete mesh;
        else
            model->poly_meshes.push_back(mesh);

        for (PolygonMesh* mesh : model->poly_meshes)
            for (const vec3f& vtx : mesh->vertices3d)
                model->bounds.extend(vtx);

        return model;
    }
}
