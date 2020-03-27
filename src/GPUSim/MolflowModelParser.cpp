//
// Created by pbahr on 15/11/2019.
//

#include "MolflowModelParser.h"
//#include "Facet_shared.h"

// debug output
#include <fstream>

#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/tuple.hpp>
#ifdef MOLFLOW
#include "Geometry_shared.h"
#include "Facet_shared.h"
#endif MOLFLOW
/*#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>*/

/*! \namespace flowgeom - Molflow Geometry code */
namespace flowgeom {

    /*inline void vector3d_to_float3(gdt::float3& t, const Vector3d& o){
        t.x = o.x;
        t.y = o.y;
        t.z = o.z;

        return;
    }*/

/*    int loadFromMolflow(Geometry& geometry, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs){
        // Loop over Vertices3

        *//*std::vector<SuperStructure> structures;

        //Facets
        for (size_t i = 0; i < geometry.sh.nbFacet; i++) { //Necessary because facets is not (yet) a vector in the interface
            Facet* fac = geometry.facets[i];
            SubprocessFacet subf;

            subf.sh = fac->sh;
            subf.indices = fac->indices;
            subf.vertices2 = fac->vertices2;

            if (subf.sh.superIdx == -1) { //Facet in all structures
                for (auto& s : structures) {
                    s.facets.push_back(subf);
                }
            }
            else {
                structures[subf.sh.superIdx].facets.push_back(subf); //Assign to structure
            }
        }*//*

        std::vector<Vector3d> geomVertices;

        for (size_t i = 0; i < geometry.sh.nbVertex; i++) {
            Vector3d vec;
            InterfaceVertex* interfaceVec = geometry.GetVertex(i);
            vec = *interfaceVec;
            geomVertices.push_back(vec);
        }

        saveFromMolflow(geomVertices, structures, wp, CDFs);

        return 0;
    }*/

    /*flowgpu::Model *loadFromMolflowSimu(const std::vector<Vector3d> &geomVertices, const std::vector<SuperStructure> &structures, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs){
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
            mesh->nbVertices = indNb;

            std::vector<float3>(geomVertices.size()).swap(mesh->vertices3d);
            mesh->poly.resize(mesh->nbFacets);
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
                newPoly.indexOffset = vertCount;
                newPoly.nbVertices = fac.indices.size();
                newPoly.stickingFactor = fac.sh.sticking;

                vector3d_to_float3(newPoly.O, fac.sh.O);
                vector3d_to_float3(newPoly.U, fac.sh.U);
                vector3d_to_float3(newPoly.V, fac.sh.V);
                vector3d_to_float3(newPoly.Nuv, fac.sh.Nuv);
                vector3d_to_float3(newPoly.nU, fac.sh.nU);
                vector3d_to_float3(newPoly.nV, fac.sh.nV);
                vector3d_to_float3(newPoly.N, fac.sh.N);

                int counter = 0;
                for(size_t vertInd : fac.indices){

                    // load vertex corresponding to index
                    float3 newVert;
                    vector3d_to_float3(newVert, geomVertices[vertInd]);

                    //newVert *= 1.0f / dot(newVert,newVert);

                    mesh->vertices3d[vertInd] = newVert;
                    mesh->indices[vertCount] = vertInd;
                    mesh->vertices2d[vertCount] = float2(fac.vertices2[counter].u, fac.vertices2[counter].v);
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
        model->nbVertices_total = 0;
        for (PolygonMesh* mesh : model->poly_meshes){
            model->nbFacets_total += mesh->nbFacets;
            model->nbVertices_total += mesh->nbVertices;
        }
        // calculate global bounds for the whole model
        for (PolygonMesh* mesh : model->poly_meshes)
            for (const float3& vtx : mesh->vertices3d)
                model->bounds.extend(vtx);


        // load some global settings (for now into model)
        model->useMaxwell = wp.useMaxwellDistribution;


        std::ofstream file( "test_geom.xml" );
        cereal::XMLOutputArchive archive( file );
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
                cereal::make_nvp("useMaxwell", model->useMaxwell) ,
                cereal::make_nvp("bounds.lower", model->bounds.lower) ,
                cereal::make_nvp("bounds.upper", model->bounds.upper)
        );

        return model;
    }*/

    struct Polygon {
    public:
        Polygon()
                : stickingFactor(-1.0), nbVertices(0), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(){
        }
        Polygon(int32_t nbOfVertices)
                : stickingFactor(-1.0), nbVertices(nbOfVertices), indexOffset(0), O(), U(), V(), Nuv(), nU(), nV(), N(){
        }
        // attributes that don't describe the geometry
        float stickingFactor;

        // variables for access to  global memory (indices, vertices)
        uint32_t nbVertices;
        uint32_t indexOffset;

        // variables for ray-plane (3d space) intersection
        Vector3d O;
        Vector3d U;
        Vector3d V;
        Vector3d Nuv;

        // normalized facet vectors
        Vector3d nU;
        Vector3d nV;
        Vector3d N;

        template<class Archive>
        void serialize(Archive & archive)
        {
            archive(
                    cereal::make_nvp("stickingFactor", stickingFactor) ,
                    cereal::make_nvp("nbVertices", nbVertices) ,
                    cereal::make_nvp("indexOffset", indexOffset) ,
                    cereal::make_nvp("O", O) ,
                    cereal::make_nvp("U", U) ,
                    cereal::make_nvp("V", V) ,
                    cereal::make_nvp("Nuv", Nuv) ,
                    cereal::make_nvp("nU", nU) ,
                    cereal::make_nvp("nV", nV) ,
                    cereal::make_nvp("N", N)
            );
        }
    };

    int saveFromMolflowTriangle(Geometry& geometry, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs){

        std::vector<std::tuple<uint32_t, uint32_t, uint32_t> > indices;
        std::vector<Vector2d> vertices2d;
        std::vector<Vector3d> vertices3d;
        std::vector<Polygon> poly;

        std::vector<Vector2d> facetProbabilities;
        std::vector<float> cdfs; // first pair value: ind%2==0 .. second pair value: ind%2==1

        uint32_t nbFacets = 0;
        uint32_t nbVertices = 0;

        int polyNb = 0;
        int indNb = 0;
        int vertNb = 0;

        polyNb = geometry.GetNbFacet();

        // Count total indices and vertices
        for (size_t i = 0; i < geometry.GetNbFacet(); i++) {
            Facet* fac = geometry.GetFacet(i);

            indNb += fac->indices.size();
            vertNb += fac->vertices2.size();
        }

        if(indNb%3 != 0){
            std::cout << "[ERROR] Amount of indices doesn't suggest triangulated geometry: "<<indNb<<"!"<<std::endl;
            return 1;
        }

        // Convert InterfaceVertices to normal Vertices
        std::vector<Vector3d> geomVertices;
        for (size_t i = 0; i < geometry.GetNbVertex(); i++) {
            Vector3d vec;
            InterfaceVertex* interfaceVec = geometry.GetVertex(i);
            vec = *interfaceVec;
            geomVertices.push_back(vec);
        }

        //for ( int s = 0; s < structures.size(); s++){
        for ( int s = 0; s <= 0; s++){ // no support for multi structures for now

            //TODO: one mesh per structure

            nbFacets = polyNb;
            nbVertices = indNb;

            std::vector<Vector3d>(geomVertices.size()).swap(vertices3d);
            poly.resize(nbFacets);
            indices.resize(nbVertices/3);
            vertices2d.resize(nbVertices);

            int size_cdf = 0;
            for(auto& cdf : CDFs){
                size_cdf += 2*cdf.size(); // size for temp bins, v bins and both pair values
            }
            cdfs.resize(size_cdf);
            facetProbabilities.resize(nbFacets);

            int32_t indexCount = 0;
            double fullOutgassing = 0;



            //for( int f = 0; f < structures[s].facets.size(); f++){
            for( int f = 0; f < geometry.GetNbFacet(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                //const Facet& fac = structures[s].facets[f];
                Facet* fac = geometry.GetFacet(f);

                //Polygon newPoly(fac.indices.size());
                poly[f] = Polygon(fac->indices.size());
                poly[f].indexOffset = indexCount;
                poly[f].nbVertices = fac->indices.size();
                if(poly[f].nbVertices != 3){
                    std::cout << "[ERROR] Facet with "<<poly[f].nbVertices<<" vertices loaded, while looking for triangles!"<<std::endl;
                    return 1;
                }
                poly[f].stickingFactor = fac->sh.sticking;

                poly[f].O = fac->sh.O;
                poly[f].U = fac->sh.U;
                poly[f].V = fac->sh.V;
                poly[f].Nuv = fac->sh.Nuv;
                poly[f].nU = fac->sh.nU;
                poly[f].nV = fac->sh.nV;
                poly[f].N = fac->sh.N;

                int counter = 0;
                for(size_t vertInd : fac->indices){

                    // load vertex corresponding to index
                    Vector3d newVert;
                    newVert = geomVertices[vertInd];

                    //newVert *= 1.0f / dot(newVert,newVert);

                    vertices3d[vertInd] = newVert;
                    //indices[indexCount] = vertInd;
                    //indexCount++;
                }
                // Set index triplet
                indices[f] = std::make_tuple(fac->indices.at(0), fac->indices.at(1), fac->indices.at(2));

                double facOutgassing = wp.latestMoment*fac->sh.outgassing / (1.38E-23*fac->sh.temperature);
                facetProbabilities[f].u = fullOutgassing; // lower border
                fullOutgassing += facOutgassing;
                facetProbabilities[f].v = fullOutgassing; // upper border
            }

            //for( int f = 0; f < structures[s].facets.size(); f++){
            for( int f = 0; f < geometry.GetNbFacet(); f++){

                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                //const SubprocessFacet& fac = structures[s].facets[f];
                facetProbabilities[f].u /= fullOutgassing; // normalize from [0,1]
                facetProbabilities[f].v /= fullOutgassing; // normalize from [0,1]
                //std::cout << "2. "<< fullOutgassing << " - ["<<mesh->facetProbabilities[f].s<<","<<mesh->facetProbabilities[f].t<<"]"<<std::endl;

            }

            int index = 0;
            for(auto& cdf : CDFs){
                for(auto& pair : cdf){
                    cdfs[index++] = pair.first; // %2==0
                    cdfs[index++] = pair.second; // %2==1
                }
            }

            // final sanity check
            if (vertices3d.empty() || indices.empty() || vertices2d.empty())
                return -1;
            /*else
                model->poly_meshes.push_back(mesh);*/
        }

        // calculate total amount of facets
        uint32_t nbFacets_total = 0;
        uint32_t nbVertices_total = 0;
        //for (PolygonMesh* mesh : model->poly_meshes){
        nbFacets_total += nbFacets;
        nbVertices_total += nbVertices;
        //}

        Vector3d lower(std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
        Vector3d upper(std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min());
        // calculate global bounds for the whole model
        //for (PolygonMesh* mesh : model->poly_meshes)
        for (const Vector3d& vtx : vertices3d){
            //bounds.extend(vtx);
            if(vtx.x > upper.x) upper.x = vtx.x;
            if(vtx.y > upper.y) upper.y = vtx.y;
            if(vtx.z > upper.z) upper.z = vtx.z;

            if(vtx.x < lower.x) lower.x = vtx.x;
            if(vtx.y < lower.y) lower.y = vtx.y;
            if(vtx.z < lower.z) lower.z = vtx.z;
        }
        // load some global settings (for now into model)
        bool useMaxwell = wp.useMaxwellDistribution;

        std::cout << "Writing parsed geometry to file: test_geom_tri.xml" << std::endl;
        std::ofstream file( "test_geom_tri.xml" );
        cereal::XMLOutputArchive archive( file );
        archive(
                cereal::make_nvp("poly", poly) ,
                cereal::make_nvp("facetProbabilities", facetProbabilities) ,
                cereal::make_nvp("cdfs", cdfs) ,
                //cereal::make_nvp("vertices2d", vertices2d) ,
                cereal::make_nvp("vertices3d", vertices3d) ,
                cereal::make_nvp("indices", indices) ,
                cereal::make_nvp("nbFacets", nbFacets) ,
                cereal::make_nvp("nbVertices", nbVertices) ,
                cereal::make_nvp("nbFacetsTotal", nbFacets_total) ,
                cereal::make_nvp("nbVerticesTotal", nbVertices_total) ,
                cereal::make_nvp("useMaxwell", useMaxwell)
        );

        return 0;
    }

    int saveFromMolflow(Geometry& geometry, const WorkerParams& wp, const std::vector<std::vector<std::pair<double,double>>> CDFs){

        std::vector<uint32_t> indices;
        std::vector<Vector2d> vertices2d;
        std::vector<Vector3d> vertices3d;
        std::vector<Polygon> poly;

        std::vector<Vector2d> facetProbabilities;
        std::vector<float> cdfs;

        uint32_t nbFacets = 0;
        uint32_t nbVertices = 0;

        int polyNb = 0;
        int indNb = 0;
        int vertNb = 0;

        polyNb = geometry.GetNbFacet();

        for (size_t i = 0; i < geometry.GetNbFacet(); i++) { //Necessary because facets is not (yet) a vector in the interface
            Facet* fac = geometry.GetFacet(i);

            indNb += fac->indices.size();
            vertNb += fac->vertices2.size();

            /*if (subf.sh.superIdx == -1) { //Facet in all structures
                for (auto& s : structures) {
                    s.facets.push_back(subf);
                }
            }
            else {
                structures[subf.sh.superIdx].facets.push_back(subf); //Assign to structure
            }*/
        }

        // Convert InterfaceVertices to normal Vertices
        std::vector<Vector3d> geomVertices;
        for (size_t i = 0; i < geometry.GetNbVertex(); i++) {
            Vector3d vec;
            InterfaceVertex* interfaceVec = geometry.GetVertex(i);
            vec = *interfaceVec;
            geomVertices.push_back(vec);
        }

        //for ( int s = 0; s < structures.size(); s++){
        for ( int s = 0; s <= 0; s++){ // no support for multi structures for now

            //TODO: one mesh per structure

            nbFacets = polyNb;
            nbVertices = indNb;

            std::vector<Vector3d>(geomVertices.size()).swap(vertices3d);
            poly.resize(nbFacets);
            indices.resize(indNb);
            vertices2d.resize(nbVertices);

            int size_cdf = 0;
            for(auto& cdf : CDFs){
                size_cdf += 2*cdf.size(); // size for temp bins, v bins and both pair values
            }
            cdfs.resize(size_cdf);
            facetProbabilities.resize(nbFacets);

            int32_t vertCount = 0;
            double fullOutgassing = 0;

            
            
            /*for( int f = 0; f < structures[s].facets.size(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                const SubprocessFacet& fac = structures[s].facets[f];*/
            for( int f = 0; f < geometry.GetNbFacet(); f++){
                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                //const Facet& fac = structures[s].facets[f];
                Facet* fac = geometry.GetFacet(f);
                
                //Polygon newPoly(fac.indices.size());
                poly[f] = Polygon(fac->indices.size());
                poly[f].indexOffset = vertCount;
                poly[f].nbVertices = fac->indices.size();
                poly[f].stickingFactor = fac->sh.sticking;

                poly[f].O = fac->sh.O;
                poly[f].U = fac->sh.U;
                poly[f].V = fac->sh.V;
                poly[f].Nuv = fac->sh.Nuv;
                poly[f].nU = fac->sh.nU;
                poly[f].nV = fac->sh.nV;
                poly[f].N = fac->sh.N;

                int counter = 0;
                for(size_t vertInd : fac->indices){

                    // load vertex corresponding to index
                    Vector3d newVert;
                    newVert = geomVertices[vertInd];

                    //newVert *= 1.0f / dot(newVert,newVert);

                    vertices3d[vertInd] = newVert;
                    indices[vertCount] = vertInd;
                    vertices2d[vertCount] = Vector2d(fac->vertices2[counter].u, fac->vertices2[counter].v);
                    counter++;
                    vertCount++;
                }
                //poly[f] = std::move(newPoly);

                double facOutgassing = wp.latestMoment*fac->sh.outgassing / (1.38E-23*fac->sh.temperature);
                facetProbabilities[f].u = fullOutgassing; // lower border
                fullOutgassing += facOutgassing;
                facetProbabilities[f].v = fullOutgassing; // upper border
            }

            //for( int f = 0; f < structures[s].facets.size(); f++){
            for( int f = 0; f < geometry.GetNbFacet(); f++){

                // Create new temp polygon and initialize number of vert and ind
                // load with proper values
                //const SubprocessFacet& fac = structures[s].facets[f];
                facetProbabilities[f].u /= fullOutgassing; // normalize from [0,1]
                facetProbabilities[f].v /= fullOutgassing; // normalize from [0,1]
                //std::cout << "2. "<< fullOutgassing << " - ["<<mesh->facetProbabilities[f].s<<","<<mesh->facetProbabilities[f].t<<"]"<<std::endl;

            }

            int index = 0;
            for(auto& cdf : CDFs){
                for(auto& pair : cdf){
                    cdfs[index++] = pair.first; // %2==0
                    cdfs[index++] = pair.second; // %2==1
                }
            }

            // final sanity check
            if (vertices3d.empty() || indices.empty() || vertices2d.empty())
                return -1;
            /*else
                model->poly_meshes.push_back(mesh);*/
        }

        // calculate total amount of facets
        uint32_t nbFacets_total = 0;
        uint32_t nbVertices_total = 0;
        //for (PolygonMesh* mesh : model->poly_meshes){
            nbFacets_total += nbFacets;
            nbVertices_total += nbVertices;
        //}

        Vector3d lower(std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
        Vector3d upper(std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min());
        // calculate global bounds for the whole model
        //for (PolygonMesh* mesh : model->poly_meshes)
        for (const Vector3d& vtx : vertices3d){
            //bounds.extend(vtx);
            if(vtx.x > upper.x) upper.x = vtx.x;
            if(vtx.y > upper.y) upper.y = vtx.y;
            if(vtx.z > upper.z) upper.z = vtx.z;

            if(vtx.x < lower.x) lower.x = vtx.x;
            if(vtx.y < lower.y) lower.y = vtx.y;
            if(vtx.z < lower.z) lower.z = vtx.z;
        }
        // load some global settings (for now into model)
        bool useMaxwell = wp.useMaxwellDistribution;

        std::cout << "Writing parsed geometry to file: test_geom.xml" << std::endl;
        std::ofstream file( "test_geom.xml" );
        cereal::XMLOutputArchive archive( file );
        archive(
                cereal::make_nvp("poly", poly) ,
                cereal::make_nvp("facetProbabilities", facetProbabilities) ,
                cereal::make_nvp("cdfs", cdfs) ,
                cereal::make_nvp("vertices2d", vertices2d) ,
                cereal::make_nvp("vertices3d", vertices3d) ,
                cereal::make_nvp("indices", indices) ,
                cereal::make_nvp("nbFacets", nbFacets) ,
                cereal::make_nvp("nbVertices", nbVertices) ,
                cereal::make_nvp("nbFacetsTotal", nbFacets_total) ,
                cereal::make_nvp("nbVerticesTotal", nbVertices_total) ,
                cereal::make_nvp("useMaxwell", useMaxwell)
        );

        return 0;
    }
}