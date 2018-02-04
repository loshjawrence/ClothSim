//
// Created by josh on 2/1/18.
//

#include "TriangleMesh.h"
#include <fstream>

template<class T, int dim>
TriangleMesh<T,dim>::TriangleMesh() : particles(), indices()
{

}


template<class T, int dim>
TriangleMesh<T, dim>::TriangleMesh(const int N) : particles(N), indices(3*(N-2))
{
    //This constructor assumes simple triangle fan construction
    AssignTriangleFanIndices();
}

template<class T, int dim>
TriangleMesh<T, dim>::~TriangleMesh()
{

}

template<class T, int dim>
void TriangleMesh<T, dim>::AssignTriangleFanIndices() {
    const int numTris = particles.pos.size() - 2;
    int fanIndex = 1;
    for(int i = 0; i < numTris; ++i) {
        const int baseIndex = 3*i;
        indices[baseIndex    ] = 0;
        indices[baseIndex + 1] = fanIndex;
        indices[baseIndex + 2] = ++fanIndex;
    }
}

template<class T, int dim>
void TriangleMesh<T, dim>::WriteObj(const std::string& filename) {
    //Open the file
    std::ofstream file;
    file.open(filename + ".obj");

    //Write the v (vertex) section
    for(u_int32_t i = 0; i < particles.pos.size(); ++i) {
        file << "\nv";
        for(int k = 0; k < dim; ++k) {
            file << " " << particles.pos[i][k];
        }
    }

    //write the f (face) section
    const int numTris = indices.size() / 3;
    for(int i = 0; i < numTris; ++i) {
        const int baseIndex = 3*i;
        file << "\nf " << indices[baseIndex    ]+1;
        file << " "  << indices[baseIndex + 1]+1;
        file << " "  << indices[baseIndex + 2]+1;
    }

    file.close();

}

template<class T, int dim>
void TriangleMesh<T, dim>::WriteObj_RandomFrames(const std::string& objFileName) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WriteObj(objFileName + std::to_string(i));
        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }

}


template class TriangleMesh<float,  2>;
template class TriangleMesh<float,  3>;
template class TriangleMesh<double, 2>;
template class TriangleMesh<double, 3>;
