//
// Created by josh on 2/1/18.
//

#include "SegmentMesh.h"
#include <fstream>

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh() : particles(), indices()
{

}

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh(const int N) : particles(N), indices(N)
{
    for(u_int32_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
}

template<class T, int dim>
SegmentMesh<T, dim>::~SegmentMesh()
{

}

template<class T, int dim>
void SegmentMesh<T,dim>::WritePoly_RandomFrames(const string& polyFileName) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WritePoly(polyFileName + to_string(i));
        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }
}

template<class T, int dim>
void SegmentMesh<T, dim>::WritePoly(const string& polyFileName) {
    //Open the file
    ofstream polyFile;
    polyFile.open(polyFileName + ".poly");

    //Write the POINTS section
    polyFile << "POINTS\n";
    for(u_int32_t i = 0; i < particles.pos.size(); ++i) {
        polyFile << i+1 << ":";
        for(int k = 0; k < dim; ++k) {
            polyFile << " " << particles.pos[i][k];
        }
        polyFile << "\n";
    }

    //Write the POLYS section, each entry is a string of connected points
    polyFile << "POLYS\n1:";
    for(u_int32_t i = 0; i < indices.size(); ++i) {
        polyFile << " " << indices[i]+1;
    }
    polyFile << "\n";

    //Write the END section
    polyFile << "END";
    polyFile.close();
}

template class SegmentMesh<float,  2>;
template class SegmentMesh<float,  3>;
template class SegmentMesh<double, 2>;
template class SegmentMesh<double, 3>;