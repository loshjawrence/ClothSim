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
void SegmentMesh<T, dim>::WriteToPoly(const string& polyFileName) {
    //Open the file
    ofstream polyFile;
    polyFile.open(polyFileName);

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