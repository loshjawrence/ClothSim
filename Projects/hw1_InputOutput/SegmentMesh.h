//
// Created by josh on 2/1/18.
//

#ifndef CISPBA_SEGMENTMESH_H
#define CISPBA_SEGMENTMESH_H
#include "Particles.h"


template<class T, int dim>
class SegmentMesh {

public:
    Particles<T, dim> particles;

    //Every adjacent pair in "indices" is a pair of indices in the particles container.
    //A pair of indices represents a single segment
    vector<u_int32_t> indices;

public:
    SegmentMesh();
    SegmentMesh(const int N);
    ~SegmentMesh();

    void WriteToPoly(const string& polyFileName);
};


#endif //CISPBA_SEGMENTMESH_H
