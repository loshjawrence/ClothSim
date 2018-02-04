//
// Created by josh on 2/1/18.
//

#ifndef CISPBA_TRIANGLEMESH_H
#define CISPBA_TRIANGLEMESH_H

#include "Particles.h"

template<class T, int dim>
class TriangleMesh {

public:
    Particles<T, dim> particles;

    //similar to index buffer object in graphics api's,
    //each set of 3 indices is a triangle face
    std::vector<u_int32_t> indices;

public:
    TriangleMesh();
    TriangleMesh(const int N);
    ~TriangleMesh();

    void AssignTriangleFanIndices();
    void WriteObj(const std::string& objFileName);
    void WriteObj_RandomFrames(const std::string& objFileName);

};


#endif //CISPBA_TRIANGLEMESH_H
