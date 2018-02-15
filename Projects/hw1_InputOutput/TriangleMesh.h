//
// Created by josh on 2/1/18.
//

#ifndef CISPBA_TRIANGLEMESH_H
#define CISPBA_TRIANGLEMESH_H

#include "Particles.h"
#include <tuple>

template<class T, int dim>
class TriangleMesh {

public:
    Particles<T, dim> particles;

    //similar to index buffer object in graphics api's,
    //each set of 3 indices is a triangle face
    std::vector<uint32_t> indices;
    std::vector<std::tuple<uint32_t, uint32_t, T>> springs;
    std::vector<std::tuple<uint32_t, uint32_t, T>> bendSprings;

    int32_t width; //num quads in x
    int32_t height;//num quads in y
//    T mass;  // kg
//    T kHat;  // N/m, young's modulus E
//    T c;     // kg/s
//    T gravity;  //m/s^2
//    T dt;    //s
    T L0;    //spring rest length, horiz, vert
    T L0_d;  //spring rest length, diag
//    T damp;  //damping coeff, gamma, b


public:
    TriangleMesh();
    TriangleMesh(const int N);
    TriangleMesh(const int X, const int Y);//cloth, quads X quads
    ~TriangleMesh();

    void Update(const T& frameTime, const uint32_t iter);
    void printForcePosVel(std::ofstream& file, const int iter, const T& time);
    void CalcForces();
    void AssignTriangleFanIndices();
    void AssignGridMeshPositionsIndices();
    void AssignSpringArray();
    void WriteObj(const std::string& objFileName);
    void WriteObj_RandomFrames(const std::string& objFileName, const uint32_t numFrames);
    void WriteObj_Sim(const std::string& objFileName, const uint32_t numFrames);

};


#endif //CISPBA_TRIANGLEMESH_H
