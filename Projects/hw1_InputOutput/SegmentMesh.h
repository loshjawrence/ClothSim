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
    std::vector< std::pair<u_int32_t, u_int32_t> > indices;
    T kHat; //young's modulus, E
    T damp;//damping coefficient, b, gamma
    T L0; //spring's rest length

public:
    SegmentMesh();
    SegmentMesh(const int N);
    ~SegmentMesh();

    void InitFromHW3Line(const std::string& inputLine);
    void SetConnections_HW3();
    void CalcForces();
    void WriteForces(const std::string& output);
    void WriteForces(std::ofstream& outputFile);
    void WritePoly_RandomFrames(const std::string& polyFileName);
    void WritePoly(const std::string& polyFileName);
};


#endif //CISPBA_SEGMENTMESH_H
