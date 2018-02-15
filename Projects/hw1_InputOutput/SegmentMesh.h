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
    std::vector< std::tuple<uint32_t, uint32_t, T> > indices;
//    T kHat; //young's modulus, E
//    T damp;//damping coefficient, b, gamma
    T L0; //spring's rest length
//    T c;
//    T dt;
//    T gravity;
//    T mass;

public:
    SegmentMesh();
    SegmentMesh(const int N);
    ~SegmentMesh();


    void InitVert(const T& height);

    void Update(const T& frameTime, const uint32_t iter);
    void WritePoly_Sim(const std::string& polyFileName, const uint32_t numFrames);
//    void InitFromHW3Line(const std::string& inputLine);
//    void SetConnections_HW3();
//    void WriteForces(const std::string& output);
//    void WriteForces(std::ofstream& outputFile);
//    void WritePoly_RandomFrames(const std::string& polyFileName);
};


#endif //CISPBA_SEGMENTMESH_H
