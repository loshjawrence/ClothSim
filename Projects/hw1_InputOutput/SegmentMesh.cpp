//
// Created by josh on 2/1/18.
//

#include "SegmentMesh.h"
#include <fstream>
#include <sstream>

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh() : particles(), indices()
{

}

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh(const int N) : particles(N), indices(N-1)
{

    particles.kHat = 20.0; //1.0
    particles.damp = 1.0; //2.0
    L0   = 1.0;   //m 1.0
    particles.dt   = 0.002; //s 0.002
    particles.mass = 0.3;  // kg 0.001
    particles.gravity = -9.81;  //m/s^2

    for(u_int32_t i = 0; i < indices.size(); ++i) {
        indices[i] = std::tuple<u_int32_t, u_int32_t, T>(i, i+1, L0);
    }
}

template<class T, int dim>
SegmentMesh<T, dim>::~SegmentMesh() {

}

template<class T, int dim>
void SegmentMesh<T, dim>::InitVert(const T& height) {
    particles.InitVert(height, L0);
}




template<class T, int dim>
void SegmentMesh<T, dim>::Update(const T& frameTime, const uint32_t iter) {
    for(T time = 0.0; time < frameTime; time += particles.dt) {
        particles.UpdateFE(indices, {});
    }
}

template<class T, int dim>
void SegmentMesh<T,dim>::WritePoly_Sim(const std::string& polyFileName, const uint32_t numFrames) {
    const T frameTime = 1.f / 24.f;
    //Open the debug file
    std::ofstream file;
    file.open("FramesOutput/Debug/segment/segmesh");
    for(uint32_t i = 1; i <= numFrames; ++i) {
        particles.WritePoly(polyFileName + std::to_string(i),indices);
        particles.WritePartio("FramesOutput/ParticleFrames/particles" + std::to_string(i));
//        particles.printForcePosVel(file, i, 0.0, indices);
        Update(frameTime, i);
    }
    file.close();
}

//
//template<class T, int dim>
//void SegmentMesh<T, dim>::WriteForces(const std::string& output) {
//    particles.WriteForces(output);
//}
//
//template<class T, int dim>
//void SegmentMesh<T, dim>::WriteForces(std::ofstream& outputFile) {
//    particles.WriteForces(outputFile);
//}

//template<class T, int dim>
//void SegmentMesh<T,dim>::WritePoly_RandomFrames(const std::string& polyFileName) {
//    const int numFrames = 120;
//    const T time_step = 1.0/24.0;
//    for(int i = 1; i <= numFrames; ++i) {
//        particles.WritePoly(polyFileName + std::to_string(i), indices);
//        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
//    }
//}

//
//template<class T, int dim>
//void SegmentMesh<T,dim>::SetConnections_HW3() {
//    indices.emplace_back(std::tuple<uint32_t, uint32_t, T>(0,1,L0));
//    indices.emplace_back(std::tuple<uint32_t, uint32_t, T>(1,2,L0));
//    indices.emplace_back(std::tuple<uint32_t, uint32_t, T>(2,3,L0));
//    indices.emplace_back(std::tuple<uint32_t, uint32_t, T>(3,0,L0));
//    indices.emplace_back(std::tuple<uint32_t, uint32_t, T>(3,1,L0));
//}

//template<>
//void SegmentMesh<double, 3>::InitFromHW3Line(const std::string &inputLine) {
//    std::istringstream iss(inputLine);
//
//    iss >> particles.kHat;
//    iss >> particles.damp;
//    iss >> L0;
//
//    //vars for holding the 3 components of a pos or vel
//    double x, y, z;
//
//    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
//    for(int i = 0; i < 4; ++i) {
//        iss>>x; iss>>y; iss>>z;
//        particles.pos.emplace_back(Eigen::Matrix<double, 3, 1>(x, y, z));
//        particles.mass = 1.0;
//        particles.springForce.emplace_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
//        particles.dampForce.emplace_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
//    }
//
//    //..followed by their 4 corresponding velocities
//    for(int i = 0; i < 4; ++i) {
//        iss>>x; iss>>y; iss>>z;
//        particles.vel.emplace_back(Eigen::Matrix<double, 3, 1>(x, y, z));
//    }
//    //setup connections between particles (predefined by hw diagram)
//    SetConnections_HW3();
//}
//
//template<>
//void SegmentMesh<float, 3>::InitFromHW3Line(const std::string &inputLine) {
//    std::istringstream iss(inputLine);
//
//    iss >> particles.kHat;
//    iss >> particles.damp;
//    iss >> L0;
//
//    //vars for holding the 3 components of a pos or vel
//    float x, y, z;
//
//    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
//    for(int i = 0; i < 4; ++i) {
//        iss>>x; iss>>y; iss>>z;
//        particles.pos.emplace_back(Eigen::Matrix<float, 3, 1>(x, y, z));
//        particles.mass=1;
//        particles.springForce.emplace_back(Eigen::Matrix<float, 3, 1>(0, 0, 0));
//        particles.dampForce.emplace_back(Eigen::Matrix<float, 3, 1>(0, 0, 0));
//    }
//
//    //..followed by their 4 corresponding velocities
//    for(int i = 0; i < 4; ++i) {
//        iss>>x; iss>>y; iss>>z;
//        particles.vel.emplace_back(Eigen::Matrix<float, 3, 1>(x, y, z));
//    }
//
//    //setup connections between particles (predefined by hw diagram)
//    SetConnections_HW3();
//}

//template class SegmentMesh<float,  2>;
//template class SegmentMesh<double, 2>;
template class SegmentMesh<float,  3>;
template class SegmentMesh<double, 3>;