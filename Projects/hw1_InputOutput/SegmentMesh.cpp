//
// Created by josh on 2/1/18.
//

#include "SegmentMesh.h"
#include <fstream>
#include <sstream>

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh() : particles(), indices(), kHat(), damp(), L0()
{

}

template<class T, int dim>
SegmentMesh<T, dim>::SegmentMesh(const int N) : particles(N), indices(N-1), kHat(), damp(), L0()
{
    for(u_int32_t i = 0; i < indices.size(); ++i) {
        indices[i] = std::pair<u_int32_t, u_int32_t>(i, i+1);
    }
}

template<class T, int dim>
SegmentMesh<T, dim>::~SegmentMesh() {

}

template<>
void SegmentMesh<double, 3>::InitFromHW3Line(const std::string &inputLine) {
    std::istringstream iss(inputLine);

    iss >> kHat;
    iss >> damp;
    iss >> L0;

    //vars for holding the 3 components of a pos or vel
    double x, y, z;

    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y; iss>>z;
        particles.pos.emplace_back(Eigen::Matrix<double, 3, 1>(x, y, z));
        particles.mass.emplace_back(1);
        particles.springForce.emplace_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
        particles.dampForce.emplace_back(Eigen::Matrix<double, 3, 1>(0, 0, 0));
    }

    //..followed by their 4 corresponding velocities
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y; iss>>z;
        particles.vel.emplace_back(Eigen::Matrix<double, 3, 1>(x, y, z));
    }
    //setup connections between particles (predefined by hw diagram)
    SetConnections_HW3();
}

template<>
void SegmentMesh<float, 3>::InitFromHW3Line(const std::string &inputLine) {
    std::istringstream iss(inputLine);

    iss >> kHat;
    iss >> damp;
    iss >> L0;

    //vars for holding the 3 components of a pos or vel
    float x, y, z;

    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y; iss>>z;
        particles.pos.emplace_back(Eigen::Matrix<float, 3, 1>(x, y, z));
        particles.mass.emplace_back(1);
        particles.springForce.emplace_back(Eigen::Matrix<float, 3, 1>(0, 0, 0));
        particles.dampForce.emplace_back(Eigen::Matrix<float, 3, 1>(0, 0, 0));
    }

    //..followed by their 4 corresponding velocities
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y; iss>>z;
        particles.vel.emplace_back(Eigen::Matrix<float, 3, 1>(x, y, z));
    }

    //setup connections between particles (predefined by hw diagram)
    SetConnections_HW3();
}


template<>
void SegmentMesh<double, 2>::InitFromHW3Line(const std::string &inputLine) {
    std::istringstream iss(inputLine);

    iss >> kHat;
    iss >> damp;
    iss >> L0;

    //vars for holding the 3 components of a pos or vel
    double x, y;

    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y;
        particles.pos.emplace_back(Eigen::Matrix<double, 2, 1>(x, y));
        particles.mass.emplace_back(1);
        particles.springForce.emplace_back(Eigen::Matrix<double, 2, 1>(0, 0));
        particles.dampForce.emplace_back(Eigen::Matrix<double, 2, 1>(0, 0));
    }

    //..followed by their 4 corresponding velocities
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y;
        particles.vel.emplace_back(Eigen::Matrix<double, 2, 1>(x, y));
    }
    //setup connections between particles (predefined by hw diagram)
    SetConnections_HW3();
}

template<>
void SegmentMesh<float, 2>::InitFromHW3Line(const std::string &inputLine) {
    std::istringstream iss(inputLine);

    iss >> kHat;
    iss >> damp;
    iss >> L0;

    //vars for holding the 3 components of a pos or vel
    float x, y;

    //there are 4 point positions in the predefined network, init mass and force to 1 and 0
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y;
        particles.pos.emplace_back(Eigen::Matrix<float, 2, 1>(x, y));
        particles.mass.emplace_back(1);
        particles.springForce.emplace_back(Eigen::Matrix<float, 2, 1>(0, 0));
        particles.dampForce.emplace_back(Eigen::Matrix<float, 2, 1>(0, 0));
    }

    //..followed by their 4 corresponding velocities
    for(int i = 0; i < 4; ++i) {
        iss>>x; iss>>y;
        particles.vel.emplace_back(Eigen::Matrix<float, 2, 1>(x, y));
    }

    //setup connections between particles (predefined by hw diagram)
    SetConnections_HW3();
}


template<class T, int dim>
void SegmentMesh<T,dim>::SetConnections_HW3() {
    indices.emplace_back(std::pair<uint32_t, uint32_t>(0,1));
    indices.emplace_back(std::pair<uint32_t, uint32_t>(1,2));
    indices.emplace_back(std::pair<uint32_t, uint32_t>(2,3));
    indices.emplace_back(std::pair<uint32_t, uint32_t>(3,0));
    indices.emplace_back(std::pair<uint32_t, uint32_t>(3,1));
}

template<class T, int dim>
void SegmentMesh<T, dim>::CalcForces() {
    //zero forces
    for(Eigen::Matrix<T,dim,1>& force : particles.springForce) {
        force.setZero();
    }
    for(Eigen::Matrix<T,dim,1>& force : particles.dampForce) {
        force.setZero();
    }

    //loop through connections and sum up damp and spring forces
    for(auto& pair: indices) {
        //accumulate damping forces
        Eigen::Matrix<T,dim,1> n12 = particles.pos[pair.first] - particles.pos[pair.second];
        Eigen::Matrix<T,dim,1> v12 = particles.vel[pair.first] - particles.vel[pair.second];
        T l = n12.norm();
        n12.normalize();//or just divide by l
        Eigen::Matrix<T,dim,1> fd1 = -damp * n12.dot(v12) * n12;
        particles.dampForce[pair.first]  +=  fd1;
        particles.dampForce[pair.second] += -fd1;

        //accumulate spring force
        Eigen::Matrix<T,dim,1> fs1 = -kHat * (l/L0 - 1) * n12;
        particles.springForce[pair.first]  +=  fs1;
        particles.springForce[pair.second] += -fs1;
    }
}

template<class T, int dim>
void SegmentMesh<T, dim>::WriteForces(const std::string& output) {
    particles.WriteForces(output);
}

template<class T, int dim>
void SegmentMesh<T, dim>::WriteForces(std::ofstream& outputFile) {
    particles.WriteForces(outputFile);
}


template<class T, int dim>
void SegmentMesh<T,dim>::WritePoly_RandomFrames(const std::string& polyFileName) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WritePoly(polyFileName + std::to_string(i));
        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }
}

template<class T, int dim>
void SegmentMesh<T, dim>::WritePoly(const std::string& polyFileName) {
    //Open the file
    std::ofstream polyFile;
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
    polyFile << "POLYS\n";
    for(u_int32_t i = 0; i < indices.size(); ++i) {
        polyFile << i+1 << ": " << indices[i].first+1 << " " << indices[i].second+1 << "\n";
    }
//    polyFile << "POLYS\n1:";
//    for(u_int32_t i = 0; i < indices.size(); ++i) {
//        polyFile << " " << indices[i].first+1;
//    }
//    polyFile << " " << indices[indices.size()-1].second+1 << "\n";

    //Write the END section
    polyFile << "END";
    polyFile.close();
}

template class SegmentMesh<float,  2>;
template class SegmentMesh<float,  3>;
template class SegmentMesh<double, 2>;
template class SegmentMesh<double, 3>;