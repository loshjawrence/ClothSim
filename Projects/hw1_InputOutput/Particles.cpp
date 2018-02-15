//
// Created by josh on 2/1/18.
//

#include "Particles.h"
#include <iostream>
#include <sstream>
#include <fstream>

template<class T, int dim>
Particles<T, dim>::Particles() : pos() , vel()
{

}
template<class T, int dim>
Particles<T, dim>::Particles(const int N) : pos(N), vel(N), springForce(N), dampForce(N)
{

}

template<class T, int dim>
Particles<T, dim>::~Particles() {

}

template<class T, int dim>
void Particles<T, dim>::RandInit(const T& posScale, const T velScale) {
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        pos[i] = Eigen::Matrix<T, dim, 1>::Random() * posScale;
        vel[i] = Eigen::Matrix<T, dim, 1>::Random() * velScale;
        vel[i].normalize();//constant speed, random direction
    }
}

template<class T, int dim>
void Particles<T, dim>::SaveFixedPositions() {
    for(u_int32_t i = 0; i < fixed.size(); ++i) {
        fixedPos.push_back(pos[fixed[i]]);
    }
}

template<class T, int dim>
void Particles<T, dim>::RandInitVel(const T velScale) {
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        vel[i] = Eigen::Matrix<T, dim, 1>::Random() * velScale;
        vel[i].normalize();//constant speed, random direction
    }
}

template<class T, int dim>
void Particles<T, dim>::InitVert(const T& height, const T& spacing) {
    for(uint32_t i = 0; i < pos.size(); ++i) {
        pos[i].setZero();
        pos[i][1] = height - i*spacing;
    }
    InitVelZero();
}

template<class T, int dim>
void Particles<T, dim>::InitVelZero() {
    for(uint32_t i = 0; i < vel.size(); ++i) {
        vel[i].setZero();
    }
}

template<class T, int dim>
void Particles<T, dim>::UpdateRandVel(const T& deltaT, const T velScale) {
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        pos[i] += vel[i] * deltaT;
        vel[i] += Eigen::Matrix<T, dim, 1>::Random() * velScale;
        vel[i].normalize();//constant speed, random direction
    }
}

template<class T, int dim>
void Particles<T, dim>::CalcForces(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs) {
    //zero forces
    for(Eigen::Matrix<T,dim,1>& force : springForce) {
        force.setZero();
    }
    for(Eigen::Matrix<T,dim,1>& force : dampForce) {
        force.setZero();
    }

    const T epsilon = 0.0000001;
    //loop through connections and sum up damp and spring forces
    for(auto& tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
        //accumulate damping forces
        Eigen::Matrix<T,dim,1> n12 = pos[std::get<0>(tup)] - pos[std::get<1>(tup)];
        Eigen::Matrix<T,dim,1> v12 = vel[std::get<0>(tup)] - vel[std::get<1>(tup)];
        T l = n12.norm();//length
        n12.normalize();//or just divide by l
        Eigen::Matrix<T,dim,1> fd1 = -damp * n12.dot(v12) * n12;
        fd1 = fd1.norm() < epsilon ? fd1.setZero() : fd1;
        dampForce[std::get<0>(tup)] +=  fd1;
        dampForce[std::get<1>(tup)] += -fd1;

        //accumulate spring force
        const T restLen = std::get<2>(tup);
        Eigen::Matrix<T,dim,1> fs1 = -kHat * (l/restLen - 1) * n12;
        fs1 = fs1.norm() < epsilon ? fs1.setZero() : fs1;
        springForce[std::get<0>(tup)]  +=  fs1;
        springForce[std::get<1>(tup)]  += -fs1;
    }
}

template<class T, int dim>
void Particles<T, dim>::UpdateFE(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs ) {

    CalcForces(springs);

    const T invM = 1.0 / mass;
    Eigen::Matrix<T,dim,1> gravityForce; gravityForce.setZero();
    gravityForce[1] = mass*gravity;

    for(uint32_t i = 0; i < pos.size(); ++i) {
        pos[i] += dt * vel[i];
        const Eigen::Matrix<T, dim, 1> f = dampForce[i] + springForce[i] + gravityForce;
        vel[i] += dt * invM * f;
    }

    //reset any fixed particles
//    for(uint32_t i = 0; i < fixedPos.size(); ++i) {
//        pos[fixed[i]] = fixedPos[i];
//        vel[fixed[i]].setZero();
//    }

//    toggle fixed points back and forth
//    //set the first time function is called, uninit when program ends
    static T time = 0.0;
    static const T pi = 3.1415926535;
    static const T mag = 2.0;
    static const T freq = 0.1;
    static const uint32_t comp = 0;
    time += dt;
    for(uint32_t i = 0; i < fixedPos.size(); ++i) {
        pos[fixed[i]] = fixedPos[i];
        pos[fixed[i]][comp] += mag*sin(pi*i+freq*2*pi*time);
        vel[fixed[i]].setZero();
        vel[fixed[i]][comp] = ( pos[fixed[i]][comp] - (fixedPos[i][comp] + mag*sin(pi*i+freq*2*pi*(time-dt))) ) / dt;
    }
}

template<class T, int dim>
void Particles<T, dim>::WriteForces(const std::string& output) {
    std::ofstream outputFile;
    outputFile.open(output);
    std::stringstream outputLine;
    for(Eigen::Matrix<T,dim,1>& force : springForce) {
        outputLine << force;
    }
    for(Eigen::Matrix<T,dim,1>& force : dampForce) {
        outputLine << force;
    }
    outputLine << "\n";
    outputFile << outputLine.rdbuf();
    outputFile.close();
}

template<class T, int dim>
void Particles<T, dim>::WriteForces(std::ofstream& outputFile) {
    //outputFile is opened and closed outside of this scope
    std::stringstream outputLine;
    for(Eigen::Matrix<T,dim,1>& force : springForce) {
        for(int i = 0; i < dim; ++i) {
            outputLine << force[i] << " ";
        }
    }
    for(Eigen::Matrix<T,dim,1>& force : dampForce) {
        for(int i = 0; i < dim; ++i) {
            outputLine << force[i] << " ";
        }
    }
    outputLine << "\n";
    outputFile << outputLine.rdbuf();
}

template<class T, int dim>
void Particles<T, dim>::WritePartio(const std::string& particleFile) {
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, vH, mH;
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    posH = parts->addAttribute("position", Partio::VECTOR, dim);
    vH = parts->addAttribute("v", Partio::VECTOR, dim);
    for (u_int32_t i = 0; i < pos.size(); ++i) {
        int idx = parts->addParticle();
        float* m = parts->dataWrite<float>(mH, idx);
        float* p = parts->dataWrite<float>(posH, idx);
        float* v = parts->dataWrite<float>(vH, idx);
        m[0] = mass;
        for (int k = 0; k < dim; k++) {
            p[k] = pos[i][k];
        }
        for (int k = 0; k < dim; k++) {
            v[k] = vel[i][k];
        }
    }


    Partio::write(std::string(particleFile + ".bgeo").c_str(), *parts);
    parts->release();
}

template<class T, int dim>
void Particles<T, dim>::WritePartio_RandomFrames(const std::string &particleFile) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WritePartio(particleFile + std::to_string(i));
        UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }
}

template<class T, int dim>
void Particles<T, dim>::WritePoly(const std::string& polyFileName,
        const std::vector<std::tuple<uint32_t, uint32_t, T>>& indices)
{
    //Open the file
    std::ofstream polyFile;
    polyFile.open(polyFileName + ".poly");

    //Write the POINTS section
    polyFile << "POINTS\n";
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        polyFile << i+1 << ":";
        for(int k = 0; k < dim; ++k) {
            polyFile << " " << pos[i][k];
        }
        polyFile << "\n";
    }

    //Write the POLYS section, each entry is a string of connected points
    polyFile << "POLYS\n";
    for(u_int32_t i = 0; i < indices.size(); ++i) {
        const uint32_t one = std::get<0>(indices[i]);
        const uint32_t two = std::get<1>(indices[i]);
        polyFile << i+1 << ": " << one+1 << " " << two+1 << "\n";
    }

    //Write the END section
    polyFile << "END";
    polyFile.close();
}

//template class Particles<float,  2>;
//template class Particles<double, 2>;
template class Particles<float,  3>;
template class Particles<double, 3>;
