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
Particles<T, dim>::Particles(const int N) :
        pos(N), vel(N), sprElastForce(N), sprDampForce(N),
        bendElastForce(N), bendDampForce(N)
{

}

template<class T, int dim>
Particles<T, dim>::~Particles() {

}

template<class T, int dim>
void Particles<T, dim>::SaveFixedPositions() {
    for(u_int32_t i = 0; i < fixed.size(); ++i) {
        fixedPos.push_back(pos[fixed[i]]);
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
void Particles<T, dim>::CalcForces(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
                                   const std::vector<std::tuple<uint32_t, uint32_t, T>>& bendSprings)
{
    //zero forces
    for(Eigen::Matrix<T,dim,1>& force : sprElastForce) {
        force.setZero();
    }
    for(Eigen::Matrix<T,dim,1>& force : sprDampForce) {
        force.setZero();
    }
    for(Eigen::Matrix<T,dim,1>& force : bendElastForce) {
        force.setZero();
    }
    for(Eigen::Matrix<T,dim,1>& force : bendDampForce) {
        force.setZero();
    }

    const T epsilon = 0.0000001;
    //loop through connections and sum up damp and spring forces
    for(auto& tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
        Eigen::Matrix<T,dim,1> n12 = pos[std::get<0>(tup)] - pos[std::get<1>(tup)];
        Eigen::Matrix<T,dim,1> v12 = vel[std::get<0>(tup)] - vel[std::get<1>(tup)];
        const T l = n12.norm();//length
        n12.normalize();//or just divide by l
        const T restLen = std::get<2>(tup);

        //damping
        Eigen::Matrix<T,dim,1> fd1 = -damp * n12.dot(v12) * n12;
        fd1 = fd1.norm() < epsilon ? fd1.setZero() : fd1;
//        fd1 = l < restLen ? fd1.setZero() : fd1;//makes it super loose
        sprDampForce[std::get<0>(tup)] +=  fd1;
        sprDampForce[std::get<1>(tup)] += -fd1;

        //accumulate spring force
        Eigen::Matrix<T,dim,1> fs1 = -kHat * (l/restLen - 1) * n12;
        fs1 = fs1.norm() < epsilon ? fs1.setZero() : fs1;
        fs1 = l < restLen ? fs1.setZero() : fs1;
        sprElastForce[std::get<0>(tup)]  +=  fs1;
        sprElastForce[std::get<1>(tup)]  += -fs1;
    }

    //loop through connections and sum up damp and spring forces
    for(auto& tup: bendSprings) {//0: first particle index, 1: second particle index, 2: rest length
        Eigen::Matrix<T,dim,1> n12 = pos[std::get<0>(tup)] - pos[std::get<1>(tup)];
        Eigen::Matrix<T,dim,1> v12 = vel[std::get<0>(tup)] - vel[std::get<1>(tup)];
        const T l = n12.norm();//length
        n12.normalize();//or just divide by l
        const T restLen = std::get<2>(tup);

        //damping
        Eigen::Matrix<T,dim,1> fd1 = -dampBend * n12.dot(v12) * n12;
        fd1 = fd1.norm() < epsilon ? fd1.setZero() : fd1;
        fd1 = l > restLen ? fd1.setZero() : fd1;
        bendDampForce[std::get<0>(tup)] +=  fd1;
        bendDampForce[std::get<1>(tup)] += -fd1;

        //accumulate spring force
        Eigen::Matrix<T,dim,1> fs1 = -kHatBend * (l/restLen - 1) * n12;
        fs1 = fs1.norm() < epsilon ? fs1.setZero() : fs1;
        fs1 = l > restLen ? fs1.setZero() : fs1;
        bendElastForce[std::get<0>(tup)]  +=  fs1;
        bendElastForce[std::get<1>(tup)]  += -fs1;
    }
}

template<class T, int dim>
void Particles<T, dim>::UpdateFE(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
                                 const std::vector<std::tuple<uint32_t, uint32_t, T>>& bendSprings)
{

    CalcForces(springs, bendSprings);

    const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> posOrig = pos;
    const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> velOrig = vel;

    ForwardEuler_explicit();

    //adjust pos,vel based on collisions
    CheckSphere();
    CheckGround();
    CheckSelf(springs, posOrig, velOrig);

    //not sure about stretch constraints, when top fold flips to other side it bounces back oddly
    //might be better off with stiff k and low mass
    AdjustForMaxMinStretch(springs, posOrig, velOrig);
//    AdjustFixedPoints();
}

template<class T, int dim>
void Particles<T, dim>::ForwardEuler_explicit() {
    const T invM = 1.0 / mass;
    for(uint32_t i = 0; i < pos.size(); ++i) {
        pos[i] += dt * vel[i];
        Eigen::Matrix<T, dim, 1> f = sprDampForce[i] + sprElastForce[i] + bendDampForce[i] + bendElastForce[i];
        f[1] += mass*gravity;
        vel[i] += dt * invM * f;
    }
}


template<class T, int dim>
bool Particles<T, dim>::AreNeighbors(const uint32_t j, const uint32_t k,
                  const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs
)
{
    for (auto &tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
        const uint32_t one = std::get<0>(tup);
        const uint32_t two = std::get<1>(tup);
        //this does the same particle check as well
        if( (j == one || j == two) && (k == one || k == two)) {return true;}
    }
    return false;
}

template<class T, int dim>
void Particles<T, dim>::CheckSelf(
        const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
        const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>>& posOrig,
        const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>>& velOrig
)
{
    //must set radius based on stretch headRoom
    const T headRoom = Particles<T, dim>::HEADROOM;
    const T maxStretch = std::get<2>(springs[0]) + headRoom;
    const T radius = maxStretch * sqrt(2.0) * 0.25 + 0.01;//prevent worst case force field falls through middle of quad
    const T minDist = 2.0*radius;
    const uint32_t passes = 1;//adjusting one will affect the others

    const uint32_t size = pos.size();
    for(uint32_t i = 0; i < passes; ++i) {
        for(uint32_t j = 0; j < size; ++j) {//the main particle
            for(uint32_t k = 0; k < size; ++k) {//the one we check against
                //maybe pass in hash that has particle to adjacent particle list
//                if(AreNeighbors(j, k, springs)) { continue; }
                if(j == k) { continue; }//already have the minmaxstretch constraint so no need to ignore neighbors

                //check if sphere fields collide, if so push apart along vec that connects them
                //subtract off vel component (or zero out?)
                Eigen::Matrix<T, dim, 1> n12 = pos[j] - pos[k];
                const T l = n12.norm();//length
                if(l > minDist) { continue; }
                n12.normalize();//or just divide by l

                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.5*(minDist - l));
                pos[j] +=  adjustPos;
                pos[k] += -adjustPos;

                vel[j] -= (vel[j].dot(n12)*n12);
                vel[k] -= (vel[k].dot(n12)*n12);//still works out, no need to negate n12
            }//k
        }//j
    }//passes
}

template<class T, int dim>
void Particles<T, dim>::AdjustForMaxMinStretch(
        const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
        const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>>& posOrig,
        const std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>>& velOrig
)
{
    const T headRoom = Particles<T, dim>::HEADROOM;
    const T maxStretch = std::get<2>(springs[0]) + headRoom;
    const T minStretch = std::get<2>(springs[0]) - headRoom;
    const uint32_t passes = 1;//adjusting one will affect teh others

    for(uint32_t i = 0; i < passes; ++i) {
        for (auto &tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
            const uint32_t one = std::get<0>(tup);
            const uint32_t two = std::get<1>(tup);
            Eigen::Matrix<T, dim, 1> n12 = pos[one] - pos[two];
            const T l = n12.norm();//length
            n12.normalize();//or just divide by l
            const T restLen = std::get<2>(tup);
            const T stretch = l / restLen;
            if (stretch < minStretch) {
                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.5*(minStretch - stretch) * restLen);
                pos[one] +=  adjustPos;
                pos[two] += -adjustPos;
            } else if( stretch > maxStretch) {
                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.5*(stretch - maxStretch) * restLen);
                pos[one] += -adjustPos;
                pos[two] +=  adjustPos;
            }
        }//tup
    }//passes

    for(uint32_t i = 0; i < pos.size(); ++i) {
        //vel that could have been used to get from pre-euler update to this adjusted post euler update
        const Eigen::Matrix<T, dim, 1> actualVel = (pos[i] - posOrig[i]) / dt;
        //diff between vel that was used and vel that could have been used
        const Eigen::Matrix<T, dim, 1> diffVel = velOrig[i] - actualVel;
        //subtract off the excess from the current vel(post euler update)
        vel[i] -= diffVel;
    }
}

template<class T, int dim>
void Particles<T, dim>::CheckSphere() {
    //move sphere back and forth
    //set the first time function is called, uninit when program ends
    static T time = 0.0;
    const T radius = 2.0;
    const T epsilon = 0.001;
    const T pi = 3.1415926535;
    const T mag = 5.0;
    const T freq = 0.1;
    const uint32_t comp = 2;
    const T phase = pi/4.0;
    time += dt;

    Eigen::Matrix<T,dim,1> spherePos(5, 5, 0);

    spherePos[comp] = mag*sin(phase+freq*2*pi*time);
//    Eigen::Matrix<T,dim,1> sphereVel(0,0,0);
//    const Eigen::Matrix<T,dim,1> spherePosOrig(5, 5, 0);
//    sphereVel[comp] = ( spherePos[comp] - (spherePosOrig[comp] + mag*sin(phase+freq*2*pi*(time-dt))) ) / dt;

    for(uint32_t i = 0; i < pos.size(); ++i) {
        Eigen::Matrix<T,dim,1> n12 = pos[i] - spherePos;
        const T l = n12.norm();//length
        n12.normalize();//or just divide by l
        if (l < radius) {
            //for pos should push slightly beyond closest point on surface
            pos[i] += n12*(epsilon + radius-l);

            //for vel, dot vel with the normal of the point of intersection
            //this will give a neg value, mult this value with the normal to get penetration vector
            //subtract off this portion from the vel (similar to gramschmidt orthonormalization)
            //ensures no vel components below the plane of intersection(glides along surface)
            //vel vector will be shorter, sorta simulates energy loss on impact
            vel[i] -= (vel[i].dot(n12)*n12);
        }
    }
}

template<class T, int dim>
void Particles<T, dim>::CheckGround() {
    for(uint32_t i = 0; i < pos.size(); ++i) {
        if (pos[i][1] < 0.0) {
            //for vel should subtract off the dot product of the vel vec with the normal of the point of intersection with the object
            //for pos should push slightly beyond closest point on surface
            pos[i][1] = vel[i][1] = 0.00;
        }
    }
}


template<class T, int dim>
void Particles<T, dim>::AdjustFixedPoints() {
////    ////reset any fixed particles
    for(uint32_t i = 0; i < fixedPos.size(); ++i) {
        pos[fixed[i]] = fixedPos[i];
        vel[fixed[i]].setZero();
    }

//    toggle fixed points back and forth
//    //set the first time function is called, uninit when program ends
//    static T time = 0.0;
//    const T pi = 3.1415926535;
//    const T mag = 1.9;
//    const T freq = 0.2;
//    const T phase = pi/4.0;
//    const uint32_t comp = 0;
//    time += dt;
//    for(uint32_t i = 0; i < fixedPos.size(); ++i) {
//        pos[fixed[i]] = fixedPos[i];
//        pos[fixed[i]][comp] += i*mag*sin(phase+freq*2*pi*time) - i*mag;
//        vel[fixed[i]].setZero();
//        vel[fixed[i]][comp] = ( pos[fixed[i]][comp] - (fixedPos[i][comp] + i*mag*sin(phase+freq*2*pi*(time-dt)) -i*mag) ) / dt;
//    }
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
void Particles<T,dim>::printForcePosVel(std::ofstream& file, const int iter, const T& time,
    std::vector<std::tuple<uint32_t, uint32_t, T>>& springs)
{
//   https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
//    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

    file << "\n\niter" << iter << ", time" << time << " : ";
    for(uint32_t i = 0; i < pos.size(); ++i) {
        file << "\npos" << i << ": " << pos[i][0] << " " << pos[i][1] << " " << pos[i][2];
        file << "\nvel" << i << ": " << vel[i][0] << " " << vel[i][1] << " " << vel[i][2];
    }

    for(uint32_t i = 0; i < springs.size(); ++i) {
        const uint32_t one = std::get<0>(springs[i]);
        const uint32_t two = std::get<1>(springs[i]);

        file << "\nsprElast " << i << "( " << one << ", " << two << " ): "
             << sprElastForce[i][0] << " " << sprElastForce[i][1] << " " << sprElastForce[i][2];

        file << "\nsprDamp  "   << i << "( " << one << ", " << two << " ): "
             << sprDampForce[i][0]   << " " << sprDampForce[i][1]   << " " << sprDampForce[i][2];

        file << "\nbendElast"   << i << "( " << one << ", " << two << " ): "
             << bendElastForce[i][0]   << " " << bendElastForce[i][1]   << " " << bendElastForce[i][2];

        file << "\nbendDamp "   << i << "( " << one << ", " << two << " ): "
             << bendDampForce[i][0]   << " " << bendDampForce[i][1]   << " " << bendDampForce[i][2];
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
//
//template<class T, int dim>
//void Particles<T, dim>::WriteForces(const std::string& output) {
//    std::ofstream outputFile;
//    outputFile.open(output);
//    std::stringstream outputLine;
//    for(Eigen::Matrix<T,dim,1>& force : sprElastForce) {
//        outputLine << force;
//    }
//    for(Eigen::Matrix<T,dim,1>& force : sprDampForce) {
//        outputLine << force;
//    }
//    outputLine << "\n";
//    outputFile << outputLine.rdbuf();
//    outputFile.close();
//}
//
//template<class T, int dim>
//void Particles<T, dim>::WriteForces(std::ofstream& outputFile) {
//    //outputFile is opened and closed outside of this scope
//    std::stringstream outputLine;
//    for(Eigen::Matrix<T,dim,1>& force : sprElastForce) {
//        for(int i = 0; i < dim; ++i) {
//            outputLine << force[i] << " ";
//        }
//    }
//    for(Eigen::Matrix<T,dim,1>& force : sprDampForce) {
//        for(int i = 0; i < dim; ++i) {
//            outputLine << force[i] << " ";
//        }
//    }
//    outputLine << "\n";
//    outputFile << outputLine.rdbuf();
//}

//
//template<class T, int dim>
//void Particles<T, dim>::UpdateRandVel(const T& deltaT, const T velScale) {
//    for(u_int32_t i = 0; i < pos.size(); ++i) {
//        pos[i] += vel[i] * deltaT;
//        vel[i] += Eigen::Matrix<T, dim, 1>::Random() * velScale;
//        vel[i].normalize();//constant speed, random direction
//    }
//}

//template<class T, int dim>
//void Particles<T, dim>::RandInitVel(const T velScale) {
//    for(u_int32_t i = 0; i < pos.size(); ++i) {
//        vel[i] = Eigen::Matrix<T, dim, 1>::Random() * velScale;
//        vel[i].normalize();//constant speed, random direction
//    }
//}

//template<class T, int dim>
//void Particles<T, dim>::RandInit(const T& posScale, const T velScale) {
//    for(u_int32_t i = 0; i < pos.size(); ++i) {
//        pos[i] = Eigen::Matrix<T, dim, 1>::Random() * posScale;
//        vel[i] = Eigen::Matrix<T, dim, 1>::Random() * velScale;
//        vel[i].normalize();//constant speed, random direction
//    }
//}

//
//template<class T, int dim>
//void Particles<T, dim>::WritePartio_RandomFrames(const std::string &particleFile) {
//    const int numFrames = 120;
//    const T time_step = 1.0/24.0;
//    for(int i = 1; i <= numFrames; ++i) {
//        WritePartio(particleFile + std::to_string(i));
//        UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
//    }
//}
//template class Particles<float,  2>;
//template class Particles<double, 2>;
template class Particles<float,  3>;
template class Particles<double, 3>;
