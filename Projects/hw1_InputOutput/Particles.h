//
// Created by josh on 2/1/18.
//

#ifndef CISPBA_PARTICLES_H
#define CISPBA_PARTICLES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Partio.h>
#include <vector>
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;



template<class T, int dim>
class Particles {

public:
    vector<Matrix<T, dim, 1>, aligned_allocator<Matrix<T,dim,1>>> pos;
    vector<Matrix<T, dim, 1>, aligned_allocator<Matrix<T,dim,1>>> vel;
    vector<T> mass;
    static constexpr T RAND_POS_SCALE = 0.5;
    static constexpr T RAND_VEL_SCALE = 0.5;

public:
    Particles();                    //create empty set of particles
    Particles(const int N);               //create N particles
    ~Particles();


    void RandInit(const T posScale, const T velScale);
    void UpdateRandVel(const T deltaT, const T velScale);

    void WritePartio(const string& particleFile);
    void WritePartio_RandomFrames(const string& particleFile);

};



#endif //CISPBA_PARTICLES_H
