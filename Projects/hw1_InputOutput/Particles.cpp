//
// Created by josh on 2/1/18.
//

#include "Particles.h"


template<class T, int dim>
Particles<T, dim>::Particles() : pos() , vel(), mass()
{

}
template<class T, int dim>
Particles<T, dim>::Particles(int N) : pos(N), vel(N), mass(N)
{

}

template<class T, int dim>
Particles<T, dim>::~Particles() {

}

template<class T, int dim>
void Particles<T, dim>::RandInit(const T posScale, const T velScale) {
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        pos[i] = Matrix<T, dim, 1>::Random() * posScale;
        vel[i] = Matrix<T, dim, 1>::Random() * velScale;
    }
}

template<class T, int dim>
void Particles<T, dim>::UpdateRandVel(const T deltaT, const T velScale) {
    for(u_int32_t i = 0; i < pos.size(); ++i) {
        pos[i] += vel[i] * deltaT;
        vel[i] = Matrix<T, dim, 1>::Random() * velScale;
    }
}


template<class T, int dim>
void Particles<T, dim>::WritePartio(const string& particleFile) {
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
        m[0] = mass[i];
        for (int k = 0; k < dim; k++) {
            p[k] = pos[i][k];
        }
        for (int k = 0; k < dim; k++) {
            v[k] = vel[i][k];
        }
    }


    Partio::write(string(particleFile + ".bgeo").c_str(), *parts);
    parts->release();
}

template<class T, int dim>
void Particles<T, dim>::WritePartio_RandomFrames(const string &particleFile) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WritePartio(particleFile + to_string(i));
        UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }
}


template class Particles<float,  2>;
template class Particles<float,  3>;
template class Particles<double, 2>;
template class Particles<double, 3>;
