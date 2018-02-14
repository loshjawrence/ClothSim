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


template<class T, int dim>
class Particles {

public:
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> pos;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> vel;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> springForce;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> dampForce;
    std::vector<T> mass;
    std::vector<uint32_t> fixed;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> fixedPos;
    static constexpr T RAND_POS_SCALE = 0.5;
    static constexpr T RAND_VEL_SCALE = 0.5;

public:
    Particles();                    //create empty set of particles
    Particles(const int N);               //create N particles
    ~Particles();


    void RandInit(const T& posScale, const T velScale);
    void RandInitVel(const T velScale);
    //initializes the segment in a vertical line up and down, first particle starting at height going down to last particle
    void InitVert(const T& height, const T& spacing);
    void SaveFixedPositions();
    void UpdateRandVel(const T& deltaT, const T velScale);
    void UpdateFE(const T& deltaT, const T& m, const T& g);
    void WriteForces(const std::string& output);
    void WriteForces(std::ofstream& outputFile);

    void WritePartio(const std::string& particleFile);
    void WritePartio_RandomFrames(const std::string& particleFile);

};



#endif //CISPBA_PARTICLES_H
