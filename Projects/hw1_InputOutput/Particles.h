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
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> sprElastForce;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> sprDampForce;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> bendElastForce;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> bendDampForce;
    std::vector<uint32_t> fixed;
    std::vector<Eigen::Matrix<T,dim,1>, Eigen::aligned_allocator<Eigen::Matrix<T,dim,1>>> fixedPos;
    T kHat; //young's modulus, E
    T damp;//damping coefficient, b, gamma
    T kHatBend; //young's modulus, E
    T dampBend;//damping coefficient, b, gamma
    T c;
    T dt;
    T gravity;
    T mass;
    static constexpr T RAND_POS_SCALE = 0.5;
    static constexpr T RAND_VEL_SCALE = 0.5;

public:
    Particles();                    //create empty set of particles
    Particles(const int N);               //create N particles
    ~Particles();


    void InitVelZero();
    //initializes the segment in a vertical line up and down, first particle starting at height going down to last particle
    void InitVert(const T& height, const T& spacing);
    void SaveFixedPositions();
    void CalcForces(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
                    const std::vector<std::tuple<uint32_t, uint32_t, T>>& bendSprings);
    void UpdateFE(const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
                  const std::vector<std::tuple<uint32_t, uint32_t, T>>& bendSprings);

    void WritePartio(const std::string& particleFile);
    void printForcePosVel(std::ofstream& file, const int iter, const T& time,
                                            std::vector<std::tuple<uint32_t, uint32_t, T>>& springs);
    void WritePoly(const std::string& polyFileName, const std::vector<std::tuple<uint32_t, uint32_t, T>>& indices);
    void CheckSphere();
    void CheckGround();
    void AdjustFixedPoints();

//    void WriteForces(const std::string& output);
//    void WriteForces(std::ofstream& outputFile);
//    void WritePartio_RandomFrames(const std::string& particleFile);
//    void UpdateRandVel(const T& deltaT, const T velScale);
//    void RandInit(const T& posScale, const T velScale);
//    void RandInitVel(const T velScale);

};



#endif //CISPBA_PARTICLES_H
