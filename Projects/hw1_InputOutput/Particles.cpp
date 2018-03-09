//
// Created by josh on 2/1/18.
//

#include "Particles.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>//MINRES

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
//    BackwardEuler_implicit(springs,bendSprings);

    CheckSelf(springs);

    //not sure about stretch constraints, when top fold flips to other side it bounces back oddly
    //might be better off with stiff k and low mass
    AdjustForMaxMinStretch(springs, posOrig, velOrig);
    CheckGround();
    AdjustFixedPoints();
}

void print3x3Mat(const Eigen::Matrix<float, 3, 3>& mat, const std::string& descr) {
    std::cout << "\n\n" << descr << ":\n" << mat << "\n-----------------------";
}

void printSparseMat(const Eigen::SparseMatrix<float>& mat, const std::string& descr) {
    Eigen::IOFormat CleanFmt(6, 0, ", ", "\n", "[", "]");
    std::cout << "\n\n" << descr << ":\n" << Eigen::MatrixXf(mat).format(CleanFmt) << "\n-----------------------";
}

template<>
void Particles<float, 3>::BackwardEuler_implicit(
        const std::vector<std::tuple<uint32_t, uint32_t, float>>& springs,
        const std::vector<std::tuple<uint32_t, uint32_t, float>>& bendSprings
) {
    //build rhs[3nx1], init to 0
    //Terms for rhs: M*V@n, M*g, F_elastic@x@n
    using T = float;
    const uint32_t dim = 3;
    const uint32_t numParticles = pos.size();
    const T invdt = 1.0 / dt;
    const T invdtdt = invdt*invdt;

    //mass*gravity vec
    Eigen::Matrix<T, dim, 1> mg; mg.setZero(); mg(1) = mass*gravity;

    //rhs vec
    Eigen::MatrixXf rhs(dim*numParticles,1); rhs.setZero();
    for(uint32_t i = 0; i < numParticles; ++i) {
        const Eigen::Matrix<T, dim, 1> sum = invdt * mass * vel[i] + mg + sprElastForce[i] + bendElastForce[i];
        for(uint32_t j = 0; j < dim; ++j) {//populate elem by elem
            rhs(dim*i+j)= sum(j);
        }
    }

    //NOTE: FULL equation that minres is solving Ax = b, where A is A, b is the rhs var and x is sigma_x in the notes
    // A * sigma_x = rhs:
    //(1/dt^2)*M -(1/dt)*G - K) * sigma_x = 1/(dt)*M*v@n + (M*g + F_elast@x@n)

    //Build SPmat A[3nx3n], init to 0, each 3x3 chunk represents and i<-->j spring interaction/entry
    Eigen::SparseMatrix<T> A(dim*numParticles, dim*numParticles); A.setIdentity();

    //add M terms to A(each entry diagonal of 3x3 chunk is the mass*(1/dt^2)), only entries i<-->j are populated where i==j (the main diag)
    A *= mass*invdtdt;

    //add G terms to A(four 3x3 entries for interaction between i<-->j equal b*n*nT, also i,i and j,j chunks are negative
    //b is the damp coeff. n*nT will result in 3x3 matrix(3x1 vec mult by 1x3 vec). must be mult by -(1/dt)*damp.
    //add K terms to A(four 3x3 entries similar to G (i,i and j,j entries are negative) but the 3x3 chunk(k_s) is defined:
    //k_s = E*(1/L0 - 1/L)(I - n*nT) + n*nT*E/L0; where E is youngs modulus, I is 3x3 identity L0 is rest length for
    //spring between nodes i and j, L is the distance between nodes i and j. n*nT will result in 3x3 matrix(3x1 vec mult by 1x3 vec
    Eigen::Matrix<T, dim, dim> I; I.setIdentity();
    for(auto& spring : springs) {
        //one needs to be < two for A update loop to work
        uint32_t one  = std::get<0>(spring);
        uint32_t two  = std::get<1>(spring);
        if(one > two) {
            const uint32_t tmp = one;
            one = two;
            two = tmp;
        }

        const T restLen     = std::get<2>(spring);
        const T E_over_L0 = kHat / restLen;

        Eigen::Matrix<T, dim, 1> n12 = pos[one] - pos[two];
        const T E_over_L = kHat / n12.norm();
        n12.normalize();

        const Eigen::Matrix<T, dim, dim> nnT = n12 * n12.transpose();
        const Eigen::Matrix<T, dim, dim> G = -invdt * damp * nnT;
        const Eigen::Matrix<T, dim, dim> K = -1.0 * ((E_over_L0 - E_over_L)*(I - nnT) + (E_over_L0*nnT));
        const Eigen::Matrix<T, dim, dim> GplusK = G + K;

        //update the four 3x3 chunks of A with G and K
        for(uint32_t i = one; i <= two; i += two-one) {//row
            for(uint32_t j = one; j <=two; j += two-one) {//col
                //update A elem by elem
                const T factor = i == j ? -1.0 : 1.0;//3x3 chunks on the diag are negative
                for(uint32_t row = 0; row < dim; ++row) {
                    for(uint32_t col = 0; col < dim; ++col) {
                        //probably best to build a list of 'triplets' and use setFromTriplets(iterBegin,iterEnd)
                        //where a triplet is a tuple<uint32_t, uint32_t, T> : row,col,value
                        A.coeffRef(i*dim+row, j*dim+col) += factor * GplusK(row,col);
                    }//col
                }//row
            }//j
        }//i
//        printSparseMat(A, "A afer " + std::to_string(one) + "," + std::to_string(two));

    }//springs

    for(auto& spring : bendSprings) {
        //one needs to be < two for A update loop to work
        uint32_t one  = std::get<0>(spring);
        uint32_t two  = std::get<1>(spring);
        if(one > two) {
            const uint32_t tmp = one;
            one = two;
            two = tmp;
        }

        const T restLen     = std::get<2>(spring);
        const T E_over_L0 = kHatBend / restLen;

        Eigen::Matrix<T, dim, 1> n12 = pos[one] - pos[two];
        const T E_over_L = kHatBend / n12.norm();
        n12.normalize();

        const Eigen::Matrix<T, dim, dim> nnT = n12 * n12.transpose();
        const Eigen::Matrix<T, dim, dim> G = -invdt * dampBend * nnT;
        const Eigen::Matrix<T, dim, dim> K = -1.0 * ((E_over_L0 - E_over_L)*(I - nnT) + (E_over_L0*nnT));
        const Eigen::Matrix<T, dim, dim> GplusK = G + K;

        //update the four 3x3 chunks of A with G and K
        for(uint32_t i = one; i <= two; i += two-one) {//row
            for(uint32_t j = one; j <=two; j += two-one) {//col
                //update A elem by elem
                const T factor = i == j ? -1.0 : 1.0;//3x3 chunks on the diag are negative
                for(uint32_t row = 0; row < dim; ++row) {
                    for(uint32_t col = 0; col < dim; ++col) {
                        //probably best to build a list of 'triplets' and use setFromTriplets(iterBegin,iterEnd)
                        //where a triplet is a tuple<uint32_t, uint32_t, T> : row,col,value
                        A.coeffRef(i*dim+row, j*dim+col) += factor * GplusK(row,col);
                    }//col
                }//row
            }//j
        }//i
    }//spring


    //call sigma_x = MINRES(A,rhs);// sigma_x is 3nx1
    //x@n+1 = x@n + sigma_x;
    //v@n+1 = sigma_x / dt;
    Eigen::MINRES< Eigen::SparseMatrix<T> > mr(A);
    const Eigen::MatrixXf sigma_x = mr.solve(rhs);
    for(uint32_t i = 0; i < numParticles; ++i) {
        Eigen::Matrix<T, dim, 1> sigma_x_respective;
        for(uint32_t k = 0; k < dim; ++k) { sigma_x_respective(k) = sigma_x(i*dim+k); }
        pos[i] += sigma_x_respective;
        vel[i] = sigma_x_respective * invdt;
    }
}


template<>
void Particles<double, 3>::BackwardEuler_implicit(
        const std::vector<std::tuple<uint32_t, uint32_t, double>>& springs,
        const std::vector<std::tuple<uint32_t, uint32_t, double>>& bendSprings
)
{
}

template<class T, int dim>
void Particles<T, dim>::ForwardEuler_explicit() {
    const T invM = 1.0 / mass;
    const T mg = mass*gravity;
    for(uint32_t i = 0; i < pos.size(); ++i) {
        pos[i] += dt * vel[i];
        Eigen::Matrix<T, dim, 1> f = sprDampForce[i] + sprElastForce[i] + bendDampForce[i] + bendElastForce[i];
        f[1] += mg;
        vel[i] += dt * invM * f;
    }
}



template<class T, int dim>
void Particles<T, dim>::CheckSelf( const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs ) {
    //must set radius based on stretch headRoom
    const T headRoom = Particles<T, dim>::HEADROOM;
    const T L0_d_MaxStretch = std::get<2>(springs[1]) * (1.0 + headRoom);//NOT A RATIO, ACTUAL dist
    const T radius = L0_d_MaxStretch * 0.25;//prevent worst case force field falls through middle of quad
    const T minDist = radius*2.1;//add a little to prevent issues (still can fall through due to iterative pass adjustments not being enough)
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
void Particles<T, dim>::CheckOtherCloth(
        const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs,
        Particles<T,dim>& otherCloth
)
{
    //must set radius based on stretch headRoom
    const T headRoom = Particles<T, dim>::HEADROOM;
    const T L0_d_MaxStretch = std::get<2>(springs[1]) * (1.0 + headRoom);//NOT A RATIO, ACTUAL dist
    const T radius = L0_d_MaxStretch * 0.25;//prevent worst case force field falls through middle of quad
    const T minDist = radius*2.1;//add a little to prevent issues (still can fall through due to iterative pass adjustments not being enough)
    const uint32_t passes = 1;//adjusting one will affect the others

    const uint32_t size = pos.size();
    for(uint32_t i = 0; i < passes; ++i) {
        for(uint32_t j = 0; j < size; ++j) {//the main particle

            for(uint32_t k = 0; k < size; ++k) {//the one in the otherCloth to check against
                //check if sphere fields collide, if so push apart along vec that connects them
                //subtract off vel component (or zero out?)
                Eigen::Matrix<T, dim, 1> n12 = pos[j] - otherCloth.pos[k];
                const T l = n12.norm();//length
                if(l > minDist) { continue; }
                n12.normalize();//or just divide by l

                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.5*(minDist - l));
                pos[j] +=  adjustPos;
                otherCloth.pos[k] += -adjustPos;
                vel[j] -= (vel[j].dot(n12)*n12);
                otherCloth.vel[k] -= (otherCloth.vel[k].dot(n12)*n12);//still works out, no need to negate n12


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
    const T headRoom = Particles<T, dim>::HEADROOM;//this is a percentage of the rest length
    const T maxStretch = 1.0+headRoom;//THIS IS A RATIO
    const T minStretch = 1.0-headRoom;
    const uint32_t passes = 6;//adjusting one will affect teh others, need more for large dt steps on backward euler

    //was added in each of the stretch if statements but has issues
//                vel[one] -= (vel[one].dot(n12)*n12);//HAS BAD ARTIFACTS, PROBABLY OVER CORRECTS (MANY SPRINGS ATTACHED TO PARTICLE)
//                vel[two] -= (vel[two].dot(n12)*n12);

    for(uint32_t i = 0; i < passes; ++i) {
        for (auto &tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
            const uint32_t one = std::get<0>(tup);
            const uint32_t two = std::get<1>(tup);
            Eigen::Matrix<T, dim, 1> n12 = pos[one] - pos[two];
            const T l = n12.norm();//length
            n12.normalize();//or just divide by l
            const T restLen = std::get<2>(tup);
            const T stretch = l / restLen;
            bool oneFixed = false;
            bool twoFixed = false;
            for(uint32_t i = 0; i < fixedPos.size(); ++i) {
                if(fixed[i] == one) {oneFixed = true;}
                if(fixed[i] == two) {twoFixed = true;}
            }
            if (stretch < minStretch) {
                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.50*(minStretch - stretch) * restLen);

                if(!oneFixed && !twoFixed) {
                    pos[one] +=  adjustPos;
                    pos[two] += -adjustPos;
//                    vel[one] -= (vel[one].dot(n12)*n12);
//                    vel[two] -= (vel[two].dot(n12)*n12);
                } else if (oneFixed && !twoFixed) {
                    pos[two] += -2.0 * adjustPos;
//                    vel[two] -= (vel[two].dot(n12)*n12);
                } else if (!oneFixed && twoFixed) {
                    pos[one] +=  2.0*adjustPos;
//                    vel[one] -= (vel[one].dot(n12)*n12);
                } else if (oneFixed && twoFixed) {
                    //do nothing
                }

//                //THE MANY UPDATES CAUSES EXPLOSIONS (prob because posOrig isnt getting updated in the 2 lines above)
//                //vel that could have been used to get from pre-euler update to this adjusted post euler update
//                Eigen::Matrix<T, dim, 1> actualVel = (pos[one] - posOrig[one]) / dt;
//                //diff between vel that was used and vel that could have been used
//                Eigen::Matrix<T, dim, 1> diffVel = velOrig[one] - actualVel;
//                //subtract off the excess from the current vel(post euler update)
//                vel[one] -= diffVel;
//                actualVel = (pos[two] - posOrig[two]) / dt;
//                diffVel = velOrig[two] - actualVel;
//                vel[two] -= diffVel;
            } else if( stretch > maxStretch) {
                const Eigen::Matrix<T, dim, 1> adjustPos = n12 * (0.50*(stretch - maxStretch) * restLen);
                if(!oneFixed && !twoFixed) {
                    pos[one] += -adjustPos;
                    pos[two] +=  adjustPos;
//                    vel[one] -= (vel[one].dot(n12)*n12);
//                    vel[two] -= (vel[two].dot(n12)*n12);
                } else if (oneFixed && !twoFixed) {
                    pos[two] += 2.0 * adjustPos;
//                    vel[two] -= (vel[two].dot(n12)*n12);
                } else if (!oneFixed && twoFixed) {
                    pos[one] += -2.0*adjustPos;
//                    vel[one] -= (vel[one].dot(n12)*n12);
                } else if (oneFixed && twoFixed) {
                    //no nothing
                }

//                //vel that could have been used to get from pre-euler update to this adjusted post euler update
//                Eigen::Matrix<T, dim, 1> actualVel = (pos[one] - posOrig[one]) / dt;
//                //diff between vel that was used and vel that could have been used
//                Eigen::Matrix<T, dim, 1> diffVel = velOrig[one] - actualVel;
//                //subtract off the excess from the current vel(post euler update)
//                vel[one] -= diffVel;
//                actualVel = (pos[two] - posOrig[two]) / dt;
//                diffVel = velOrig[two] - actualVel;
//                vel[two] -= diffVel;
            }
        }//tup
    }//passes

    //can cause oscillations, vs the gram schmidt method in the two ifs above but they correct as well as this
    //causes reflections when hitting the ground if performed after ground update (or any other objects)
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
void Particles<T, dim>::CheckSphere(const Eigen::Matrix<T,dim,1>& spherePos, const T radius) {
    for(uint32_t i = 0; i < pos.size(); ++i) {
        Eigen::Matrix<T,dim,1> n12 = pos[i] - spherePos;
        const T l = n12.norm();//length
        n12.normalize();//or just divide by l
        if (l < radius) {
            //for pos should push slightly beyond closest point on surface
            pos[i] += n12*(0.001 + radius-l);

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

//template<class T, int dim>
//bool Particles<T, dim>::AreNeighbors(const uint32_t j, const uint32_t k,
//                  const std::vector<std::tuple<uint32_t, uint32_t, T>>& springs
//)
//{
//    for (auto &tup: springs) {//0: first particle index, 1: second particle index, 2: rest length
//        const uint32_t one = std::get<0>(tup);
//        const uint32_t two = std::get<1>(tup);
//        //this does the same particle check as well
//        if( (j == one || j == two) && (k == one || k == two)) {return true;}
//    }
//    return false;
//}

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
