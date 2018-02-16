//
// Created by josh on 2/1/18.
//

#include "TriangleMesh.h"
#include <fstream>

template<class T, int dim>
TriangleMesh<T,dim>::TriangleMesh() : particles(), indices()
{

}

template<class T, int dim>
TriangleMesh<T, dim>::TriangleMesh(const int N) : particles(N), indices(3*(N-2))
{
    //This constructor assumes simple triangle fan construction
    AssignTriangleFanIndices();
}

template<class T, int dim>
TriangleMesh<T, dim>::TriangleMesh(const int X, const int Y) :
        particles((X+1)*(Y+1)), indices(3*2*(X*Y)), springs(4*X*Y + X + Y), bendSprings((Y+1)*(X-1)+(X+1)*(Y-1)), width(X), height(Y)
{
    //This constructor assumes a segmented quad sheet (piece of cloth), of X by Y quads
    particles.kHat = 10.0; //1.0, 10
    particles.damp = 1.0; //2.0
    particles.kHatBend = 10.0; //1.0
    particles.dampBend = 1.0; //2.0
    L0   = 1.0;   //m 1.0
    L0_d = L0 * sqrt(2.0);//m
    particles.c    = 0.01;   // kg/s 0.01
    particles.dt   = 0.002; //s 0.002, 0.012, 0.009
    particles.mass = 0.03;  // kg 0.001, 0.005
    particles.gravity = -9.81;  //m/s^2
    AssignGridMeshPositionsIndices();
    AssignSpringArray();
    particles.InitVelZero();
};

template<class T, int dim>
TriangleMesh<T, dim>::~TriangleMesh()
{

}

template<>
void TriangleMesh<float, 3>::AssignGridMeshPositionsIndices() {

    //vertices
    //Bot Left Corner: 0,0 ; Top Right Corner: width*L0, height*L0
    for (int32_t y = 0; y <= height; ++y) {
        for (int32_t x = 0; x <= width; ++x) {
            const Eigen::Matrix<float, 1,1> randZ = Eigen::Matrix<float,1,1>::Random()*0.05;
            particles.pos[x + (width+1)*y] = Eigen::Matrix<float, 3, 1>(x*L0, y*L0, randZ[0]);
//            std::cout << "\n\nPos(" << x <<"," << y << ") : ( " << x*L0 << ", " << y*L0 << ")";
        }
    }

    //indexing COUNTER CLOCKWISE
    //  -------
    //  |   2 |
    //  | *   |
    //  | 1 * |
    //  |_____|
    //
    uint32_t index = 0;
    uint32_t first, second, third;
    for (int32_t y = 0; y < height; ++y) {
        for (int32_t x = 0; x < width; ++x) {
            first = y*(width+1) + x;            //bot left
            second = first + 1;                 //bot right
            third = first + 1 + width;          //top left

//            std::cout << "\n\nf: " << first << " s: " << second << " t: " << third;
            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;

            first = third;                    //top left
            //second = second;                  //bot right
            third = first + 1;               //bot left

            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;
//            std::cout << "\nf: " << first << " s: " << second << " t: " << third;
        }
    }
}

template<>
void TriangleMesh<double, 3>::AssignGridMeshPositionsIndices() {
    //vertices
    //Bot Left Corner: 0,0 ; Top Right Corner: width*L0, height*L0
    for (int32_t y = 0; y <= height; ++y) {
        for (int32_t x = 0; x <= width; ++x) {
            const Eigen::Matrix<double, 1,1> randZ = Eigen::Matrix<double,1,1>::Random()*0.05;
            particles.pos[x + (width+1)*y] = Eigen::Matrix<double, 3, 1>(x*L0, y*L0, randZ[0]);
//            std::cout << "\n\nPos(" << x <<"," << y << ") : ( " << x*L0 << ", " << y*L0 << ")";
        }
    }

    //indexing COUNTER CLOCKWISE
    //  -------
    //  |   2 |
    //  | *   |
    //  | 1 * |
    //  |_____|
    //
    uint32_t index = 0;
    uint32_t first, second, third;
    for (int32_t y = 0; y < height; ++y) {
        for (int32_t x = 0; x < width; ++x) {
            first = y*(width+1) + x;            //bot left
            second = first + 1;                 //bot right
            third = first + 1 + width;          //top left

//            std::cout << "\n\nf: " << first << " s: " << second << " t: " << third;
            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;

            first = third;                    //top left
            //second = second;                  //bot right
            third = first + 1;               //bot left

            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;
//            std::cout << "\nf: " << first << " s: " << second << " t: " << third;
        }
    }

}

template<class T, int dim>
void TriangleMesh<T, dim>::AssignSpringArray() {

    //lower left corner is 0,0
    //go by quads, adding bottom, left and cross springs
    //
    //  |*   *
    //  | * *
    //  |  *
    //  | * *
    //  |*___*
    //
    uint32_t LL, UL, LR, UR;//lower left, upper left, lower right, upper right
    uint32_t index = 0;
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            LL = x + (width+1)*y;
            LR = LL+1;
            UL = LL + width + 1;
            UR = UL + 1;
//            std::cout << "\n\nquad" << index/4;
//            std::cout << "\nbot:"   << LL << "," << LR;
//            std::cout << "\nToUR:"  << LL << "," << UR;
//            std::cout << "\nleft:"  << LL << "," << UL;
//            std::cout << "\nToLR:"  << UL << "," << LR;
            springs[index++] = std::tuple<uint32_t, uint32_t, T>(LL, LR, L0);       //bottom spring
            springs[index++] = std::tuple<uint32_t, uint32_t, T>(LL, UR, L0_d);     //ToUR spring
            springs[index++] = std::tuple<uint32_t, uint32_t, T>(LL, UL, L0);       //left spring
            springs[index++] = std::tuple<uint32_t, uint32_t, T>(UL, LR, L0_d);     //ToLR spring
        }
    }

    //go back and add outside right and top springs (outside right and top edges of cloth)
    // _________
    //          |
    //          |
    //          |
    //          |
    //          |
    //          |
    //right edge
    for(int y = 0; y < height; ++y) {
        LR = (width+1)*y + width;
        UR = LR + width + 1;
//        std::cout << "\nrightedge" << y << ": " << LR << "," << UR;
        springs[index++] = std::tuple<uint32_t, uint32_t, T>(LR, UR, L0);
    }
    //top edge
    for(int x = 0; x < width; ++x) {
        UL = (width+1)*height + x;
        UR = UL + 1;
//        std::cout << "\ntopedge" << x << ": " << UL << "," << UR;
        springs[index++] = std::tuple<uint32_t, uint32_t, T>(UL, UR, L0);
    }


    //assign the horizontal bend springs
    // //(attempts to keep cloth locally flat, pushes adjacent faces away from each other)
    //they connect every other particle together(0,2), (1,3) etc
    index = 0;
    for(int y = 0; y <= height; ++y) {
        for(int x = 0; x < width-1; ++x) {
            LL = x + (width+1)*y;
            LR = LL + 2;
//            std::cout << "\nhorizBendSpring" << ": " << LL << "," << LR;
            bendSprings[index++] = std::tuple<uint32_t, uint32_t, T>(LL, LR, 2.0*L0);
        }
    }
    //vertical bend springs
    for(int x = 0; x <= width; ++x) {
        for(int y = 0; y < height-1; ++y) {
            LL = x + (width+1)*y;
            UL = LL + 2*(width+1);
//            std::cout << "\nvertiBendSpring" << ": " << LL << "," << UL;
            bendSprings[index++] = std::tuple<uint32_t, uint32_t, T>(LL, UL, 2.0*L0);
        }
    }


}

template<class T, int dim>
void TriangleMesh<T, dim>::AssignTriangleFanIndices() {
    const int numTris = particles.pos.size() - 2;
    int fanIndex = 1;
    for(int i = 0; i < numTris; ++i) {
        const int baseIndex = 3*i;
        indices[baseIndex    ] = 0;
        indices[baseIndex + 1] = fanIndex;
        indices[baseIndex + 2] = ++fanIndex;
    }
}

template<class T, int dim>
void TriangleMesh<T, dim>::WriteObj(const std::string& filename) {
    //Open the file
    std::ofstream file;
    file.open(filename + ".obj");

    //Write the v (vertex) section
    for(u_int32_t i = 0; i < particles.pos.size(); ++i) {
        file << "\nv";
        for(int k = 0; k < dim; ++k) {
            file << " " << particles.pos[i][k];
        }
    }

    //write the f (face) section
    const int numTris = indices.size() / 3;
    for(int i = 0; i < numTris; ++i) {
        const int baseIndex = 3*i;
        file << "\nf " << indices[baseIndex    ]+1;
        file << " "  << indices[baseIndex + 1]+1;
        file << " "  << indices[baseIndex + 2]+1;
    }

    file.close();

}

template<class T, int dim>
void TriangleMesh<T, dim>::Update(const T& frameTime, const uint32_t iter) {
    for(T time = 0.0; time < frameTime; time += particles.dt) {
//        CalcForces();
//        particles.UpdateFE(dt, mass, gravity);
        particles.UpdateFE(springs,bendSprings);
    }
}

template<class T, int dim>
void TriangleMesh<T,dim>::WriteObj_Sim(const std::string& polyFileName, const uint32_t numFrames) {
    const T frameTime = 1.f / 24.f;
    //Open the debug file
    std::ofstream file;
    file.open("FramesOutput/Debug/triangle/trimesh");

    for(uint32_t i = 1; i <= numFrames; ++i) {
        const std::string itoa = std::to_string(i);
        WriteObj(polyFileName + itoa);
        particles.WritePartio("FramesOutput/ParticleFrames/particles" + itoa);
        particles.WritePoly("FramesOutput/SegmentMeshFrames/segmentMesh" + itoa, springs);
//        particles.printForcePosVel(file, i, 0.0,springs);
        Update(frameTime, i);
    }
    file.close();
}


//template<class T, int dim>
//void TriangleMesh<T, dim>::WriteObj_RandomFrames(const std::string& objFileName, const uint32_t numFrames) {
//    const T time_step = 1.0/24.0;
//    for(uint32_t i = 1; i <= numFrames; ++i) {
//        WriteObj(objFileName + std::to_string(i));
//        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
//    }
//}

//template class TriangleMesh<float,  2>;
//template class TriangleMesh<double, 2>;
template class TriangleMesh<float,  3>;
template class TriangleMesh<double, 3>;
