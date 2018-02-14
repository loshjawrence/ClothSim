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
        particles((X+1)*(Y+1)), indices(3*2*(X*Y)), springs(4*X*Y + X + Y), width(X), height(Y)
{
    //This constructor assumes a segmented quad sheet (piece of cloth), of X by Y quads
    mass = 0.001;  // kg
    kHat = 1.0;
    c    = 0.01;   // kg/s
    grav = -9.81;  //m/s^2
    dt   = 0.002; //s
    L0   = 1.0;   //m
    L0_d = L0 * sqrt(2.0);//m
    damp = 2.0;
    AssignGridMeshPositionsIndices();
    AssignSpringArray();
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
            particles.pos[x + (width+1)*y] = Eigen::Matrix<float, 3, 1>(x*L0, y*L0, 0.0);
            std::cout << "\n\nPos(" << x <<"," << y << ") : ( " << x*L0 << ", " << y*L0 << ")";
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

            std::cout << "\n\nf: " << first << " s: " << second << " t: " << third;
            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;

            first = third;                    //top left
            //second = second;                  //bot right
            third = first + 1;               //bot left

            indices[index++] = first;
            indices[index++] = second;
            indices[index++] = third;
            std::cout << "\nf: " << first << " s: " << second << " t: " << third;
        }
    }
}

template<>
void TriangleMesh<double, 3>::AssignGridMeshPositionsIndices() {
    //vertices
    //Bot Left Corner: 0,0 ; Top Right Corner: width*L0, height*L0
    for (int32_t y = 0; y <= height; ++y) {
        for (int32_t x = 0; x <= width; ++x) {
            particles.pos[x + (width+1)*y] = Eigen::Matrix<double, 3, 1>(x*L0, y*L0, 0.0);
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
//
//template<>
//void TriangleMesh<float, 2>::AssignGridMeshPositionsIndices() {
//}
//
//template<>
//void TriangleMesh<double, 2>::AssignGridMeshPositionsIndices() {
//}

template<class T, int dim>
void TriangleMesh<T, dim>::AssignSpringArray() {

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
            std::cout << "\n\nquad" << index/4;
            std::cout << "\nbot:"   << LL << "," << LR;
            std::cout << "\nToUR:"  << LL << "," << UR;
            std::cout << "\nleft:"  << LL << "," << UL;
            std::cout << "\nToLR:"  << UL << "," << LR;
            springs[index++] = std::pair<uint32_t, uint32_t>(LL, LR);   //bottom spring
            springs[index++] = std::pair<uint32_t, uint32_t>(LL, UR);   //ToUR spring
            springs[index++] = std::pair<uint32_t, uint32_t>(LL, UL);   //left spring
            springs[index++] = std::pair<uint32_t, uint32_t>(UL, LR);   //ToLR spring
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
        std::cout << "\nrightedge" << y << ": " << LR << "," << UR;
        springs[index++] = std::pair<uint32_t, uint32_t>(LR, UR);
    }
    //top edge
    for(int x = 0; x < width; ++x) {
        UL = (width+1)*height + x;
        UR = UL + 1;
        std::cout << "\ntopedge" << x << ": " << UL << "," << UR;
        springs[index++] = std::pair<uint32_t, uint32_t>(UL, UR);
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
void TriangleMesh<T, dim>::WriteObj_RandomFrames(const std::string& objFileName) {
    const int numFrames = 120;
    const T time_step = 1.0/24.0;
    for(int i = 1; i <= numFrames; ++i) {
        WriteObj(objFileName + std::to_string(i));
        particles.UpdateRandVel(time_step, Particles<T,dim>::RAND_VEL_SCALE);
    }

}


//template class TriangleMesh<float,  2>;
//template class TriangleMesh<double, 2>;
template class TriangleMesh<float,  3>;
template class TriangleMesh<double, 3>;
