#include "Particles.h"
#include "SegmentMesh.h"
#include "TriangleMesh.h"
#include <iostream>
#include <fstream>

using T = float;
constexpr int dim = 3;

void hw3Processor(const std::string& inputFile, const std::string& outputFile);
void WriteObj_Sim();
void WriteFrame(std::vector<TriangleMesh<T, dim>>& triMeshes, const uint32_t frame, std::ofstream& dbfile);

int main(int argc, char* argv[]) {
//    Particles<T,dim> particles(10);
//    particles.RandInit(Particles<T, dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);
//    particles.WritePartio_RandomFrames("FramesOutput/ParticleFrames/particles");

//    SegmentMesh<T, dim> segmentMesh(10);
//    segmentMesh.particles.RandInit(Particles<T,dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);
//    segmentMesh.WritePoly_RandomFrames("FramesOutput/SegmentMeshFrames/segmentMesh");

//    TriangleMesh<T, dim> TriangleMesh(10);
//    TriangleMesh.particles.RandInit(Particles<T,dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);
//    TriangleMesh.WriteObj_RandomFrames("FramesOutput/TriangleMeshFrames/triangleMesh");

    //JIGGLY CLOTH
//    const uint32_t size = 10;
//    TriangleMesh<T, dim> triangleMesh(size,size);
//    const uint32_t topleft = (size+1)*(size);
//    const uint32_t topright = topleft+size;
//    triangleMesh.particles.fixed = {topleft, topright};
//    triangleMesh.particles.SaveFixedPositions();
//    triangleMesh.WriteObj_Sim("FramesOutput/TriangleMeshFrames/triangleMesh", 240);

    //
    WriteObj_Sim();

    //JIGGLY SINGLE LINK SPRING
//    SegmentMesh<T, dim> segmentMesh(10);
//    segmentMesh.InitVert(10.0);
//    segmentMesh.particles.fixed = {0};
//    segmentMesh.particles.SaveFixedPositions();
//    segmentMesh.WritePoly_Sim("FramesOutput/SegmentMeshFrames/segmentMesh", 240);

    //HOMEWORK 3, no longer works b/c calc forces moved to particles
//    hw3Processor("hw3Input/test_input.txt", "hw3Input/mytest_output.txt");
//    hw3Processor("hw3Input/actual_input.txt", "hw3Input/myactual_output.txt");
    return 0;
}
void WriteObj_Sim() {
    //make array of cloth
    const uint32_t size = 10;
    const uint32_t numTriMeshes = 1;
    std::vector<TriangleMesh<T, dim>> triMeshes;
    for (uint32_t i = 0; i < numTriMeshes; ++i) {
        triMeshes.emplace_back(TriangleMesh<T, dim>(size, size, ORI::HORIZ, Eigen::Matrix<T,dim,1>(0,0,i*1.0)));
        const uint32_t topleft = (size + 1) * (size);
        const uint32_t topright = topleft + size;
        triMeshes[i].particles.fixed = {topleft, topright};
        triMeshes[i].particles.SaveFixedPositions();
    }

    //Open the debug file
    std::ofstream dbfile; dbfile.open("FramesOutput/Debug/triangle/trimesh");

    //setup the sphere
    struct sphere {
        T radius;
        Eigen::Matrix<T,dim,1> pos;
        void Update(const T globalTime) {
            const T pi = 3.1415926535;
            const T mag = 5.0;
            const T freq = 0.1;
            const uint32_t comp = 2;
            const T phase = pi / 4.0;
            pos[comp] = mag * sin(phase + freq * 2 * pi * globalTime);
        }
    };
    sphere theSphere = {2.0, Eigen::Matrix<T,dim,1>(5,5,0)};
    T globalTime = 0.0;

    //setup time vars
    const T dt = triMeshes[0].particles.dt;
    const T frameTime = 1.0/24.0;

    for(uint32_t frame = 1; frame <= 240; ++frame) {

        WriteFrame(triMeshes, frame, dbfile);

        //take tiny steps within the frame for phys sim
        for(T time = 0.0; time < frameTime; time += dt) {
            globalTime +=dt;

            //update ball
            theSphere.Update(globalTime);

            //update cloth
            for(uint32_t j = 0; j < numTriMeshes; ++j) {
                triMeshes[j].particles.UpdateFE(triMeshes[j].springs, triMeshes[j].bendSprings);
            }

            //collision checks
            for(uint32_t j = 0; j < numTriMeshes; ++j) {
                //check sphere
                triMeshes[j].particles.CheckSphere(theSphere.pos, theSphere.radius);

                //check other meshes
                for (uint32_t k = 0; k < numTriMeshes; ++k) {
                    if(j == k) {continue;}
                    triMeshes[j].particles.CheckOtherCloth(triMeshes[k].springs, triMeshes[k].particles);
                }//k
            }//j


        }//time
    }//frame
    dbfile.close();
}

void WriteFrame(std::vector<TriangleMesh<T, dim>>& triMeshes, const uint32_t frame, std::ofstream& dbfile) {
    for (uint32_t j = 0; j < triMeshes.size(); ++j) {
        const std::string frNumString = std::to_string(frame);
        const std::string triNumString = std::to_string(j);
        triMeshes[j].WriteObj("FramesOutput/TriangleMeshFrames/triangleMesh" + triNumString + "__" + frNumString);
        triMeshes[j].particles.WritePartio("FramesOutput/ParticleFrames/particles" + triNumString + "__" + frNumString);
        triMeshes[j].particles.WritePoly( "FramesOutput/SegmentMeshFrames/segmentMesh" + triNumString + "__" + frNumString, triMeshes[j].springs);
        //triMeshes[j].particles.printForcePosVel(dbfile, frame, 0.0,triMeshes[j].springs);
    }
}

//void hw3Processor(const std::string& input, const std::string& output) {
//    std::ifstream inputFile;
//    inputFile.open(input);
//
//    std::ofstream outputFile;
//    outputFile.open(output);
//
//    std::string inputLine;
//
//    std::vector<SegmentMesh<T,dim>> hw3networks;
//    if (inputFile.is_open()) {
//        while ( getline(inputFile, inputLine)) {
//            //lines have predefined format and points have predefined connections
//            //first 3 numbers are
//            //make a segment mesh with a specific connection based on diagram
//            SegmentMesh<T,dim> network; network.InitFromHW3Line(inputLine);
//            hw3networks.push_back(network);
//            hw3networks.back().CalcForces();
//        }
//    }
//
//    for(SegmentMesh<T,dim>& segmesh : hw3networks) {
//        //dump forces to file
//        segmesh.WriteForces(outputFile);
//    }
//
//    outputFile.close();
//    inputFile.close();
//}
