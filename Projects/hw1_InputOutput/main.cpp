#include "Particles.h"
#include "SegmentMesh.h"
#include "TriangleMesh.h"
#include <iostream>
#include <fstream>

using T = double;
constexpr int dim = 3;

void hw3Processor(const std::string& inputFile, const std::string& outputFile);

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

//    TriangleMesh<T, dim> TriangleMesh(2,2);
//    TriangleMesh.particles.RandInitVel(Particles<T,dim>::RAND_VEL_SCALE);
//    TriangleMesh.WriteObj_RandomFrames("FramesOutput/TriangleMeshFrames/triangleMesh");

    SegmentMesh<T, dim> segmentMesh(4);
    segmentMesh.InitVert(10.0);
    segmentMesh.particles.fixed = {0};
    segmentMesh.particles.SaveFixedPositions();
    segmentMesh.WritePoly_Sim("FramesOutput/SegmentMeshFrames/segmentMesh", 300);

//    hw3Processor("hw3Input/test_input.txt", "hw3Input/mytest_output.txt");
//    hw3Processor("hw3Input/actual_input.txt", "hw3Input/myactual_output.txt");


    return 0;
}

void hw3Processor(const std::string& input, const std::string& output) {
    std::ifstream inputFile;
    inputFile.open(input);

    std::ofstream outputFile;
    outputFile.open(output);

    std::string inputLine;

    std::vector<SegmentMesh<T,dim>> hw3networks;
    if (inputFile.is_open()) {
        while ( getline(inputFile, inputLine)) {
            //lines have predefined format and points have predefined connections
            //first 3 numbers are
            //make a segment mesh with a specific connection based on diagram
            SegmentMesh<T,dim> network; network.InitFromHW3Line(inputLine);
            hw3networks.push_back(network);
            hw3networks.back().CalcForces();
        }
    }

    for(SegmentMesh<T,dim>& segmesh : hw3networks) {
        //dump forces to file
        segmesh.WriteForces(outputFile);
    }

    outputFile.close();
    inputFile.close();
}
