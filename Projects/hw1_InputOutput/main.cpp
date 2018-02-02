#include "Particles.h"
#include "SegmentMesh.h"
#include "TriangleMesh.h"


using T = double;
constexpr int dim = 3;

int main(int argc, char* argv[]) {
//    Particles<T,dim> particles(10);
//    particles.RandInit(Particles<T, dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);

//    const std::string file="FramesOutput/ParticleFrames/particles";
//    particles.WritePartio_RandomFrames(file);

//    SegmentMesh<T, dim> segmentMesh(10);
//    segmentMesh.particles.RandInit(Particles<T,dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);
//    segmentMesh.WritePoly_RandomFrames("FramesOutput/SegmentMeshFrames/segmentMesh");

    TriangleMesh<T, dim> TriangleMesh(10);
    TriangleMesh.particles.RandInit(Particles<T,dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);
    TriangleMesh.WriteObj_RandomFrames("FramesOutput/TriangleMeshFrames/triangleMesh");

    return 0;
}
