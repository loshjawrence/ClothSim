#include "Particles.h"


using T = double;
constexpr int dim = 3;

int main(int argc, char* argv[]) {
//    std::string file="test.pda";
    Particles<T,dim> particles(10);
    particles.RandInit(Particles<T, dim>::RAND_POS_SCALE, Particles<T,dim>::RAND_VEL_SCALE);

    const std::string file="particles";
    particles.WritePartio_RandomFrames(file);

    return 0;
}
