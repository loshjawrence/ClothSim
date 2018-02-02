#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <vector>
#include <string>
#include <iostream>

using T = double;
constexpr int dim = 3;


int main(int argc, char* argv[])
{
    std::vector<Eigen::Matrix<T,dim,1>> hello(1);
    Eigen::Matrix<T,dim,1> x;
    x << (T).1, (T).2, (T).3;
    hello[0] << (T).1, (T).2, (T).3;
    std::cout << x.transpose() << std::endl;
    std::cout << hello[0].transpose() << std::endl;

    return 0;
}
