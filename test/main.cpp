#include "polytope/DynamicStability.h"

#include <iostream>

using namespace equilib;

int main()
{
    init_library();
    vector_t x(2); x << 2, 3;
    vector_t y(2); y << 4, 5;
    vector_t fz(2); fz << 8., 8.;
    vector_t nu(2); nu << 0.5, 1;
    matrix_t cone = friction_polytopes(x, y, fz, nu);

    //T_transform_t transform(8,4);
    // assert all vectors are different
    for(int i = 0; i < cone.cols(); ++i)
    {
        for(int j = i+1; j < cone.cols(); ++j)
        {
            if((cone.block<6,1>(0,i) - cone.block<6,1>(0,j)).norm()
                    < 2 * std::numeric_limits<double>::epsilon())
            {
                std::cout << "found two equal vectors at col indexes " << i << " and " << j << std::endl;
                return -1;
            }
        }
    }
    std::cout << "computed contact cone \n" << cone << std::endl;
    return 0;
}
