#include "polytope/DynamicStability.h"

#include <iostream>

using namespace equilib;

int main()
{
    init_library();
    vector_t x(2); x << 2, 0.5;
    vector_t y(2); y << 3, 0.5;
    vector_t fz(2); fz << 40., 40.;
    vector_t nu(2); nu << 0.5, 0.5;
    matrix_t cone = V_all(x, y, fz, nu);
//std::cout << "computed contact cone \n" << cone << std::endl;
    // assert all vectors are different
    for(int i = 0; i < 37; ++i)
    {
        for(int j = i+1; j < 37; ++j)
        {
            if((cone.block<6,1>(0,i) - cone.block<6,1>(0,j)).norm()
                    < 2 * std::numeric_limits<double>::epsilon())
            {
                std::cout << "found two equal vectors at col indexes " << i << " and " << j << std::endl;
                std::cout << cone.block<6,1>(0,i) << "\n \n" << cone.block<6,1>(0,j) << std::endl;
                return -1;
            }
        }
    }

    T_transform_t contact = Eigen::Matrix <value_type, 4, 4>::Identity();
    contact.block<3,1>(0,3) = vector3_t(1,1,1);
    T_transform_t contacts = Eigen::Matrix <value_type, 8, 4>::Zero();
    contacts.block<4,4>(0,0) = contact;
    contact.block<3,1>(0,3) = vector3_t(1,-1,1);
    contacts.block<4,4>(4,0) = contact;
    std::cout << "contacts " << contacts << std::endl;
//std::cout << "computed projection matrix \n" << A_stance(contacts) << std::endl;

    ProjectedCone* c = U_stance(contacts,nu,fz,x,y);

    ProjectedCone* c2 = U_stance(contacts.block<4,4>(0,0),
                                nu.head(1),fz.head(1),x.head(1),y.head(1));
    std::cout << "inequalities c\n" << c->HRepresentation() << std::endl;
    std::cout << "inequalities c2\n" << c2->HRepresentation() << std::endl;
    vector3_t gravity(0,0,-9.81);
    value_type mass(40);
    vector3_t  goodP(1,1,1);
    vector3_t  badP(-1,1,3);

    std::cout << "good position verified ? " << c2->IsValid(goodP,gravity,mass) << std::endl;
    std::cout << "bad position verified ? " << c2->IsValid(badP,gravity,mass) << std::endl;

    return 0;
}
