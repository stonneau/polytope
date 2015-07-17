#include "polytope/stability_margin.h"

#include "cdd/cddmp.h"
#include "cdd/setoper.h"
#include "cdd/cddtypes.h"
#include "cdd/cdd.h"

#include <iostream>
#include <memory>

using namespace polytope;

// defining timer
#ifdef WIN32
#include <windows.h>
double get_time()
{
    LARGE_INTEGER t, f;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&f);
    return (double)t.QuadPart/(double)f.QuadPart;
}

#else
#include <sys/time.h>
#include <sys/resource.h>

double get_time()
{
    struct timeval t;
    struct timezone tzp;
    gettimeofday(&t, &tzp);
    return t.tv_sec + t.tv_usec*1e-6;
}
#endif

double unifRand(value_type a, value_type b)
{
    return (b-a)*(rand() / value_type(RAND_MAX)) + a;
}

double unifRand(value_type b)
{
    return (b)*(rand() / value_type(RAND_MAX));
}

typedef Eigen::AngleAxis<value_type> angle_axis_t;

struct Problem
{
    Problem()
        :x(vector_t(2))
        ,y(vector_t(2))
        ,nu(vector_t(2))
        ,positions(vector_t(6))
        ,contacts(Eigen::Matrix <value_type, 6, 3>::Zero())
    {}
    vector_t x;
    vector_t y;
    vector_t nu;
    vector3_t xyz;
    vector_t positions;
    T_rotation_t contacts;
};

void timeGeneration()
{
    int nbProblems = 10000;
    std::vector<Problem> problems;
    problems.resize(nbProblems);
    for(int i=0; i < nbProblems; ++i)
    {
        Problem& p = problems[i];
        for(int cid =0; cid <2; ++cid)
        {
            p.x(cid) = unifRand(2, 2);
            p.y(cid) = unifRand(2, 2);
            p.nu(cid) = unifRand(0.1,1.);
            //p.nu(cid) = 1;
            p.xyz = vector3_t(unifRand(M_PI),unifRand(M_PI),unifRand(M_PI));
            Eigen::Matrix <value_type, 3, 3> rot = //-Eigen::Matrix <value_type, 3, 3>::Identity();
                    Eigen::Matrix <value_type, 3, 3>(
                      angle_axis_t(p.xyz[0], vector3_t::UnitX())
                    * angle_axis_t(p.xyz[1], vector3_t::UnitY())
                    * angle_axis_t(p.xyz[2], vector3_t::UnitZ()));
            vector3_t translation;
            for(int tid = 0; tid < 3; ++tid)
            {
                translation[tid] = 0; //unifRand(3.);
                translation[tid] = unifRand(1,3.);
            }
            p.contacts.block<3,3>(cid*3,0) = rot;
            p.positions.segment<3>(cid*3) = translation;
        }
    }

    double time = get_time();
    int saved(0); int failed(0);
    int i = 0;
    for(std::vector<Problem>::iterator cit = problems.begin();
        cit != problems.end(); ++cit, ++i)
    {
        //std::auto_ptr<const ProjectedCone> cone (U_stance(cit->contacts.block<3,3>(0,0), cit->positions.head(3), cit->nu.head(1), cit->x.head(1), cit->y.head(1)));
        std::auto_ptr<const ProjectedCone> cone (U_stance(cit->contacts, cit->positions, cit->nu, cit->x, cit->y));
        if(!cone.get())
        {
            /*std::cout << "rot x y z \n" << cit->xyz * 180 / M_PI << std::endl;
            std::cout << "positions " << cit->positions.head(3) << std::endl;
            std::cout << "friction " << cit->nu.head(1) << std::endl;
            cit->nu[0] = cit->nu[0] + 0.2;*/
            //const ProjectedCone* coneOffset = U_stance(cit->contacts.block<3,3>(0,0), cit->positions.head(3), cit->nu.head(1), cit->x.head(1), cit->y.head(1));
            cit->nu[0] = cit->nu[0] + 0.2;
            const ProjectedCone* coneOffset = U_stance(cit->contacts, cit->positions, cit->nu, cit->x, cit->y);
            if(coneOffset)
            {
                ++saved;
                delete coneOffset;
            }
            else
            {
                ++failed;
            }
            //dd_debug = false;
            //break;
        }
    }
    time = get_time() - time;
    std::cout << "average time to generate two contact problem (one projection): \n"
              << time / (double)nbProblems << std::endl;
    std::cout << "number of instable cones:" << (saved+failed) << std::endl;
    std::cout << "percentage of instable cones:" << (double)(saved+failed) * 100 / (double)nbProblems << "%" << std::endl;
    std::cout << "percentage of saved instable cones:" << (double)(saved) * 100 / (double)(saved+failed) << "%" << std::endl;
}

int creationTest()
{
    int ret = 0;
    vector_t x(2); x << 1, 1;
    vector_t y(2); y << 1, 1;
    vector_t nu(2); nu << 1, 0.5;
    T_rotation_t contact = Eigen::Matrix <value_type, 3,3>::Identity();
    T_rotation_t contacts = Eigen::Matrix <value_type, 6, 3>::Zero();
    contacts.block<3,3>(0,0) = contact;
    contacts.block<3,3>(3,0) = contact;
    vector_t position(6); position << 0,0,0,3,0,0;
    // assert all vectors are different
    matrix_t cone = V_all(x.head(1), y.head(1), nu.head(1)).transpose();
    std::cout << "cone \n" << cone.transpose() << std::endl;
    /*for(int i = 0; i < 37; ++i)
    {
        for(int j = i+1; j < 32; ++j)
        {
            if((cone.block<1,6>(i,0) - cone.block<1,6>(j,0)).norm()
                    < 2 * std::numeric_limits<double>::epsilon())
            {
                std::cout << "found two equal vectors at col indexes " << i << " and " << j << std::endl;
                std::cout << cone.block<1,6>(i,0) << "\n \n" << cone.block<1,6>(j,0) << std::endl;
                ret = -1;
            }
        }
    }*/

    const ProjectedCone* c = U_stance(contacts,position, nu,x,y);
    const ProjectedCone* c1 = U_stance(contacts.block<3,3>(0,0), position.head<3>(),
                                nu.head(1),x.head(1),y.head(1));
    const ProjectedCone* c2 = U_stance(contacts.block<3,3>(3,0), position.tail<3>(),
                                nu.tail(1),x.tail(1),y.tail(1));
    vector3_t gravity(0,0,-9.81);
    value_type mass(40);
    vector3_t  goodP(0,0,3.);
    vector3_t  goodP2(1,0,0);
    vector3_t  badP(1.5,0,3.); // valid only in double support

    std::cout << "good position verified ? " << c1->IsValid(goodP2,gravity,mass) << std::endl;
    std::cout << "good position verified ? " << c1->IsValid(goodP,gravity,mass) << std::endl;
    std::cout << "bad position verified ? " << c1->IsValid(badP,gravity,mass) << std::endl;
    std::cout << "bad position verified ? " << c2->IsValid(badP,gravity,mass) << std::endl;
    std::cout << "good position verified ? " << c->IsValid(badP,gravity,mass) << std::endl;
    delete c;
    delete c1;
    delete c2;
    return ret;
}

int main()
{
    srand (time(NULL));
    int ret = 0;
    init_library();
    ret+= creationTest();
    timeGeneration();
    Eigen::Matrix<value_type,2,1> nu, x, y;
    nu << 2, 4; x << 3, 6; y << 5, 10;

    /*Eigen::Matrix<value_type,1,1> nu, x, y;
    nu << 1; x << 1; y << 1;*/

    std::cout << V_all(nu,x,y) << std::endl;
    release_library();
    return ret;
}
