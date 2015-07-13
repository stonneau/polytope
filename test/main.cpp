#include "polytope/DynamicStability.h"

#include <iostream>

using namespace equilib;

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
        ,fz(vector_t(2))
        ,nu(vector_t(2))
        ,contacts(Eigen::Matrix <value_type, 8, 4>::Zero())
    {}
    vector_t x;
    vector_t y;
    vector_t fz;
    vector_t nu;
    T_transform_t contacts;
};

void timeGeneration()
{
    int nbProblems = 100;
    std::vector<Problem> problems;
    problems.resize(nbProblems);
    for(int i=0; i < nbProblems; ++i)
    {
        Problem& p = problems[i];
        for(int cid =0; cid <2; ++cid)
        {
            p.x(cid) = unifRand(0.3, 2);
            p.y(cid) = unifRand(0.3, 2);
            p.fz(cid) = 400;
            p.nu(cid) = unifRand(0.1,1.);
            Eigen::Matrix <value_type, 3, 3> rot = Eigen::Matrix3d(
                     angle_axis_t(unifRand(M_PI), vector3_t::UnitX())
                    *angle_axis_t(unifRand(M_PI), vector3_t::UnitY())
                    *angle_axis_t(unifRand(M_PI), vector3_t::UnitZ()));
            vector3_t translation;
            for(int tid = 0; tid < 3; ++tid)
            {
                translation[tid] = unifRand(3.);
            }
            p.contacts.block<3,3>(cid*4,0) = rot;
            p.contacts.block<3,1>(cid*4,3) = translation;
        }
    }

    double time = get_time();
    int i = 0;
    for(std::vector<Problem>::const_iterator cit = problems.begin();
        cit != problems.end(); ++cit, ++i)
    {
        //ProjectedCone(V_all(cit->nu, cit->fz, cit->x, cit->y).transpose() * A_stance(cit->contacts));
        ProjectedCone(V_all(cit->nu.head(1), cit->fz.head(1), cit->x.head(1), cit->y.head(1)).transpose() * A_stance(cit->contacts.block<4,4>(0,0)));
        //delete cone;
    }
    time = get_time() - time;
    std::cout << "average time to generate two contact problem (one projection): \n"
              << time / (double)nbProblems << std::endl;
}


int creationTest()
{
    int ret = 0;
    vector_t x(2); x << 1, 1;
    vector_t y(2); y << 1, 1;
    vector_t fz(2); fz << 400., 400.;
    vector_t nu(2); nu << 1, 0.5;
    matrix_t cone = V_all(x, y, fz, nu);
    T_transform_t contact = Eigen::Matrix <value_type, 4, 4>::Identity();
    contact.block<3,1>(0,3) = vector3_t(0,0,0);
    T_transform_t contacts = Eigen::Matrix <value_type, 8, 4>::Zero();
    contacts.block<4,4>(0,0) = contact;
    contact.block<3,1>(0,3) = vector3_t(3,0,0);
    contacts.block<4,4>(4,0) = contact;
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
                ret = -1;
            }
        }
    }

    //const ProjectedCone* c = U_stance(contacts,nu,fz,x,y);
    const ProjectedCone* c2 = U_stance(contacts.block<4,4>(0,0),
                                nu.head(1),fz.head(1),x.head(1),y.head(1));
    vector3_t gravity(0,0,-9.81);
    value_type mass(40);
    vector3_t  goodP(0,0,3.);
    vector3_t  goodP2(0,0,0);
    vector3_t  badP(2,0,3.);

    std::cout << "good position verified ? " << c2->IsValid(goodP2,gravity,mass) << std::endl;
    std::cout << "good position verified ? " << c2->IsValid(goodP,gravity,mass) << std::endl;
    std::cout << "bad position verified ? " << c2->IsValid(badP,gravity,mass) << std::endl;
    //std::cout << "good position verified ? " << c->IsValid(badP,gravity,mass) << std::endl;
    return ret;
}

int main()
{
    srand (time(NULL));
    int ret = 0;
    init_library();
    ret+= creationTest();
    timeGeneration();
    timeGeneration();
    release_library();
    return ret;
}
