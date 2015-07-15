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
        ,fz(vector_t(2))
        ,nu(vector_t(2))
        ,positions(vector_t(6))
        ,contacts(Eigen::Matrix <value_type, 6, 3>::Zero())
    {}
    vector_t x;
    vector_t y;
    vector_t fz;
    vector_t nu;
    vector_t positions;
    T_rotation_t contacts;
};

void timeGeneration()
{
    int nbProblems = 1000;
    std::vector<Problem> problems;
    problems.resize(nbProblems);
    for(int i=0; i < nbProblems; ++i)
    {
        Problem& p = problems[i];
        for(int cid =0; cid <2; ++cid)
        {
            p.x(cid) = unifRand(2, 2);
            p.y(cid) = unifRand(2, 2);
            p.fz(cid) = 1;
            p.nu(cid) = unifRand(0.1,1.);
            Eigen::Matrix <value_type, 3, 3> rot = Eigen::Matrix <value_type, 3, 3>(
                      angle_axis_t(unifRand(M_PI), vector3_t::UnitX())
                    * angle_axis_t(unifRand(M_PI), vector3_t::UnitY())
                    * angle_axis_t(unifRand(M_PI), vector3_t::UnitZ()));
            vector3_t translation;
            for(int tid = 0; tid < 3; ++tid)
            {
                translation[tid] = unifRand(3.);
            }
            p.contacts.block<3,3>(cid*3,0) = rot;
            p.positions.segment<3>(cid*3) = translation;
        }
    }

    double time = get_time();
    int i = 0;
    for(std::vector<Problem>::const_iterator cit = problems.begin();
        cit != problems.end(); ++cit, ++i)
    {
        //std::auto_ptr<const ProjectedCone> cone(U_stance(cit->contacts.block<3,3>(0,0), cit->positions.head(3), cit->nu.head(1), cit->fz.head(1), cit->x.head(1), cit->y.head(1)));
        std::auto_ptr<const ProjectedCone> cone (U_stance(cit->contacts, cit->positions, cit->nu, cit->fz, cit->x, cit->y));
    }
    time = get_time() - time;
    std::cout << "average time to generate two contact problem (one projection): \n"
              << time / (double)nbProblems << std::endl;
}

/*void timeGenerationSparse()
{
    int nbProblems = 1000;
    std::vector<Problem> problems;
    problems.resize(nbProblems);
    for(int i=0; i < nbProblems; ++i)
    {
        Problem& p = problems[i];
        for(int cid =0; cid <2; ++cid)
        {
            p.x(cid) = unifRand(2, 2);
            p.y(cid) = unifRand(2, 2);
            p.fz(cid) = 1;
            p.nu(cid) = unifRand(0.1,1.);
            Eigen::Matrix <value_type, 3, 3> rot = Eigen::Matrix <value_type, 3, 3>(
                      angle_axis_t(unifRand(M_PI), vector3_t::UnitX())
                    * angle_axis_t(unifRand(M_PI), vector3_t::UnitY())
                    * angle_axis_t(unifRand(M_PI), vector3_t::UnitZ()));
            vector3_t translation;
            for(int tid = 0; tid < 3; ++tid)
            {
                translation[tid] = unifRand(3.);
            }
            p.contacts.block<3,3>(cid*3,0) = rot;
            p.positions.segment<3>(cid*3) = translation;
        }
    }

    double time = get_time();
    int i = 0;
    for(std::vector<Problem>::const_iterator cit = problems.begin();
        cit != problems.end(); ++cit, ++i)
    {
        std::auto_ptr<const ProjectedCone> cone (U_stanceSparse(cit->contacts, cit->positions, cit->nu, cit->fz, cit->x, cit->y));
    }
    time = get_time() - time;
    std::cout << "average time to generate two contact problem (one projection, sparse): \n"
              << time / (double)nbProblems << std::endl;
}*/


int creationTest()
{
    int ret = 0;
    vector_t x(2); x << 1, 1;
    vector_t y(2); y << 1, 1;
    vector_t fz(2); fz << 1., 1.;
    vector_t nu(2); nu << 1, 0.5;
    T_rotation_t contact = Eigen::Matrix <value_type, 3,3>::Identity();
    T_rotation_t contacts = Eigen::Matrix <value_type, 6, 3>::Zero();
    contacts.block<3,3>(0,0) = contact;
    contacts.block<3,3>(3,0) = contact;
    vector_t position(6); position << 0,0,0,3,0,0;
//std::cout << "computed contact cone \n" << cone << std::endl;
    // assert all vectors are different
    matrix_t cone = V_all(x.head(1), y.head(1), fz.head(1), nu.head(1)).transpose();
    std::cout << "cone \n" << cone.transpose() << std::endl;
    for(int i = 0; i < 37; ++i)
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
    }

    const ProjectedCone* c2 = U_stance(contacts.block<3,3>(0,0), position.head<3>(),
                                nu.head(1),fz.head(1),x.head(1),y.head(1));
    vector3_t gravity(0,0,-9.81);
    value_type mass(40);
    vector3_t  goodP(0,0,3.);
    vector3_t  goodP2(1,0,0);
    vector3_t  badP(2,0,3.);

    std::cout << "good position verified ? " << c2->IsValid(goodP2,gravity,mass) << std::endl;
    std::cout << "good position verified ? " << c2->IsValid(goodP,gravity,mass) << std::endl;
    std::cout << "bad position verified ? " << c2->IsValid(badP,gravity,mass) << std::endl;
    const ProjectedCone* c = U_stance(contacts,position, nu,fz,x,y);
    std::cout << "good position verified ? " << c->IsValid(badP,gravity,mass) << std::endl;
    delete c;
    delete c2;
    return ret;
}

dd_MatrixPtr FromEigen (const matrix_t& input, dd_ErrorType *Error)
{
    // b - Ax > 0
    dd_rowrange m_input = (dd_rowrange)(input.rows());
    dd_colrange d_input = (dd_colrange)(input.cols());
    dd_MatrixPtr M=NULL;
    dd_rowrange i;
    dd_colrange j;
    dd_RepresentationType rep=dd_Inequality;
    mytype value;
    dd_NumberType NT = dd_Real;
    dd_init(value);
    M=dd_CreateMatrix(m_input, d_input);
    M->representation=rep;
    M->numbtype=NT;
    for (i = 0; i < m_input; i++)
    {
        for (j = 0; j < d_input; j++)
        {
            dd_set_d(value, input(i,j));
        dd_set(M->matrix[i][j],value);
        }
    }
  dd_clear(value);
  return M;
}

/*typedef Eigen::Matrix<value_type, 6, 1> vector6_t;
matrix_t raysFromInequalities(value_type nu, value_type x, value_type y, bool bound=false)
{
    // vector is fx, fy, fz, tx, ty, tz
    //friction cone on x
    matrix_t b_A = matrix_t::Zero(11 + (bound? 1 : 0 ),1+6);
    // - fx + ufz >= 0
    b_A(0,1) = -1; b_A(0,3) = nu;
    // fx + ufz >= 0
    b_A(1,1) = 1; b_A(1,3) = nu;
    // fy + ufz >= 0
    b_A(2,2) = 1; b_A(2,3) = nu;
    // -fy + ufz >= 0
    b_A(3,2) = -1; b_A(3,3) = nu;
    // fz >= 0
    b_A(4,3) = 1;
    // tx + yfz >= 0
    b_A(5,4) = 1; b_A(5,3) = y;
    // -tx + yfz >= 0
    b_A(6,4) = -1; b_A(6,3) = y;
    // -ty + xfz >= 0
    b_A(7,5) = -1; b_A(7,3) = x;
    // ty + xfz >= 0
    b_A(8,5) = 1; b_A(8,3) = x;
    // tz + ufz >= 0
    b_A(9,6) = 1; b_A(9,3) = nu;
    // -tz + ufz >= 0
    b_A(10,6) = -1; b_A(10,3) = nu;
    if(bound)
    {
        b_A(11,0) = 400; b_A(11,3) = -1;
    }

    dd_ErrorType error = dd_NoError;
    dd_MatrixPtr M = FromEigen (b_A, &error);
    if(error == dd_NoError)
    {
        dd_PolyhedraPtr poly = dd_DDMatrix2Poly(M, &error);

        if(error != dd_NoError)
        {
            std::cout << "marche pas " << std::endl;
        }
        dd_MatrixPtr G = dd_CopyGenerators(poly);
        // now copy that to eigen matrix
        int realRowSize = 0;
        matrix_t V = matrix_t::Zero(G->rowsize, G->colsize-1);
        for(int i=1; i <= G->rowsize; i++)
        {
            if(*(G->matrix[i-1][0]) == 0)
            {
                for(int j=2; j <= G->colsize; j++)
                {
                    V(realRowSize, j-2) = (value_type)(*(G->matrix[i-1][j-1]));
                }
                ++realRowSize;
            }
        }
        V = V.block(0,0,realRowSize, V.cols());

        realRowSize = 0;
        Eigen::MatrixXd Vertex = Eigen::MatrixXd::Zero(G->rowsize, G->colsize-1);
        for(int i=1; i <= G->rowsize; i++)
        {
            if(*(G->matrix[i-1][0]) == 1)
            {
                for(int j=2; j <= G->colsize; j++)
                {
                    Vertex(realRowSize, j-2) = (value_type)(*(G->matrix[i-1][j-1]));
                }
                ++realRowSize;
            }
        }
        Vertex = Vertex.block(0,0,realRowSize, Vertex.cols());
        std::cout << "generatrices \n" <<  V << std::endl;
        std::cout << "vertex \n" <<  Vertex << std::endl;

        {
            dd_MatrixPtr Ab = dd_CopyInequalities(poly);
            matrix_t A = matrix_t ((int)Ab->rowsize, (int)Ab->colsize-1);
            vector_t b = vector_t ((int)Ab->rowsize);
            for(int i=0; i < Ab->rowsize; ++i)
            {
                b(i) = (value_type)(*(Ab->matrix[i][0]));
                for(int j=1; j < Ab->colsize; ++j)
                {
                    A(i, j-1) = -(value_type)(*(Ab->matrix[i][j]));
                }
            }
            std::cout << "A \n" <<  A << std::endl;
            std::cout << "b \n" <<  b << std::endl;
            dd_FreeMatrix(G);dd_FreeMatrix(M);dd_FreePolyhedra(poly);
            dd_FreeMatrix(Ab);
            return V;
        }
    }
}

dd_MatrixPtr FromEigenV (const cref_matrix_t& input)
{
    dd_MatrixPtr M=NULL;
    dd_rowrange i;
    dd_colrange j;
    dd_rowrange m_input = (dd_rowrange)(input.rows());
    dd_colrange d_input = (dd_colrange)(7);
    dd_RepresentationType rep=dd_Generator;
    mytype value;
    dd_NumberType NT = dd_Real;;
    dd_init(value);

    M=dd_CreateMatrix(m_input, d_input);
    M->representation=rep;
    M->numbtype=NT;

    for (i = 1; i <= input.rows(); i++)
    {
        dd_set_d(value, 0);
        dd_set(M->matrix[i-1][0],value);
        for (j = 2; j <= d_input; j++)
        {
          dd_set_d(value, input(i-1,j-2));
          dd_set(M->matrix[i-1][j-1],value);
        }
    }
  dd_clear(value);
  return M;
}*/


/*matrix_t V_alltmp(value_type friction, value_type f_z_max,
                     value_type x, value_type y )
{
    const int c_dim = 6;
    const int num_gen = 32;
    const int nbContacts = 1;
    matrix_t cones = matrix_t::Zero(c_dim*nbContacts, num_gen*nbContacts);
    int colid = 0;
    for(int contact = 0; contact < nbContacts; ++contact, colid += num_gen)
    {
        const value_type f_z_max_c = 1;
        const value_type nu_f_z = f_z_max_c * friction;
        const value_type x_f_z =  x * f_z_max_c;
        const value_type y_f_z =  y * f_z_max_c;
        Eigen::Matrix <value_type, 6, 1> plus;
        plus[0] = nu_f_z;
        plus[1] = nu_f_z;
        plus[2] = f_z_max_c;
        plus[3] = y_f_z;
        plus[4] = x_f_z;
        plus[5] = nu_f_z;

        Eigen::Matrix <value_type, 6, 1> minus = - plus; minus[2] = f_z_max_c;
        std::size_t numblocks = 32;
        std::size_t length = 1;
        // generating linear combination of all non 0 points
        // first row alternates positive and negative values
        // second row switches signs every two columns
        // last row switches signs after 16 colums
        // this allows to generate all combinations
        for(int rowid=c_dim*contact; rowid<c_dim*contact+c_dim; ++rowid)
        {
            for(std::size_t j = 0; j< numblocks; j+=2)
            {
                for(std::size_t col = length*j; col < length*(j+1); ++col)
                {
                    cones(rowid,colid+col) = plus(rowid%c_dim);
                }
                for(std::size_t col = length*(j+1); col < length*(j+2); ++col)
                {
                    cones(rowid,colid+col) = minus(rowid%c_dim);
                }
            }
            if(rowid % c_dim != 1)
            {
                numblocks /= 2;
                length *= 2;
            }
        }
        // adding non null forces with 0 moment
        //cones.block<3,4>(c_dim*contact,33 + colid) = cones.block<3,4>(c_dim*contact, colid);
    }
    dd_MatrixPtr V_ = FromEigenV(cones);
    dd_ErrorType error;
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(V_, &error);
    dd_MatrixPtr Ab = dd_CopyInequalities(poly);
    matrix_t A = matrix_t ((int)Ab->rowsize, (int)Ab->colsize-1);
    vector_t b = vector_t ((int)Ab->rowsize);
    for(int i=0; i < Ab->rowsize; ++i)
    {
        b(i) = (double)(*(Ab->matrix[i][0]));
        for(int j=1; j < Ab->colsize; ++j)
        {
            A(i, j-1) = -(double)(*(Ab->matrix[i][j]));
        }
    }
    std::cout << "A \n" <<  A << std::endl;
    std::cout << "b \n" <<  b << std::endl;
    return cones;
}*/

int main()
{
    srand (time(NULL));
    int ret = 0;
    init_library();

    /*matrix_t v1 = V_alltmp(1,1,1,1);
    std::cout << "va_all" << v1 << std::endl;
    matrix_t v2 = raysFromInequalities(1,1,1).transpose();

    for(int i=0; i< v1.cols(); ++i)
    {
        vector6_t v16 = v1.block<6,1>(0,i);
        bool found = false;
        for (int j=0; j< v2.cols(); ++j )
        {
            vector6_t v26 = v2.block<6,1>(0,j);
            if((v16 - v26).norm() < std::numeric_limits<double>::epsilon() * 2)
            {
                found = true;
                break;
            }
        }
        if(!found) { std::cout << "not found v1 col " << i <<std::endl; }
    }

    for(int i=0; i< v2.cols(); ++i)
    {
        vector6_t v16 = v2.block<6,1>(0,i);
        bool found = false;
        for (int j=0; j< v1.cols(); ++j )
        {
            vector6_t v26 = v1.block<6,1>(0,j);
            if((v16 - v26).norm() < 1)
            {
                found = true;
                break;
            }
        }
        if(!found) { std::cout << "not found v2 col " << i <<std::endl; }
    }

   ret+= creationTest();*/
   timeGeneration();
   //timeGenerationSparse();
   /*  //timeGenerationProjection();

    value_type nu = 0.5;
    value_type x = 2;value_type y = 3;
    raysFromInequalities(1,1,1);



    /*raysFromInequalities(1,1,1,true);
    raysFromInequalities(0.5,1,1);
    raysFromInequalities(2,x,y);
    raysFromInequalities(nu,4,y);
    raysFromInequalities(nu,x,6);*/
    release_library();
    return ret;
}
