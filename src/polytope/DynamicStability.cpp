#include "polytope/DynamicStability.h"

#include "cdd/cddmp.h"
#include "cdd/setoper.h"
#include "cdd/cddtypes.h"
#include "cdd/cdd.h"

#include <Eigen/Sparse>

#include <vector>
#include <iostream>

using namespace std;
namespace equilib
{

typedef Eigen::Matrix <value_type, 6, 1> vector6_t;
typedef Eigen::Matrix <value_type, 3, 3> rotation_t;
typedef Eigen::Matrix <value_type, 4, 4> transform_t;

typedef Eigen::Triplet<value_type> T;

const vector3_t X(1,0,0);
const vector3_t Y(0,1,0);
const vector3_t Z(0,0,1);

// number of generators per contact
const int num_gen = 37;
// contact dimension
const int c_dim = 6;

rotation_t skew(cref_vector3_t& x)
{
    rotation_t res = rotation_t::Zero();
    res(0,1) = - x(2); res(0,2) =   x(1);
    res(1,0) =   x(2); res(1,2) = - x(0);
    res(2,0) = - x(1); res(2,1) =   x(0);
    return res;
}

// TODO TEST faster with sparse matrix despite init ?

matrix_t A_stance(cref_T_transform_t contacts)
{
    int nbContacts = contacts.rows() / 4;
    assert(contacts.rows() %4 == 0);
    //Eigen::SparseMatrix<value_type> mat(nbContacts*6,nbContacts*6);
    //reserve non zero values
    //mat.reserve(ReserveSparse());
    matrix_t mat(nbContacts*c_dim, nbContacts*c_dim);
    for(int i = 0; i< nbContacts; ++i)
    {
        const rotation_t mRi = -contacts.block<3,3>(4*i,0);
        mat.block(i*c_dim,i*c_dim,3,3) = mRi;
        mat.block(i*c_dim+3,i*c_dim+3,3,3) = mRi;
        mat.block(i*c_dim+3,i*c_dim,3,3) = skew(contacts.block<3,1>(4*i,3)) *mRi;
    }
    //mat.makeCompressed();
    //std::cout << "mat" << mat << std::endl;
    return mat;
}

matrix_t V_all(cref_vector_t friction, cref_vector_t f_z_max,
                     cref_vector_t x, cref_vector_t y )
{
    const int nbContacts = x.rows();
    assert(y.rows() == nbContacts && f_z_max.rows() == nbContacts &&
           friction.rows() == nbContacts);
    matrix_t cones = matrix_t::Zero(c_dim*nbContacts, num_gen*nbContacts);
    int colid = 0;
    for(int contact = 0; contact < nbContacts; ++contact, colid += num_gen)
    {
        const value_type& f_z_max_c = f_z_max[contact];
        const value_type nu_f_z = f_z_max_c * friction[contact];
        const value_type x_f_z =  x[contact] * f_z_max_c;
        const value_type y_f_z =  y[contact] * f_z_max_c;
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
        cones.block<3,4>(c_dim*contact,33 + colid) = cones.block<3,4>(c_dim*contact, colid);
    }
    //std::cout << " v_ all \n" << cones << "\n end v_all" << std::endl;
    return cones;
}

smatrix_t V_allsparse(cref_vector_t friction, cref_vector_t f_z_max,
                     cref_vector_t x, cref_vector_t y )
{
    const int nbContacts = x.rows();
    assert(y.rows() == nbContacts && f_z_max.rows() == nbContacts &&
           friction.rows() == nbContacts);
    smatrix_t cones = smatrix_t(c_dim*nbContacts, num_gen*nbContacts);
    cones.reserve(Eigen::VectorXi::Constant(num_gen*nbContacts,6));
    int colid = 0;
    for(int contact = 0; contact < nbContacts; ++contact, colid += num_gen)
    {
        const value_type& f_z_max_c = f_z_max[contact];
        const value_type nu_f_z = f_z_max_c * friction[contact];
        const value_type x_f_z =  x[contact] * f_z_max_c;
        const value_type y_f_z =  y[contact] * f_z_max_c;
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
                    cones.insert(rowid,colid+col) = plus(rowid%c_dim);
                }
                for(std::size_t col = length*(j+1); col < length*(j+2); ++col)
                {
                    cones.insert(rowid,colid+col) = minus(rowid%c_dim);
                }
            }
            if(rowid % c_dim != 1)
            {
                numblocks /= 2;
                length *= 2;
            }
        }
        // adding non null forces with 0 moment
        for(int r = 0; r<3; ++r)
        {
            for(int c = 0; c<4; ++c)
            {
                cones.insert(r+c_dim*contact,c+33 + colid) = cones.coeff(r+c_dim*contact,c+ colid);
            }
        }
        //cones.block<3,4>(c_dim*contact,33 + colid) = cones.block<3,4>(c_dim*contact, colid);
    }
    return cones;
}

void init_library()
{
    dd_set_global_constants();//dd_debug = true;
}

void release_library()
{
    // TODO
    //dd_free_global_constants();
}

dd_MatrixPtr FromEigen (const cref_matrix_t& input)
{
    //std::cout << "matrice \n " << input << "wtf"  << std::endl;
    int nbContacts = input.cols() / c_dim;
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

    int rowOffset = 0; int colOffset = 0;
    for(int contact = 0; contact < nbContacts;
        ++contact, rowOffset+=num_gen, colOffset+=c_dim)
    {
        for (i = 1; i <= num_gen; i++)
        {
            dd_set_d(value, 1);
            dd_set(M->matrix[rowOffset+i-1][0],value);
            for (j = 2; j <= d_input; j++)
            {
              dd_set_d(value, input(rowOffset+i-1,colOffset+j-2));
              dd_set(M->matrix[rowOffset+i-1][j-1],value);
            }
        }
    }
  dd_clear(value);
  return M;
}

struct PImpl
{
    PImpl(cref_matrix_t vRepresentation)
    {
        V_ = FromEigen(vRepresentation);
        dd_ErrorType error = dd_NoError;
        H_= dd_DDMatrix2Poly(V_, &error);
        if(error != dd_NoError)
        {
            std::cout << "numerical inst " << std::endl;
        }
        //else
        {
            b_A = dd_CopyInequalities(H_);
            A = matrix_t ((int)b_A->rowsize, (int)b_A->colsize-1);
            b = vector_t ((int)b_A->rowsize);
            for(int i=0; i < b_A->rowsize; ++i)
            {
                b(i) = (double)(*(b_A->matrix[i][0]));
                for(int j=1; j < b_A->colsize; ++j)
                {
                    A(i, j-1) = -(double)(*(b_A->matrix[i][j]));
                }
            }
        }
    }
    ~PImpl()
    {
        dd_FreeMatrix(V_); dd_FreePolyhedra(H_); dd_FreeMatrix(b_A);
    }

    // V representation of polytope
    dd_MatrixPtr V_;
    // H representation of polytope
    dd_MatrixPtr b_A;
    dd_PolyhedraPtr H_;
    matrix_t A;
    vector_t b;
    bool hComputed_;
};

ProjectedCone::ProjectedCone(cref_matrix_t vRepresentation)
    : pImpl_(new PImpl(vRepresentation))
    , A(pImpl_->A)
    , b(pImpl_->b)
{
    // NOTHING
}

bool ProjectedCone::IsValid(cref_vector3_t p_com, const cref_vector3_t gravity, const value_type& m) const
{
    vector6_t wrench;
    wrench.block<3,1>(0,0)= m*gravity;
    wrench.block<3,1>(3,0)= p_com.cross(m*gravity);
    matrix_t res = A * wrench - b;
    for(int i =0; i< res.rows(); ++i)
    {
        if(res(i) > 0)
        {
             return false;
        }
    }
    return true;
}

const ProjectedCone *U_stance(cref_T_transform_t contacts,
                                     cref_vector_t friction,cref_vector_t f_z_max,
                                     cref_vector_t x, cref_vector_t y)
{
    matrix_t projection = V_all(friction, f_z_max, x, y).transpose() * A_stance(contacts);
    return  new ProjectedCone(projection);
}
}
