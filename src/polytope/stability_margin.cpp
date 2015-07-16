#include "polytope/stability_margin.h"

#include "cdd/cddmp.h"
#include "cdd/setoper.h"
#include "cdd/cddtypes.h"
#include "cdd/cdd.h"

#include <Eigen/Sparse>

#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std;
namespace polytope
{

typedef Eigen::Matrix <value_type, 6, 1> vector6_t;

const vector3_t X(1,0,0);
const vector3_t Y(0,1,0);
const vector3_t Z(0,0,1);

// number of generators per contact
const int num_gen = 32;
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

matrix_t A_stance(cref_T_rotation_t contacts, cref_vector_t positions)
{
    int nbContacts = contacts.rows() / 3;
    assert(contacts.rows() % 3 == 0 && contacts.rows() == positions.rows());
    matrix_t mat = matrix_t::Zero(c_dim, nbContacts*c_dim);
    for(int i = 0; i< nbContacts; ++i)
    {
        const cref_rotation_t& mRi = contacts.block<3,3>(3*i,0);
        mat.block<3,3>(0,i*c_dim)   = -mRi;
        mat.block<3,3>(3,i*c_dim)   = -skew(positions.segment<3>(3*i)) *mRi;
        mat.block<3,3>(3,i*c_dim+3) = -mRi;
    }
    return mat;
}

matrix_t V_all(cref_vector_t friction, cref_vector_t x,
                     cref_vector_t y )
{
    const int nbContacts = x.rows();
    assert(y.rows() == nbContacts && friction.rows() == nbContacts);
    matrix_t cones = matrix_t::Zero(c_dim*nbContacts, num_gen*nbContacts);
    int colid = 0;
    for(int contact = 0; contact < nbContacts; ++contact, colid += num_gen)
    {
        const value_type nu_f_z = friction[contact];
        const value_type x_f_z =  x[contact];
        const value_type y_f_z =  y[contact];
        Eigen::Matrix <value_type, 6, 1> plus;
        plus[0] = nu_f_z;
        plus[1] = nu_f_z;
        plus[2] = 1;
        plus[3] = y_f_z;
        plus[4] = x_f_z;
        plus[5] = nu_f_z;

        Eigen::Matrix <value_type, 6, 1> minus = - plus; minus[2] = 1;
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
    }
    return cones;
}

void init_library()
{
    dd_set_global_constants();dd_debug = false;
}

void release_library()
{
    //dd_free_global_constants();
}

dd_MatrixPtr FromEigen (const cref_matrix_t& input)
{
    if (dd_debug)
    {
        std::cout << "from eigen input " << input << std::endl;
    }
    dd_debug = false;
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
}

struct PImpl
{
    PImpl(cref_matrix_t vRepresentation)
        : init_(false)
        , v(vRepresentation) {}

    ~PImpl()
    {
        dd_FreeMatrix(V_); dd_FreePolyhedra(H_); dd_FreeMatrix(b_A);
    }

    bool init()
    {
        V_ = FromEigen(v);
        dd_ErrorType error = dd_NoError;
        /*dd_rowset redset,impl_linset;
        dd_rowindex newpos;
        dd_MatrixCanonicalize(&V_, &impl_linset, &redset, &newpos, &error);
        if(error != dd_NoError)
        {
            std::cout << ("can not reduce matrix") << std::endl;
        }
        fprintf(stdout, "\nRedundant rows: ");
        set_fwrite(stdout, redset);
        fprintf(stdout, "\n");
        set_free(redset);
        set_free(impl_linset);
        free(newpos);
        error = dd_NoError;*/
        H_= dd_DDMatrix2Poly(V_, &error);
        if(error != dd_NoError)
        {
            if(dd_debug)
                std::cout << ("numerical instability in cddlib. ill formed polytope") << std::endl;
        }
        else
        {
            init_= true;
            b_A = dd_CopyInequalities(H_);
            // get equalities and add them as complementary inequality constraints
            long elem;
            std::vector<long> eq_rows;
            for(elem=1;elem<=b_A->linset[0];++elem)
            {
                if (set_member(elem,b_A->linset))
                   eq_rows.push_back(elem);
            }
            int rowsize = (int)b_A->rowsize;
            A = matrix_t (rowsize + eq_rows.size(), (int)b_A->colsize-1);
            b = vector_t (rowsize + eq_rows.size());
            for(int i=0; i < rowsize; ++i)
            {
                b(i) = (value_type)(*(b_A->matrix[i][0]));
                for(int j=1; j < b_A->colsize; ++j)
                {
                    A(i, j-1) = -(value_type)(*(b_A->matrix[i][j]));
                }
            }
            int i = 0;
            for(std::vector<long int>::const_iterator cit = eq_rows.begin();
                cit != eq_rows.end(); ++cit, ++i)
            {
                b(rowsize + i) = -b((int)(*cit));
                A(rowsize + i) = -A((int)(*cit));
            }
        }
        return init_;
    }

private:
    bool init_;

public:
    matrix_t v;
    // V representation of polytope
    dd_MatrixPtr V_;
    // H representation of polytope
    dd_MatrixPtr b_A;
    dd_PolyhedraPtr H_;
    matrix_t A;
    vector_t b;
};

ProjectedCone::ProjectedCone(const PImpl * pImpl)
    : pImpl_(pImpl)
    , V(pImpl_->v)
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

const ProjectedCone* fromGenerators(cref_matrix_t vRep)
{
    PImpl* pImpl = new PImpl(vRep);
    if(pImpl->init())
    {
        return new ProjectedCone(pImpl);
    }
    return 0;
}

const ProjectedCone* U_stance(cref_T_rotation_t contacts, cref_vector_t positions,
                              cref_vector_t friction, cref_vector_t x, cref_vector_t y)
{
    return fromGenerators((A_stance(contacts, positions)* V_all(friction, x, y)).transpose());
}
} // namespace polytope
