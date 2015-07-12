#include "polytope/DynamicStability.h"

#include "cddmp.h"
#include "setoper.h"
#include "cddtypes.h"
#include "cdd.h"

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
    matrix_t mat(nbContacts*6, nbContacts*6);
    for(int i = 0; i< nbContacts; ++i)
    {
        const rotation_t mRi = -contacts.block<3,3>(4*i,0);
        mat.block(i*6,i*6,3,3) = mRi;
        mat.block(i*6+3,i*6+3,3,3) = mRi;
        mat.block(i*6+3,i*6,3,3) = skew(contacts.block<3,1>(4*i,3)) *mRi;
    }
    //mat.makeCompressed();
    return mat;
}

matrix_t V_all(cref_vector_t friction, cref_vector_t f_z_max,
                     cref_vector_t x, cref_vector_t y )
{
    const int nbContacts = x.rows();
    assert(y.rows() == nbContacts && f_z_max.rows() == nbContacts &&
           friction.rows() == nbContacts);
    matrix_t cones = matrix_t::Zero(6*nbContacts, 37*nbContacts);
    int colid = 0;
    for(int contact = 0; contact < nbContacts; ++contact, colid += 37)
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
        for(int rowid=6*contact; rowid<6*contact+6; ++rowid)
        {
            for(std::size_t j = 0; j< numblocks; j+=2)
            {
                for(std::size_t col = length*j; col < length*(j+1); ++col)
                {
                    cones(rowid,colid+col) = plus(rowid%6);
                }
                for(std::size_t col = length*(j+1); col < length*(j+2); ++col)
                {
                    cones(rowid,colid+col) = minus(rowid%6);
                }
            }
            if(rowid % 6 != 1)
            {
                numblocks /= 2;
                length *= 2;
            }
        }
        // adding non null forces with 0 moment
        cones.block<3,4>(6*contact,33 + colid) = cones.block<3,4>(6*contact, colid);
    }
    return cones;
}

void init_library()
{
    dd_set_global_constants();
}

/*dd_MatrixPtr FromEigen (const cref_matrix_t& input)
{
    dd_MatrixPtr M=NULL;
    dd_rowrange i;
    dd_colrange j;
    dd_rowrange m_input = (dd_rowrange)(input.rows());
    dd_colrange d_input = (dd_colrange)(input.cols() + 1);
    dd_RepresentationType rep=dd_Generator;
    mytype value;
    dd_NumberType NT = dd_Real;;
    dd_init(value);

    M=dd_CreateMatrix(m_input, d_input);
    M->representation=rep;
    M->numbtype=NT;

    for (i = 1; i <= m_input; i++)
    {
        dd_set_d(value, 1);
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

// need to transpose matrix...
dd_MatrixPtr FromEigen (const cref_matrix_t& input)
{
    int nbContacts = input.cols() / 6;
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
        ++contact, rowOffset+=37, colOffset+=6)
    {
        for (i = 1; i <= 37; i++)
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
    PImpl(const dd_MatrixPtr V, const dd_PolyhedraPtr H)
        : V_(V), H_(H), hComputed_(false)
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
    ~PImpl()
    {
        dd_FreeMatrix(V_); dd_FreePolyhedra(H_); dd_FreeMatrix(b_A);
    }

    bool IsValid(const vector6_t& wrench)
    {
        matrix_t res = A * wrench - b;
        for(int i =0; i< res.rows(); ++i)
        {
            if(res(i) > 0)
            {
                std::cout << "invalid " << res << std::endl;
                 return false;
            }
        }
        return true;
    }

    // V representation of polytope
    dd_MatrixPtr V_;
    // H representation of polytope
    dd_MatrixPtr b_A;
    dd_PolyhedraPtr H_;
    matrix_t A;
    matrix_t b;
    bool hComputed_;
};

ProjectedCone::ProjectedCone(cref_matrix_t vRepresentation)
{
    dd_ErrorType error = dd_NoError;
    dd_MatrixPtr Vpc = FromEigen(vRepresentation);
    dd_PolyhedraPtr H = dd_DDMatrix2Poly(Vpc, &error);
    if(error != dd_NoError)
    {
        std::cout << "numerical inst " << std::endl;
    }
    pImpl_.reset(new PImpl(Vpc, H));
}

bool ProjectedCone::IsValid(cref_vector3_t p_com, const cref_vector3_t& gravity, const value_type& mass) const
{
    vector6_t wrench;
    wrench.block<3,1>(0,0)= p_com;
    wrench.block<3,1>(3,0)= (mass * p_com).cross(gravity);
    return pImpl_->IsValid(wrench);
}

matrix_t ProjectedCone::HRepresentation() const
{
    matrix_t res(pImpl_->A.rows(), pImpl_->A.cols()+1);
    res.block(0,0,pImpl_->A.rows(),pImpl_->A.cols()) = pImpl_->A;
    res.block(0,pImpl_->A.cols(),pImpl_->A.rows(),1) = pImpl_->b;
    return res;
}

ProjectedCone* U_stance(cref_T_transform_t contacts,
                                     cref_vector_t friction,cref_vector_t f_z_max,
                                     cref_vector_t x, cref_vector_t y)
{

    matrix_t projection = A_stance(contacts) *  V_all(friction, f_z_max, x, y);
    return new ProjectedCone(projection.transpose());
}


//dd_MatrixPtr Hc = dd_CopyInequalities(polyH);

/*
dd_MatrixPtr FromEigen (const matrix_t& b, const matrix_t& input, dd_ErrorType *Error)
{
    // b - Ax > 0
    dd_rowrange m_input = (dd_rowrange)(input.rows());
    dd_colrange d_input = (dd_colrange)(input.cols() + 1);
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

    for (i = 1; i <= m_input; i++)
    {
        dd_set_d(value, b(i-1,0));
        dd_set(M->matrix[i-1][0],value);
        for (j = 2; j <= d_input; j++)
        {
        #if defined GMPRATIONAL
            *Error=dd_NoRealNumberSupport;
            goto _L99;
        #else
            dd_set_d(value, -input(i-1,j-2));
        #endif
        dd_set(M->matrix[i-1][j - 1],value);
        }
    }
  dd_clear(value);
  return M;
}

dd_MatrixPtr FromEigen (const matrix_td& input, dd_ErrorType *Error)
{
    dd_MatrixPtr M=NULL;
    dd_rowrange i;
    dd_colrange j;
    dd_rowrange m_input = (dd_rowrange)(input.rows());
    dd_colrange d_input = (dd_colrange)(input.cols() + 1);
    dd_RepresentationType rep=dd_Generator;
    mytype value;
    dd_NumberType NT = dd_Real;;
    dd_init(value);

    M=dd_CreateMatrix(m_input, d_input);
    M->representation=rep;
    M->numbtype=NT;

    for (i = 1; i <= m_input; i++)
    {
        dd_set_d(value, 1);
        dd_set(M->matrix[i-1][0],value);
        for (j = 2; j <= d_input; j++)
        {
          dd_set_d(value, input(i-1,j-2));
          dd_set(M->matrix[i-1][j - 1],value);
        }
    }
  dd_clear(value);
  return M;
}*/

/*
// human motions analysis and simulation based on a
// general criterion of stability
void ComputePolytope(const T_transform_t& contactTransforms,
                     const vector_t &friction, const vector_t& f_z_max)
{
    std::size_t nbContacts = contactTransforms.rows() / 4;
    // for n contacts, we have 4n transforms
    assert(contactTransforms.rows() % 4 == 0);
    matrix_t C(3, (nbContacts) * 3);
    matrix_t A(3, (nbContacts) * 3);
    // Init A with identities
    Eigen::Matrix3d id = Eigen::Matrix3d::Identity();
    for(int k = 0; k < nbContacts + nbGrasps; ++k)
    {
        A.block<3,3>(0, 3*k) = id;
    }
    Eigen::MatrixXd Af, Ag;
    Eigen::VectorXd bc, bg;
    int afRows, afCols, agRows, agCols;
    agRows = nbGrasps * 6;
    agCols = nbGrasps * 3;
    afRows = (nbContacts > 0) ?  (nbContacts * 4 +1)  : 0;
    afCols = nbContacts * 3;
    if(nbContacts > 0)
    {
        Af = Eigen::MatrixXd::Zero(afRows, afCols);
        int acIndex = Af.rows()-1;
        bc = Eigen::VectorXd::Zero(afRows);
        bc(acIndex) = flimit;

        for(int i=0; i< nbContacts; ++i)
        {
            C.block<3,3>(0,3*i) = skew(contactTransforms[i].block<3,1>(0, 3));
            const Rotation& Ri = contactTransforms[i].block<3,3>(0, 0);
            const Vector vi = Ri*Y;
            const Vector si = Ri*X; // inversing Z and Y because Y up
            const Vector ti = Ri*Z;
            Eigen::Matrix4d beta_i;
            beta_i.block<1,3>(0,0) = -(friction * vi + si).transpose();
            beta_i.block<1,3>(1,0) = -(friction * vi - si).transpose();
            beta_i.block<1,3>(2,0) = -(friction * vi + ti).transpose();
            beta_i.block<1,3>(3,0) = -(friction * vi - ti).transpose();
            Af.block<4,3>(4*i,3*i) = beta_i.block<4,3>(0,0);
            Af(acIndex, 3*i + 1) = 1.; // +2 ??? Z Y
        }
    }
    if(nbGrasps > 0)
    {
        //Compute Ag
        Ag = Eigen::MatrixXd::Zero(agRows, agCols);
        bg = Eigen::VectorXd::Zero(agRows);
        for(int i=0; i< nbGrasps; ++i)
        {
            C.block<3,3>(0,3*(i+nbContacts)) = skew(graspTransforms[i].block<3,1>(0, 3));
            Ag(6*i  ,3*i)   = 1;
            Ag(6*i+1,3*i)   =-1;
            Ag(6*i+2,3*i+1) = 1;
            Ag(6*i+3,3*i+1) =-1;
            Ag(6*i+4,3*i+2) = 1;
            Ag(6*i+5,3*i+2) =-1;
            bg(6*i)   = maxGraspingForces(3*i);
            bg(6*i+1) = maxGraspingForces(3*i);
            bg(6*i+2) = maxGraspingForces(3*i+1);
            bg(6*i+3) = maxGraspingForces(3*i+1);
            bg(6*i+4) = maxGraspingForces(3*i+2);
            bg(6*i+5) = maxGraspingForces(3*i+2);
        }
    }
    //Eigen::MatrixXd Acg = Eigen::MatrixXd::Zero(Af.rows() + Ag.rows(), Af.cols() + Ag.cols());
    Eigen::MatrixXd Acg = Eigen::MatrixXd::Zero(afRows + agRows, afCols + agCols);
    Eigen::MatrixXd bcg(afRows + agRows, 1) ;
    if(nbContacts>0)
    {
        Acg.block(0,0,afRows,afCols) = Af;
        bcg.block(0,0,afRows,1) = bc;
    }
    if(nbGrasps>0)
    {
        Acg.block(afRows,afCols,agRows,agCols) = Ag;
        bcg.block(afRows,0,agRows,1) = bg;
    }
    //Project and compute H and h

    dd_set_global_constants();
    dd_ErrorType error;
    dd_MatrixPtr matrix = FromEigen(bcg, Acg, &error);
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &error);
    dd_MatrixPtr G = dd_CopyGenerators(poly);
    // now copy that to eigen matrix
    int realRowSize = 0;
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(G->rowsize, G->colsize-1);
    for(int i=1; i <= G->rowsize; i++)
    {
        if(*(G->matrix[i-1][0]) == 1)
        {
            for(int j=2; j <= G->colsize; j++)
            {
                V(realRowSize, j-2) = (double)(*(G->matrix[i-1][j-1]));
            }
            ++realRowSize;
        }
    }
    V = V.block(0,0,realRowSize, V.cols());
    Eigen::MatrixXd Pac(A.rows() + C.rows(), A.cols());
    Pac.block(0,0,A.rows(),A.cols()) = A;
    Pac.block(A.rows(),0,C.rows(),C.cols()) = C;
    Eigen::MatrixXd Vp = V * Pac.transpose();


    dd_MatrixPtr Vpc = FromEigen(Vp, &error);
    if(error != dd_NoError)
    {
        std::cout << "numerical inst " << std::endl;
        return false;
    }
    dd_PolyhedraPtr polyH = dd_DDMatrix2Poly(Vpc, &error);
    if(error != dd_NoError)
    {
        std::cout << "numerical inst " << std::endl;
        return false;
    }
    dd_MatrixPtr Hc = dd_CopyInequalities(polyH);
    Eigen::MatrixXd H(Hc->rowsize, Hc->colsize-1);
    Eigen::MatrixXd h(Hc->rowsize, 1);
    for(int i=1; i <= Hc->rowsize; i++)
    {
        h(i-1,0)= (double)(*(Hc->matrix[i-1][0]));
        for(int j=2; j <= Hc->colsize; j++)
        {
            H(i-1, j-2) = -(double)(*(Hc->matrix[i-1][j-1]));
        }
    }

    //Eigen::Vector3d gravity(0,0,-9.81);
    Eigen::Vector3d gravity(0,-9.81,0);
    Eigen::Vector3d W = mass * acceleration - mass * gravity;
    Eigen::MatrixXd H1 = H.block(0,0,H.rows(),3);
    Eigen::MatrixXd H2 = H.block(0,3,H.rows(),3);
    Eigen::MatrixXd res = (H1 + H2 * skew(comLocation)) * W - h;
    //now reconvert to H rep
    dd_FreeMatrix(G);dd_FreeMatrix(Hc);dd_FreeMatrix(matrix);dd_FreePolyhedra(poly);dd_FreePolyhedra(polyH);
    dd_free_global_constants();
    for(int i =0; i< res.rows();++i)
    {
        if(res(i)>1)
        {
            return false;
        }
    }
    return true;
}


// human motions analysis and simulation based on a
// general criterion of stability
// ajoute acceleration, mass et localisation
double ResidualRadius(const T_Transform &contactTransforms, const T_Transform &graspTransforms, const Eigen::VectorXd& maxGraspingForces,
                      const Eigen::Vector3d& acceleration, const Eigen::Vector3d& comLocation,
                      const float mass, float friction, double flimit)
{
    //Compute Ac
    // Ac = [Af, Av]'
        // Compute Af matches contact forces
        // Af = diag(B_1,B_N)
    int nbContacts = contactTransforms.size();
    int nbGrasps = graspTransforms.size();
    Eigen::MatrixXd C(3, (nbContacts + nbGrasps) * 3);
    Eigen::MatrixXd A(3, (nbContacts + nbGrasps) * 3);
    // Init A with identities
    Eigen::Matrix3d id = Eigen::Matrix3d::Identity();
    for(int k = 0; k < nbContacts + nbGrasps; ++k)
    {
        A.block<3,3>(0, 3*k) = id;
    }
    Eigen::MatrixXd Af, Ag;
    Eigen::VectorXd bc, bg;
    int afRows, afCols, agRows, agCols;
    agRows = nbGrasps * 6;
    agCols = nbGrasps * 3;
    afRows = (nbContacts > 0) ?  (nbContacts * 4 +1)  : 0;
    afCols = nbContacts * 3;
    if(nbContacts > 0)
    {
        Af = Eigen::MatrixXd::Zero(afRows, afCols);
        int acIndex = Af.rows()-1;
        bc = Eigen::VectorXd::Zero(afRows);
        bc(acIndex) = flimit;

        for(int i=0; i< nbContacts; ++i)
        {
            C.block<3,3>(0,3*i) = skew(contactTransforms[i].block<3,1>(0, 3));
            const Rotation& Ri = contactTransforms[i].block<3,3>(0, 0);
            const Vector vi = Ri*Z;
            const Vector si = Ri*X;
            const Vector ti = Ri*Y;
            Eigen::Matrix4d beta_i;
            beta_i.block<1,3>(0,0) = -(friction * vi + si).transpose();
            beta_i.block<1,3>(1,0) = -(friction * vi - si).transpose();
            beta_i.block<1,3>(2,0) = -(friction * vi + ti).transpose();
            beta_i.block<1,3>(3,0) = -(friction * vi - ti).transpose();
            Af.block<4,3>(4*i,3*i) = beta_i.block<4,3>(0,0);
            Af(acIndex, 3*i + 2) = 1.;
        }
    }
    if(nbGrasps > 0)
    {
        //Compute Ag
        Ag = Eigen::MatrixXd::Zero(agRows, agCols);
        bg = Eigen::VectorXd::Zero(agRows);
        for(int i=0; i< nbGrasps; ++i)
        {
            C.block<3,3>(0,3*(i+nbContacts)) = skew(graspTransforms[i].block<3,1>(0, 3));
            Ag(6*i  ,3*i)   = 1;
            Ag(6*i+1,3*i)   =-1;
            Ag(6*i+2,3*i+1) = 1;
            Ag(6*i+3,3*i+1) =-1;
            Ag(6*i+4,3*i+2) = 1;
            Ag(6*i+5,3*i+2) =-1;
            bg(6*i)   = maxGraspingForces(3*i);
            bg(6*i+1) = maxGraspingForces(3*i);
            bg(6*i+2) = maxGraspingForces(3*i+1);
            bg(6*i+3) = maxGraspingForces(3*i+1);
            bg(6*i+4) = maxGraspingForces(3*i+2);
            bg(6*i+5) = maxGraspingForces(3*i+2);
        }
    }
    //Eigen::MatrixXd Acg = Eigen::MatrixXd::Zero(Af.rows() + Ag.rows(), Af.cols() + Ag.cols());
    Eigen::MatrixXd Acg = Eigen::MatrixXd::Zero(afRows + agRows, afCols + agCols);
    Eigen::MatrixXd bcg(afRows + agRows, 1) ;
    if(nbContacts>0)
    {
        Acg.block(0,0,afRows,afCols) = Af;
        bcg.block(0,0,afRows,1) = bc;
    }
    if(nbGrasps>0)
    {
        Acg.block(afRows,afCols,agRows,agCols) = Ag;
        bcg.block(afRows,0,agRows,1) = bg;
    }
    //Project and compute H and h

    dd_set_global_constants();
    dd_ErrorType error;
    dd_MatrixPtr matrix = FromEigen(bcg, Acg, &error);
    dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &error);
    dd_MatrixPtr G = dd_CopyGenerators(poly);
    // now copy that to eigen matrix
    int realRowSize = 0;
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(G->rowsize, G->colsize-1);
    for(int i=1; i <= G->rowsize; i++)
    {
        if(*(G->matrix[i-1][0]) == 1)
        {
            for(int j=2; j <= G->colsize; j++)
            {
                V(realRowSize, j-2) = (double)(*(G->matrix[i-1][j-1]));
            }
            ++realRowSize;
        }
    }
    V = V.block(0,0,realRowSize, V.cols());
    Eigen::MatrixXd Pac(A.rows() + C.rows(), A.cols());
    Pac.block(0,0,A.rows(),A.cols()) = A;
    Pac.block(A.rows(),0,C.rows(),C.cols()) = C;
    Eigen::MatrixXd Vp = V * Pac.transpose();


    dd_MatrixPtr Vpc = FromEigen(Vp, &error);
    dd_PolyhedraPtr polyH = dd_DDMatrix2Poly(Vpc, &error);
	if(error != dd_NoError)
	{
		return -1000;
	}
    dd_MatrixPtr Hc = dd_CopyInequalities(polyH);
    Eigen::MatrixXd H(Hc->rowsize, Hc->colsize-1);
    Eigen::MatrixXd h(Hc->rowsize, 1);
    for(int i=1; i <= Hc->rowsize; i++)
    {
        h(i-1,0)= (double)(*(Hc->matrix[i-1][0]));
        for(int j=2; j <= Hc->colsize; j++)
        {
            H(i-1, j-2) = -(double)(*(Hc->matrix[i-1][j-1]));
        }
    }
    //Eigen::Vector3d gravity(0,0,-9.81);
    Eigen::Vector3d gravity(0,-9.81,0);
    Eigen::Vector3d W = mass * acceleration - mass * gravity;
    Eigen::MatrixXd H1 = H.block(0,0,H.rows(),3);
    Eigen::MatrixXd H2 = H.block(0,3,H.rows(),3);
    H = (H1 + H2 * skew(comLocation));
    Eigen::MatrixXd res = H * W - h;
    //now reconvert to H rep
    dd_FreeMatrix(G);dd_FreeMatrix(Hc);dd_FreeMatrix(matrix);dd_FreePolyhedra(poly);dd_FreePolyhedra(polyH);
    dd_free_global_constants();
    for(int i =0; i< res.rows();++i)
    {
        if(res(i)>0) return -1;
    }
    for(int i=0; i< H.rows(); ++i)
	{
		H.block(i,0,1,H.cols()) = H.block(i,0,1,H.cols()) / H.block(i,0,1,H.cols()).norm();
		h(i) = h(i) / H.block(i,0,1,H.cols()).norm();
	}
	return (h - H * W).minCoeff();
}*/
}
