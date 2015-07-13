/**
* \file Test.h
* \brief Test for equilibirum functions
* \author Steve T.
* \version 0.1
* \date 29/07/2014
*
*/
#ifndef _STRUCT_TEST
#define _STRUCT_TEST

#include "polytope/config.h"

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/src/Core/util/Macros.h>
#include <memory>

namespace equilib
{
// if C++ 11, can use template typedef instead
#ifdef USE_FLOAT
typedef float value_type;
#else
typedef double value_type;
#endif

/*37 extreme points for cone (f_max constant) => 5 points which can take two non zero values (2⁵5)
+ 5 0 moment points ==> 2^5 +5)*/
typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1> vector_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 4> T_transform_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
typedef Eigen::SparseMatrix<value_type> smatrix_t;

//define Eigen ref if available
#if EIGEN_VERSION_AT_LEAST(3,2,0)
typedef const Eigen::Ref<const vector_t>& cref_vector_t;
typedef const Eigen::Ref<const vector3_t>& cref_vector3_t;
typedef const Eigen::Ref<const matrix_t>& cref_matrix_t;
typedef const Eigen::Ref<const smatrix_t>& cref_smatrix_t;
typedef const Eigen::Ref<const T_transform_t>& cref_T_transform_t;
#else
typedef const vector_t& cref_vector_t;
typedef const vector3_t& cref_vector3_t;
typedef const matrix_t& cref_matrix_t;
typedef const smatrix_t& cref_smatrix_t;
typedef const T_transform_t& cref_T_transform_t;
#endif


//must imperatively call this method to init polytope
// library !!!!!!!!!!!!!!!
void init_library();
void release_library();

// private implementation to
// encapsulate target polytope library
struct PImpl;
struct POLYTOPE_DLLAPI ProjectedCone
{
    // generate H and V representation in target
    // library format
    explicit ProjectedCone(cref_matrix_t vRepresentation);
    // wether the current wrench is achievable (static equilibirum test)
    bool IsValid(cref_vector3_t p_com, const cref_vector3_t gravity, const value_type& mass) const;

private:
    std::auto_ptr<const PImpl> pImpl_;

public:
    const matrix_t& A;
    const vector_t& b;
};

const ProjectedCone* U_stance(cref_T_transform_t contacts,
                                     cref_vector_t friction,cref_vector_t f_z_max,
                                     cref_vector_t x, cref_vector_t y);


// TODO test perfo for Sparse matrix for v_all and a_stance
/*compute V representation analiticaly for a a given set of contacts
assumes contact is rectangle, x and y describe half length of rectangle sides*/
matrix_t V_all (cref_vector_t friction,
                cref_vector_t f_z_max, cref_vector_t x,
                cref_vector_t y);

smatrix_t V_allsparse(cref_vector_t friction,
                cref_vector_t f_z_max, cref_vector_t x,
                cref_vector_t y);

/*compute projection matrices into GICW*/
matrix_t POLYTOPE_DLLAPI A_stance(cref_T_transform_t contacts);
} //namespace planner
#endif //equilib
