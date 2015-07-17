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

namespace polytope
{
//#define USE_FLOAT 1;
#ifdef USE_FLOAT
typedef float value_type;
#else
typedef double value_type;
#endif

typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
typedef Eigen::Matrix <value_type, 3, 3> rotation_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1> vector_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 3> T_rotation_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
//typedef Eigen::SparseMatrix<value_type> smatrix_t;

//define Eigen ref if available
#if EIGEN_VERSION_AT_LEAST(3,2,0)
typedef const Eigen::Ref<const vector_t>    & cref_vector_t;
typedef const Eigen::Ref<const rotation_t>  & cref_rotation_t;
typedef const Eigen::Ref<const matrix_t>    & cref_matrix_t;
typedef const Eigen::Ref<const vector3_t>   & cref_vector3_t;
typedef const Eigen::Ref<const T_rotation_t>& cref_T_rotation_t;
#else
typedef const vector_t    & cref_vector_t;
typedef const rotation_t  & cref_rotation_tt;
typedef const matrix_t    & cref_matrix_t;
typedef const vector3_t   & cref_vector3_t;
typedef const T_rotation_t& cref_T_rotation_t;
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
    // wether the current wrench is achievable (static equilibirum test)
    // TODO: lp instead of inequalities to validate a point ? this avoids the conversion
    bool IsValid(cref_vector3_t p_com, const cref_vector3_t gravity, const value_type& mass) const;

private:
    std::auto_ptr<const PImpl> pImpl_;

public:
    const matrix_t& V;
    const matrix_t& A;
    const vector_t& b;

    friend const ProjectedCone* fromGenerators(cref_matrix_t vRep);

private:
    ProjectedCone(const PImpl*);
    ProjectedCone(const ProjectedCone&); // no copy constructor
};

const ProjectedCone* fromGenerators(cref_matrix_t vRep);

/*16 generator rays per contact. If we set x = min(x,y),
we then only have 8 generators*/
matrix_t V_all (cref_vector_t frictions,
                cref_vector_t xs,
                cref_vector_t ys);

/*compute projection matrices into GICW. Positional factors of contact wrenches (rotation and position
of contacts are computed in world frame and stacked in a 6 X (6*nbContacts)  as follows:
[-R_0              0 ... -R_n              0]
[-(^(p_0))R_0   -R_0 ... -(^(p_n))R_n   -R_n]
*/
matrix_t A_stance(cref_T_rotation_t contacts, cref_vector_t positions);

/*compute polytope correspding to the stability margin given by the gravito inertial wrenchs*/
const ProjectedCone* U_stance(cref_T_rotation_t contacts, cref_vector_t positions,
                              cref_vector_t friction,cref_vector_t x, cref_vector_t y);

} //namespace planner
#endif //equilib
