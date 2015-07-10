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

#include <vector>
#include <Eigen/Dense>

namespace equilib
{

typedef double value_type;
//typedef Eigen::Matrix <value_type, 3, 1> vector3_t;
//typedef Eigen::Matrix <value_type, 3, 3> rotation_t;
//typedef Eigen::Matrix <value_type, 4, 4> transform_t;
//typedef Eigen::Matrix <value_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
/*32 extreme points for cone (5 non constant points+ 0 ==> 2^5 +1)*/
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1> vector_t;
typedef Eigen::Matrix <value_type, Eigen::Dynamic, 4> T_transform_t;
typedef Eigen::Matrix <value_type, 6, 33> contact_cone_v_t;

// defining references
typedef const Eigen::Ref<const vector_t> EConstVector_t;
typedef const Eigen::Ref<const T_transform_t> EConstT_transform_t;

/*compute V representation analiticaly for a single contact*/
/*assumes contact is rectangle, x and y describe half length of rectangle sides*/
contact_cone_v_t ComputePolytope(const value_type& x, const value_type& y,
                                 const value_type& f_z_max, const value_type& friction);

/*void ComputePolytope(const EConstT_transform_t& contactTransforms,
                     const EConstVector_t& friction, const EConstVector_t& f_z_max);*/
} //namespace planner
#endif //equilib
