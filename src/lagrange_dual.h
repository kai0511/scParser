#ifndef LAGRANGE_DUAL_HPP
#define LAGRANGE_DUAL_HPP


#include "../inst/include/SR2_types.h"

mat lagrange_dual(const mat& gram, const mat& cdtx, const int& reg /*= 1*/, 
                  const int& max_iter /*= 100*/, const double& tol /*= 1e-5*/);

#endif