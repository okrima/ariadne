
#include "numeric/numeric.hpp"
#include "numeric/real.hpp"
#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"
#include "function/function.decl.hpp"
#include "algebra/tensor.hpp"
//#include "typedefs.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

namespace Ariadne
{
    struct Parameter2D
    {
        FloatMP length;    //Length of string
        FloatMP damping;   //damping
    }; 

    // Solving the one dimensional pde
    Tensor<3, FloatMP> pde_2Dsolver(std::function<FloatMP(FloatMP, FloatMP)> &phi0, std::function<FloatMP(FloatMP, FloatMP, FloatMP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny);
    
    // Set initial condition
    Tensor<3, FloatMP> setIC(Tensor<3, FloatMP>& uts, std::function<FloatMP(FloatMP, FloatMP)> &phi0, SizeType Nx, SizeType Ny, Array<FloatMP> spacePointX, Array<FloatMP> spacePointY);

    // Create the linspace
    //Array<FloatMP> linspace(FloatMP L, SizeType n);

}