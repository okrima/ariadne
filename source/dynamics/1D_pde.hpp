//#include "ariadne.hpp"

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
    struct Parameter1D
    {
        FloatMP length;    //Length of string
        FloatMP tension;
        FloatMP mass;
        FloatMP frequency; //Frequency of oscillation
        FloatMP wavelength;    // n*L
        FloatMP damping;   //damping
        FloatMP CourantNumber; //Courant Number 
        FloatMP c; //velocity
    };    


    // Solving the one dimensional pde
    Tensor<2, FloatMP> pde_1Dsolver(std::function<FloatMP(FloatMP)> &phi0, std::function<FloatMP(FloatMP, FloatMP)>& source, Parameter1D& stringParameter, SizeType Nx);
    
    // Set initial condition
    Tensor<2, FloatMP> setIC(Tensor<2, FloatMP>& uts, std::function<FloatMP(FloatMP)> &phi0, SizeType Nx, Array<FloatMP> spacePoint);

    // Create the linspace
    Array<FloatMP> linspace(FloatMP L, SizeType n);

}   // namespace Ariadne