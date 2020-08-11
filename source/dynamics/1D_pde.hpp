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
/*
    struct pde_solution_1D
    {
        Array<Real> nextV;  //Next solution
        Array<Real> currV;  //Current solution
        Array<Real> prevV;  //Previous solution
    };
*/
    struct string1D
    {
        FloatMP length;    //Length of string
        FloatMP tension;
        FloatMP massPerUnit;
        FloatMP frequency; //Frequency of oscillation
        FloatMP wavelength;    // n*L
        FloatMP damping;   //damping
        FloatMP CourantNumber; //Courant Number 
    };    


    // Solving the one dimensional pde
    Tensor<2, FloatMP> pde_1Dsolver(std::function<FloatMP(FloatMP)> &phi0, string1D& stringParameter, SizeType Nx);
    
    // Set initial condition
    Tensor<2, FloatMP> setIC(Tensor<2, FloatMP>& uts, std::function<FloatMP(FloatMP)> &phi0, SizeType Nx, Array<FloatMP> spacePoint);

    // Create the linespace
    Array<FloatMP> linspace(FloatMP L, SizeType n);

}   // namespace Ariadne