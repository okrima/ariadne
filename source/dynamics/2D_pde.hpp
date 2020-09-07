
#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "algebra/tensor.hpp"

namespace Ariadne
{
    struct Parameter2D
    {
        FloatDP length;    //Length of string
        FloatDP damping;   //damping
    }; 

    // Solving the one dimensional pde
    //Tensor<3, FloatDP> pde_2Dsolver(std::function<FloatDP(FloatDP, FloatDP)> &phi0, std::function<FloatDP(FloatDP, FloatDP, FloatDP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny);
    
    // Set initial condition
    //Tensor<3, FloatDP> setIC(Tensor<3, FloatDP>& uts, std::function<FloatDP(FloatDP, FloatDP)> &phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY);

    //Array<FloatDP> linspace2D(FloatDP L, SizeType n);

class pde2D
{
private:
    // Set initial condition
    Tensor<3, FloatDP> setIC(Tensor<3, FloatDP>& uts, std::function<FloatDP(FloatDP, FloatDP)> &phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY);

    Array<FloatDP> linspace2D(FloatDP L, SizeType n);
public:
    pde2D();
    ~pde2D();

    // Solving the one dimensional pde
    Tensor<3, FloatDP> pde_2Dsolver(std::function<FloatDP(FloatDP, FloatDP)> &phi0, std::function<FloatDP(FloatDP, FloatDP, FloatDP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny);
    
};
pde2D::pde2D()
{
}

pde2D::~pde2D()
{
}



}