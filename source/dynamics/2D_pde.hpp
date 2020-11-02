#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "algebra/tensor.hpp"

namespace Ariadne
{
    struct Parameter2D
    {
        FloatDP length;    //Length 
        FloatDP damping;   //damping
    }; 


class pde2D
{
private:
    // Set initial condition
    //Tensor<2, FloatDP> setIC2D(std::function<FloatDP(FloatDP, FloatDP)> phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY);
    Tensor<2, FloatDP> setIC2D(std::function<FloatDP(FloatDP, FloatDP)>& phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY);
    Array<FloatDP> linspace2D(FloatDP L, SizeType n);
    
public:
    pde2D();
    ~pde2D();

    // Solving the one dimensional pde
    Tensor<3, FloatDP> pde_2Dsolver(std::function<FloatDP(FloatDP, FloatDP)>& phi0, std::function<FloatDP(FloatDP, FloatDP, FloatDP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny);  

};
pde2D::pde2D()
{
}

pde2D::~pde2D()
{
}




}