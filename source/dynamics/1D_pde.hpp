#include "numeric/numeric.hpp"
#include "utility/array.hpp"
#include "algebra/tensor.hpp"


namespace Ariadne
{
    struct Parameter1D
    {
        FloatDP length;    //Length of string
        FloatDP tension;
        FloatDP mass;
        //FloatDP frequency; //Frequency of oscillation
        //FloatDP wavelength;    // n*L
        FloatDP damping;   //damping
        FloatDP CourantNumber; //Courant Number 
    };    

    // Solving the one dimensional pde
    //Tensor<2, FloatDP> pde_1Dsolver(std::function<FloatDP(FloatDP)> &phi0, std::function<FloatDP(FloatDP, FloatDP)>& source, Parameter1D& stringParameter, SizeType Nx);
    
    // Set initial condition
    //Tensor<2, FloatDP> setIC(Tensor<2, FloatDP>& uts, std::function<FloatDP(FloatDP)> &phi0, SizeType Nx, Array<FloatDP> spacePoint);

    //Array<FloatDP> linspace1D(FloatDP L, SizeType n);

class pde1D
{
private:
    
    Tensor<2, FloatDP> setIC(Tensor<2, FloatDP>& uts, std::function<FloatDP(FloatDP)> &phi0, SizeType Nx, Array<FloatDP> spacePoint);
    Array<FloatDP> linspace1D(FloatDP L, SizeType n);  
public:
    pde1D();
    ~pde1D();
    Tensor<2, FloatDP> pde_1Dsolver(std::function<FloatDP(FloatDP)> &phi0, std::function<FloatDP(FloatDP, FloatDP)>& source, Parameter1D& stringParameter, SizeType Nx);
};
pde1D::pde1D()
{
}
pde1D::~pde1D()
{
}



}   // namespace Ariadne

