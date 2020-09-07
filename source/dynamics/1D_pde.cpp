#include "1D_pde.hpp"

//using namespace std;
namespace Ariadne{

    Array<FloatDP> pde1D::linspace1D(FloatDP L, SizeType n)
    {
        Array<FloatDP> linspaced(n);
        if (n == 0)
            return linspaced;
        if (n == 1)
        {
            linspaced[0] = L;
            return linspaced;
        }
        FloatDP delta = L/(n - 1);
        for (SizeType i = 0; i < (n - 1); i++)
        {
            linspaced[i] = (0 + delta*i);
        }
        linspaced[n - 1] = L;
        
        return linspaced;
    }

    // Set initial condition
    Tensor<2, FloatDP> pde1D::setIC(Tensor<2, FloatDP>& uts, std::function<FloatDP(FloatDP)> &phi0, SizeType Nx, Array<FloatDP> spacePoint)
    {
        for (SizeType i = 0; i < Nx; i++)
        {
            //solution.currV[i] = phi0(spacePoint[i]);
            uts[{i, 0}] = phi0(spacePoint[i]);
        }  
        return uts;
    }

    //Solving one dimensional pde
    Tensor<2, FloatDP> pde1D::pde_1Dsolver(std::function<FloatDP(FloatDP)>& phi0, std::function<FloatDP(FloatDP, FloatDP)>& source, Parameter1D& stringParameter, SizeType Nx)
    {
        FloatDP c = sqrt((stringParameter.tension/(stringParameter.mass/stringParameter.length)));
        //Real c = stringParameter.frequency*stringParameter.wavelength;
        //Real c = 342;
        FloatDP T = 0.125;

        FloatDP C2 = pow(stringParameter.CourantNumber, 2);

        Array<FloatDP> space = linspace1D(stringParameter.length, Nx);

        FloatDP dx = space[1] - space[0];
        FloatDP dt = stringParameter.CourantNumber*dx/c;
        FloatDP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();

        Array<FloatDP> time = linspace1D(T, Ntime);

        Tensor<2, FloatDP> uts({Nx, Ntime}, 0);

        //FloatDP k = 2*pi_opp()/stringParameter.wavelength;

        std::cout << "Computing";

        uts = setIC(uts, phi0, Nx, space);

        std::cout << " .";
        // Set first Time step
        SizeType n = 0;
        uts[{0, n}] = 0;
        uts[{Nx - 1, n}] = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            uts[{i, n+1}] = uts[{i, n}] - 0.5*C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}])
                + pow(dt,2)*source(space[i], time[n]);
        }// first time step

        std::cout << " .";
        // For each time step - main loop
        for (n = 1; n < Nt.get_d(); n++)
        {
            uts[{0, n}] = 0;
            uts[{Nx - 1, n}] = 0;
            // Update inner points
            for (SizeType i = 1; i < Nx - 1; i++)
            {
                uts[{i, n+1}] = (1/(1+0.5*stringParameter.damping*dt))*(2*uts[{i, n}] - uts[{i, n-1}] + (stringParameter.damping*dt/2)*uts[{i, n-1}] +
                    C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}])) + 1000*pow(dt, 2)*source(space[i], time[n]);        
            }
        }//main loop
        std::cout << " .";
        return uts;
    }
        
}// namespace Ariadne