#include "1D_pde.hpp"

//using namespace std;

namespace Ariadne{

    Array<FloatMP> linspace(FloatMP L, SizeType n)
    {
        Array<FloatMP> linspaced(n);

        if (n == 0)
            return linspaced;
        if (n == 1)
        {
            linspaced[0] = L;
            return linspaced;
        }

        FloatMP delta = L/(n - 1);

        for (SizeType i = 0; i < (n - 1); i++)
        {
            linspaced[i] = (0 + delta*i);
        }
        linspaced[n - 1] = L;
        
        return linspaced;
    }
    // Set initial condition
    Tensor<2, FloatMP> setIC(Tensor<2, FloatMP>& uts, std::function<FloatMP(FloatMP)> &phi0, SizeType Nx, Array<FloatMP> spacePoint)
    {
        for (SizeType i = 0; i < Nx; i++)
        {
            //solution.currV[i] = phi0(spacePoint[i]);
            uts[{i, 0}] = phi0(spacePoint[i]);
        }  
        return uts;
    }

    //Solving one dimensional pde
    Tensor<2, FloatMP> pde_1Dsolver(std::function<FloatMP(FloatMP)>& phi0, std::function<FloatMP(FloatMP, FloatMP)>& source, Parameter1D& stringParameter, SizeType Nx)
    {
        FloatMP c = sqrt((stringParameter.tension/(stringParameter.mass/stringParameter.length)));
        //Real c = stringParameter.frequency*stringParameter.wavelength;
        //Real c = 342;
        FloatMP T = 0.125;

        FloatMP C2 = pow(stringParameter.CourantNumber, 2);

        Array<FloatMP> space = linspace(stringParameter.length, Nx);

        FloatMP dx = space[1] - space[0];
        FloatMP dt = stringParameter.CourantNumber*dx/c;
        FloatMP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();

        Array<FloatMP> time = linspace(T, Ntime);

        Tensor<2, FloatMP> uts({Nx, Ntime}, 0);

        //Array<Real> time = linspace(Nt*dt, Nt);
        //DoublePrecision pr;
        FloatMP k = 2*pi_opp()/stringParameter.wavelength;

        std::cout << "Computing";

        uts = setIC(uts, phi0, Nx, space);

        std::cout << " .";
        // Set first Time step
        SizeType n = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            uts[{i, n+1}] = uts[{i, n}] - 0.5*C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}])
                + pow(dt,2)*source(space[i], time[n]);
        }// first time step

        std::cout << " .";
        // For each time step - main loop
        for (n = 1; n < Nt.get_d(); n++)
        {
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