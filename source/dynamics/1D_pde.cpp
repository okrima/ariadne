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

    //Formattare la soluzione in modo da avere un array formato da tante tuple di vettori, in cui una
    //tupla sia la soluzione al primo istante di tempo, la seconda la soluzione al secondo istante 
    //di tempo e così via. Usare Tensor<2, Array<double> per definire tempo e spazio e usare un Array<double>
    //per definire il vettore nello spazio.

    //Solving one dimensional pde
    Tensor<2, FloatMP> pde_1Dsolver(std::function<FloatMP(FloatMP)>& phi0, string1D& stringParameter, SizeType Nx)
    {
        FloatMP c = sqrt((stringParameter.tension/stringParameter.massPerUnit));
        //Real c = stringParameter.frequency*stringParameter.wavelength;
        //Real c = 342;
        FloatMP T = 0.5;

        FloatMP C2 = pow(stringParameter.CourantNumber, 2);

        Array<FloatMP> space = linspace(stringParameter.length, Nx);

        FloatMP dx = space[1] - space[0];
        FloatMP dt = stringParameter.CourantNumber*dx/c;
        FloatMP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();

        Tensor<2, FloatMP> uts({Nx, Ntime}, 0);

        //Array<Real> time = linspace(Nt*dt, Nt);
        //DoublePrecision pr;
        FloatMP k = 2*pi_opp()/stringParameter.wavelength;

        uts = setIC(uts, phi0, Nx, space);

        // Set Boundary Condition
        //uts[{0, 0}] = 0;
        //uts[{Nx, 0}] = 0;
        //solution.nextV[0] = 0;
        //solution.nextV[Nx.get_d() - 1] = 0;

        // Set first Time step
        SizeType n = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            uts[{i, n+1}] = uts[{i, n}] - 0.5*C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}]);
/*
            solution.nextV[i] = solution.currV[i] //+ dt*velIC
                            - 0.5_q*C2*(solution.currV[i-1] - 2*solution.currV[i] + solution.currV[i+1]);
                            // + 0.5_q*pow(dt,2)*source(space[i], time[n])
*/
        }//first time step

/*
        solution.nextV[0] = 0;
        solution.nextV[Nx.get_d() - 1] = 0;
        
        //switch variables

        solution.prevV = solution.currV;
        solution.currV = solution.nextV;
*/
        // For each time step - main loop
        for (n = 1; n < Nt.get_d(); n++)
        {
            // Update inner points
            for (SizeType i = 1; i < Nx - 1; i++)
            {
                uts[{i, n+1}] = (1/(1+0.5*stringParameter.damping*dt))*(2*uts[{i, n}] - uts[{i, n-1}] + (stringParameter.damping*dt/2)*uts[{i, n-1}] +
                    C2*(uts[{i-1, n}] - 2*uts[{i, n}] + uts[{i+1, n}]));
/*
                solution.nextV[i] = (1/(1+0.5_q*stringParameter.damping*dt))*(-(solution.prevV[i])+1.1_q*solution.currV[i] + (stringParameter.damping*dt/2)*solution.prevV[i] +
                    C2*(solution.currV[i-1] - 2*solution.currV[i] + solution.currV[i+1]));
                    //+ pow(dt,2)*source(space[i], time[n]))
*/          
            }
/*
            solution.nextV[0] = 0;
            solution.nextV[Nx.get_d() - 1] = 0;

            //switch variables
            solution.prevV = solution.currV;
            solution.currV = solution.nextV;
*/
        }//main loop

        return uts;
    }
        
}// namespace Ariadne