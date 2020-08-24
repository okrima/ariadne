#include "2D_pde.hpp"
#include "1D_pde.hpp"

namespace Ariadne
{
/*    
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
*/
    Tensor<3, FloatMP> setIC(Tensor<3, FloatMP>& uts, std::function<FloatMP(FloatMP, FloatMP)> &phi0, SizeType Nx, SizeType Ny, Array<FloatMP> spacePointX, Array<FloatMP> spacePointY)
    {
        for (SizeType i = 1; i < Nx -1; i++)
        {
            for (SizeType j = 1; j < Ny -1; j++)
            {
                uts[{i, j, 0}] = phi0(spacePointX[i], spacePointY[j]);
            }
        }
        return uts;
    }

    Tensor<3, FloatMP> pde_2Dsolver(std::function<FloatMP(FloatMP, FloatMP)> &phi0, std::function<FloatMP(FloatMP, FloatMP, FloatMP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny)
    {
        auto spaceX = linspace(firstDim.length, Nx);    //Mesh point x
        auto spaceY = linspace(secondDim.length, Ny);   //Mesh point y

        auto dx = spaceX[1] - spaceX[0];
        auto dy = spaceY[1] - spaceY[0];
        auto dx2 = pow(dx, 2);
        auto dy2 = pow(dy, 2);

        FloatMP c = 10.0;
        FloatMP T = 20;

        FloatMP stability_limit = (1/c)*(1/sqrt(1/dx2 + 1/dy2));
        FloatMP dt = 0.01;
        //ARIADNE_ASSERT(dt < stability_limit);
        if (dt <= 0)
        {
            FloatMP safetyFact = -dt;
            dt = safetyFact*stability_limit;
        }
        

        FloatMP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();

        auto time = linspace(T, Ntime);                 //Mesh point t

        FloatMP Cx = (c*dt/dx);                         //Courant Numbers
        FloatMP Cy = (c*dt/dy);
        FloatMP Cx2 = pow(Cx, 2);
        FloatMP Cy2 = pow(Cy, 2);

        Tensor<3, FloatMP> uts({Nx, Ny, Ntime}, 0);

        std::cout << "Start Computing";
        //Set boundary condition
    /*
        for (SizeType i = 0; i < Nx; i++)
        {
            uts[{i, 0}] = 0;
            uts[{i, Ny - 1}] = 0;
        }
        for (SizeType i = 0; i < Ny; i++)
        {
            uts[{0, i}] = 0;
            uts[{Nx - 1, i}] = 0;
        }
    */
        uts = setIC(uts, phi0, Nx, Ny, spaceX, spaceY);
        std::cout << " .";
        //Compute first time step
        SizeType n = 0;
        for (SizeType i = 1; i < Nx - 1; i++)
        {
            for (SizeType j = 1; j < Ny - 1; j++)
            {
                auto u_xx = uts[{i-1, j, n}] - 2*uts[{i, j, n}] + uts[{i+1, j, n}];
                auto u_yy = uts[{i, j-1, n}] - 2*uts[{i, j, n}] + uts[{i, j+1, n}];
                uts[{i, j, n+1}] = uts[{i, j, n}] + 0.5*Cx2*u_xx + 0.5*Cy2*u_yy + pow(dt, 2)*source(spaceX[i], spaceY[j], time[n]);
                //uts[{i, j, n+1}] += dt*VelIC(spaceX[i], spaceY[j]);
            }
        }
        std::cout << " .";
    /*
        //Set boundary condition
        for (SizeType i = 0; i < Nx; i++)
        {
            uts[{i, 0}] = 0;
            uts[{i, Ny - 1}] = 0;
        }
        for (SizeType i = 0; i < Ny; i++)
        {
            uts[{0, i}] = 0;
            uts[{Nx - 1, i}] = 0;
        }
    */
        //Compute each time step
        for (n = 1; n < Ntime; n++)
        {
            for (SizeType i = 1; i < Nx - 1; i++)
            {
                for (SizeType j = 1; j < Ny - 1; j++)
                {
                    auto u_xx = uts[{i-1, j, n}] - 2*uts[{i, j, n}] + uts[{i+1, j, n}];
                    auto u_yy = uts[{i, j-1, n}] - 2*uts[{i, j, n}] + uts[{i, j+1, n}];
                    uts[{i, j, n+1}] = 2*uts[{i, j, n}] - uts[{i, j, n-1}] + Cx2*u_xx + Cy2*u_yy + pow(dt, 2)*source(spaceX[i], spaceY[j], time[n]);
                }
                
            }
        }
        std::cout << " ." << std::endl;
    return uts;
    }
}//namespace Ariadne