#include "2D_pde.hpp"

namespace Ariadne
{
    Array<FloatDP> pde2D::linspace2D(FloatDP L, SizeType n)
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

    Tensor<3, FloatDP> pde2D::setIC(Tensor<3, FloatDP>& uts, std::function<FloatDP(FloatDP, FloatDP)> &phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY)
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

    Tensor<3, FloatDP> pde2D::pde_2Dsolver(std::function<FloatDP(FloatDP, FloatDP)> &phi0, std::function<FloatDP(FloatDP, FloatDP, FloatDP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny)
    {
        auto spaceX = linspace2D(firstDim.length, Nx);    //Mesh point x
        auto spaceY = linspace2D(secondDim.length, Ny);   //Mesh point y

        auto dx = spaceX[1] - spaceX[0];
        auto dy = spaceY[1] - spaceY[0];
        auto dx2 = pow(dx, 2);
        auto dy2 = pow(dy, 2);

        FloatDP c = 10.0;
        FloatDP T = 5;

        FloatDP stability_limit = (1/c)*(1/sqrt(1/dx2 + 1/dy2));
        FloatDP dt = 0.01;
        if (dt <= 0)
        {
            FloatDP safetyFact = -dt;
            dt = safetyFact*stability_limit;
        }
        
        FloatDP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();

        auto time = linspace2D(T, Ntime);                 //Mesh point t

        FloatDP Cx = (c*dt/dx);                         //Courant Numbers
        FloatDP Cy = (c*dt/dy);
        FloatDP Cx2 = pow(Cx, 2);
        FloatDP Cy2 = pow(Cy, 2);

        Tensor<3, FloatDP> uts({Nx, Ny, Ntime}, 0);

        std::cout << "Start Computing";

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