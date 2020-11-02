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

    Tensor<2, FloatDP> pde2D::setIC2D(std::function<FloatDP(FloatDP, FloatDP)>& phi0, SizeType Nx, SizeType Ny, Array<FloatDP> spacePointX, Array<FloatDP> spacePointY)
    {
        Tensor<2, FloatDP> u({Nx, Ny}, 0);
        for (SizeType i = 0; i < Nx; i++)
        {
            for (SizeType j = 0; j < Ny; j++)
            {
                u[{i, j}] = phi0(spacePointX[i], spacePointY[j]);
            }
        }
        return u;
    }


Tensor<3, FloatDP> pde2D::pde_2Dsolver(std::function<FloatDP(FloatDP, FloatDP)>& phi0, std::function<FloatDP(FloatDP, FloatDP, FloatDP)>& source, Parameter2D& firstDim, Parameter2D& secondDim, SizeType Nx, SizeType Ny)
    {
        auto spaceX = linspace2D(firstDim.length, Nx);    //Mesh point x
        auto spaceY = linspace2D(secondDim.length, Ny);   //Mesh point y

        auto dx = spaceX[1] - spaceX[0];
        auto dy = spaceY[1] - spaceY[0];

        FloatDP c = 10.0;
        FloatDP T = 5;

        FloatDP dt = 0.01;
        
        FloatDP Nt = round(T/dt);
        SizeType Ntime = Nt.get_d();
        auto time = linspace2D(T, Ntime);                 //Mesh point t

        FloatDP Cx = (c*dt/dx);                         //Courant Numbers
        FloatDP Cy = (c*dt/dy);
        FloatDP Cx2 = pow(Cx, 2);
        FloatDP Cy2 = pow(Cy, 2);

        Tensor<2, FloatDP> u({Nx, Ny}, 0);
        Tensor<3, FloatDP> uts({Nx, Ny, Ntime}, 0);

        std::cout << "Start Computing";
        u = setIC2D(phi0, Nx, Ny, spaceX, spaceY);
        for (SizeType x = 0; x < Nx; x++)
        {
            for (SizeType y = 0; y < Ny; y++)
            {
                uts[{x, y, 0}] = u[{x, y}];
            }   
        }

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