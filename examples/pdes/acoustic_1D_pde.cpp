#include "numeric/numeric.hpp"
#include "numeric/real.hpp"
#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"

#include "output/gnuplot.hpp"
#include "dynamics/1D_pde.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

//#include "gnuplot-iostream.h"

using namespace Ariadne;
using namespace std;
//using  std::ofstream;

int main(int argc, const char** argv) {

    SizeType Nx = 41;                // # point in space

    //String tempStr;

    //tempStr = to_string(Nx.get_d() - 1);

    //Real Zero = 0;

    //pde_solution_1D solution;
    string1D stringModel;

    //Solution variable
    //Array<Real> u(Nx.get_d(), 0);   // Next value
    //Array<Real> u_1(Nx.get_d(), 0); // Current value
    //Array<Real> u_2(Nx.get_d(), 0); // Past value

    //solution.nextV = u;
    //solution.currV = u_1;
    //solution.prevV = u_2;

    //Array<double> sol1D(Nx.get_d());

    stringModel.length = 1.0;  // Length
    stringModel.tension = 100000;
    stringModel.massPerUnit = 1;
    FloatMP x0 = 0.85*stringModel.length;  // Point of max amplitube - Triangular IC
    stringModel.frequency = 440; // Frequency
    stringModel.wavelength = 2*stringModel.length;
    //Real c = string.frequency*string.wavelength;
    FloatMP amp = 0.8;
    stringModel.damping = 100;

    //Real T = string.wavelength/c;
    stringModel.CourantNumber = 0.8;
    //Real C2 = pow(string.CourantNumber, 2);

    FloatMP k = 2*pi_opp()/stringModel.wavelength;

/*
    Real L = 1.0_q;   // m - Lenght
    //Real x0 = 0.8_q*L;          // Triangular wave - point of max amplitude
    Real freq = 220;        // Oscillation frequency
    Real omega = 2*pi*freq;
    Real wavelength = 2*L;  //wavelength
    Real c = freq*wavelength;   // m/s - wave velocity
    Real amp = 0.5_q;           // amplitude
    Real damp = 15000;
*/
    //Real numPeriod = 1;
    //Real T = (2*pi/omega);
    //Real T = wavelength/c;            

    //Real C = 0.1_q;        // Courant Number <= 1
    //Real C2 = pow(C, 2);

    //Array<Real> space = linspace(string.length, Nx);    // Mesh in space

    //Real dx = space[1] - space[0];
    //Real dt = string.CourantNumber*dx/c;
    //Real Nt = round(T/dt);
    //SizeType _Nt = Nt.get_d();              
    //Array<Real> time = linspace(Nt*dt, Nt);

    // Set function initial position condition
    //std::function<FloatMP(FloatMP)> PosInitSteady = [&](FloatMP x){return sin(k*x);};  //Steady-state wave
    
    std::function<FloatMP(FloatMP)> PosInitTri = [&](FloatMP x){ //Triangular
        if (x.get_d() < x0.get_d())
            return amp*(x/x0);
        else
            return amp*(stringModel.length - x)/(stringModel.length - x0);
    };

    auto data = pde_1Dsolver(PosInitTri, stringModel, Nx);
  

    //auto VelInitCond = [&](Real x){return 0;};            // Set function initial velocity position
    //auto source = [&](Real x, Real t){return 0;};

    Gnuplot gp;
    GnuplotCanvas canvas = GnuplotCanvas();

    Image2D image;

    canvas.setTerminal(gp, image, _gif, "string-TimeEvol");
    canvas.setMultiplot(gp, false);
    canvas.setTitle(gp, "String Evolution with triangular IC");
    canvas.setXLabel(gp, "x - Space");
    canvas.setYLabel(gp, "Amplitude");
    canvas.setRange2D(image, 0, Nx, -1, 1);

    canvas.plotTensor2D(gp, image, data);

/*
    Gnuplot gp;
   
    gp << "set term gif animate\n";
    gp << "set output 'animate.gif'\n";
    gp << "set xlabel 'X'\n";
    gp << "set ylabel 'Amplitude'\n";

    ofstream outdata;
    outdata.open("data.dat");
    if (!outdata)
    {
        std::cerr << "Error: file couldn't be opened" << std::endl;
        exit(1);     
    }
    
    // Set the initial condition
    for (SizeType i = 0; i < Nx.get_d(); i++)
    {
        solution.currV[i] = PosInitTri(space[i]);   // Set initial condition - Steady-state 
        //u_1[i] = PosInitTri(space[i]);
        sol1D[i] = solution.currV[i].get_d();

        //outdata << u_1[i].get_d() << std::endl;
    }

    //sleep(1);
    //gp << "plot 'data.dat' \n";
    //outdata.close();

    //gp << "plot [0:" << tempStr << "][-1:1] '-' with linespoints ls 1\n";
    //gp.send1d(sol1D);
    //gp.flush();

    // Set Boundary Condition
    u[0] = Zero;
    u[Nx.get_d() - 1] = Zero;

    sol1D[0] = Zero.get_d();
    sol1D[Nx.get_d() - 1] = Zero.get_d();

    // Set first time step
    SizeType n = 0;
    for (SizeType i = 1; i < Nx.get_d() - 1; i++)
    {
        u[i] = u_1[i] //+ dt*VelInitCond(space[i]) 
            - 0.5_q*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]);
                //+ 0.5_q*pow(dt, 2)*source(space[i], time[n]);
        sol1D[i] = u[i].get_d();
    }
    u[0] = Zero;
    u[Nx.get_d() - 1] = Zero;

    sol1D[0] = Zero.get_d();
    sol1D[Nx.get_d() - 1] = Zero.get_d();

    gp << "plot [0:" << tempStr << "][-1:1] '-' with linespoints ls 1\n";
    gp.send1d(sol1D);

    gp.flush();
    gp << "plot '-' \n";
    gp.send1d(sol1D);

    // Before next step switch variables

    u_2 = u_1;
    u_1 = u;

    // For each time step
    for (n = 1; n < Nt.get_d(); n++)
    {
        // Update the inner points
        for (SizeType i = 1; i < Nx.get_d() - 1; i++)
        {
            u[i] = (1/(1+0.5_q*damp*dt))*(-(u_2[i]) + 1.1_q*u_1[i] + damp*dt/2*u_2[i] +
                C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]));
                    //+ pow(dt, 2)*source(space[i], time[n]));
            //std::cout << u[i].get_d() << "\n";

            sol1D[i] = u[i].get_d();
        }

        //cout << endl;
        // Insert Boundary Condition
        u[0] = Zero;
        u[Nx.get_d() - 1] = Zero;
        sol1D[0] = Zero.get_d();
        sol1D[Nx.get_d() - 1] = Zero.get_d();

        gp << "plot '-' [0:100]\n";
        gp.send1d(sol1D);
        gp.flush();
        
        gp << "plot [0:" << tempStr << "][-1:1] '-' with linespoints ls 1\n";
        gp.send1d(sol1D);
        sleep(0.2);

        // Before next step switch variables
        u_2 = u_1;
        u_1 = u;
       
        // Plotting w/ GNUPLOT
        gp << "plot '-' with lines 'ez' \n";
        gp.send1d(u);
        gp.flush();

    }

    //gp << "set output\n";
*/
    return 0;

}