#include "numeric/numeric.hpp"
#include "numeric/real.hpp"
#include "algebra/algebra.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"

#include "output/gnuplot.hpp"
#include "dynamics/2D_pde.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

using namespace Ariadne;
using namespace std;

int main(){

    SizeType Nx = 10;   //Mesh 1° dim
    SizeType Ny = 10;   //Mesh 2° dim

    Parameter2D firstDim, secondDim;

    firstDim.length = 10;
    secondDim.length = 10;

    std::function<FloatMP(FloatMP,FloatMP)> IC = [&](FloatMP x, FloatMP y){return sin(pi_opp()*2*x/firstDim.length)*sin(pi_opp()*3*y/secondDim.length);};
    std::function<FloatMP(FloatMP,FloatMP)> IC_Gauss = [&](FloatMP x, FloatMP y){return (exp(-0.5*pow(x-firstDim.length/2, 2)-0.5*pow(y-secondDim.length/2, 2)));};

    std::function<FloatMP(FloatMP, FloatMP, FloatMP)> source = [&](FloatMP x, FloatMP y, FloatMP t){return 0;};

    auto data = pde_2Dsolver(IC, source, firstDim, secondDim, Nx+1, Ny+1);

    Gnuplot gp;
    GnuplotCanvas canvas;

    Image3D image;
    _Line3D line3D;
    line3D.style = surface3D;

    canvas.setTerminal(gp, image, _gif, "TimeEvol");
    canvas.setMultiplot(gp, false);
    canvas.setTitle(gp, "Evolution");
    canvas.setXLabel(gp, "x - Space");
    canvas.setYLabel(gp, "y - Space");
    canvas.setZLabel(gp, "Amplitude");
    canvas.setRange3D(image, Nx, Ny, 1);
    canvas.setLineStyle(image, line3D);
    canvas.set3DPalette(gp, image, true);
    //canvas.setMap(gp);

    canvas.plotTensor3D(gp, image, data, "temp.dat");

    return 0;
}

