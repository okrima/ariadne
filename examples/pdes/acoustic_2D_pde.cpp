#include "numeric/numeric.hpp"

#include "algebra/vector.hpp"
#include "utility/array.hpp"

#include "output/gnuplot.hpp"
#include "dynamics/2D_pde.hpp"

#include "utility/gnuplot-iostream.h"

using namespace Ariadne;
using namespace std;

int main(){

    SizeType Nx = 10;   //Mesh 1° dim
    SizeType Ny = 10;   //Mesh 2° dim

    pde2D solver2D = pde2D();

    Parameter2D firstDim, secondDim;

    firstDim.length = 10;
    secondDim.length = 10;

    std::function<FloatDP(FloatDP,FloatDP)> IC = [&](FloatDP x, FloatDP y){return sin(pi_opp()*2*x/firstDim.length)*sin(pi_opp()*3*y/secondDim.length);};
    std::function<FloatDP(FloatDP,FloatDP)> IC_Gauss = [&](FloatDP x, FloatDP y){return (exp(-0.5*pow(x-firstDim.length/2, 2)-0.5*pow(y-secondDim.length/2, 2)));};

    std::function<FloatDP(FloatDP, FloatDP, FloatDP)> source = [&](FloatDP x, FloatDP y, FloatDP t){return 0;};

    auto data = solver2D.pde_2Dsolver(IC_Gauss, source, firstDim, secondDim, Nx+1, Ny+1);

    Gnuplot gp;
    GnuplotCanvas canvas;

    Image3D image;
    _Range3D range3D;
    _Line3D line3D;
    line3D.style = surface3D;

    canvas.setTerminal(gp, _gif, "TimeEvol3D");
    canvas.setMultiplot(gp, false);
    canvas.setTitle(gp, "Evolution");
    canvas.setXLabel(gp, "x - Space");
    canvas.setYLabel(gp, "y - Space");
    canvas.setZLabel(gp, "Amplitude");
    canvas.setRange3D(range3D, Nx, Ny, 1);
    canvas.setLineStyle(image, line3D);
    canvas.set3DPalette(gp, image, -1, 1, 0.2, true);
    //canvas.setMap(gp);

    canvas.plotTensor3D(gp, image, range3D, data);

    return 0;
}

