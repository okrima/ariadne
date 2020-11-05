#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "utility/array.hpp"

#include "output/gnuplot.hpp"
#include "dynamics/1D_pde.hpp"

#include "utility/gnuplot-iostream.h"

using namespace Ariadne;
using namespace std;

int main() {

    SizeType Nx = 10;             // # point in space

    pde1D solver1D = pde1D();

    //pde_solution_1D solution;
    Parameter1D stringModel;

    stringModel.length = 1.0;  // Length
    stringModel.tension = 100000;
    stringModel.mass = 4;
    FloatDP x0 = 0.85*stringModel.length;  // Point of max amplitube - Triangular IC
    FloatDP wavelength = 2*stringModel.length;
    //Real c = string.frequency*string.wavelength;
    FloatDP amp = 0.8;
    stringModel.damping = 100;

    FloatDP c = sqrt((stringModel.tension/(stringModel.mass/stringModel.length)));

    FloatDP frequency = c/wavelength; // Frequency

    //Real T = string.wavelength/c;
    stringModel.CourantNumber = 0.8;
    //Real C2 = pow(string.CourantNumber, 2);

    FloatDP k = 2*pi_opp()/wavelength;
    FloatDP omega = 2*pi_opp()*frequency;

    // Set function initial position condition
    //std::function<FloatDP(FloatDP)> PosInitSteady = [&](FloatDP x){return sin(k*x);};  //Steady-state wave
   
    std::function<FloatDP(FloatDP)> PosInitTri = [&](FloatDP x){ //Triangular
        if (x.get_d() < x0.get_d())
            return amp*(x/x0);
        else
            return amp*(stringModel.length - x)/(stringModel.length - x0);
    };

    std::function<FloatDP(FloatDP, FloatDP)> source = [&](FloatDP x, FloatDP t){return (sin(k*x)*cos(omega*t));};

    auto data = solver1D.pde_1Dsolver(PosInitTri, source, stringModel, Nx+1);
  

    //auto VelInitCond = [&](Real x){return 0;};            // Set function initial velocity position


    Gnuplot gp = Gnuplot("tee acoustic_1D.gnu | gnuplot -persist");

    GnuplotCanvas canvas = GnuplotCanvas();

    Image2D image;
    _Range2D range;

    canvas.setTerminal(gp, _gif, "string-TimeEvol");
    
    canvas.setMultiplot(gp, false);
    canvas.setTitle(gp, "String Evolution with triangular IC");
    canvas.setXLabel(gp, "x - Space");
    canvas.setYLabel(gp, "Amplitude");
    canvas.setRange2D(range, 0, Nx+1, -1, 1);

    //Array<FloatDP> space = linspace(stringModel.length, Nx+1);

    canvas.plotTensor2D(gp, image, range, data);

    return 0;
}