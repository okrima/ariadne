/***************************************************************************
 *            test_gnuplot.cpp
 * 
 ****************************************************************************
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../test.hpp"

#include "numeric/numeric.hpp"

#include "utility/gnuplot-iostream.h"
#include "output/gnuplot.hpp"
#include "dynamics/1D_pde.hpp"
#include "dynamics/2D_pde.hpp"

using namespace Ariadne;
using namespace std;

class TestGnuplot
{
    public:
        void test()
        {
            ARIADNE_TEST_CALL(LinFun());
            ARIADNE_TEST_CALL(defaultMultiExp());
            ARIADNE_TEST_CALL(defaultSincFunc());
            ARIADNE_TEST_CALL(stringAnimation());
            ARIADNE_TEST_CALL(gauss3D());
            ARIADNE_TEST_CALL(gauss3DAnimation());
        }

        void LinFun()
        {
            Gnuplot gp = Gnuplot("tee test_gnuplot-LinFunc.gnu | gnuplot -persist");   //commands dump and make graphics
            GnuplotCanvas canvas = GnuplotCanvas();

            Image2D img1, img2, img3, img4;
            SizeType nImg = 4;
            Array<Image2D> imgList(nImg);
            
            SizeType dim = 50;

            Array<FloatDP> data1(dim);
            Array<FloatDP> data2(dim);
            Array<FloatDP> data3(dim);
            Array<FloatDP> data4(dim);
            Array<Array<FloatDP>> dataList(nImg);

            std::function<FloatDP(FloatDP, FloatDP)> fun = [](FloatDP x, FloatDP c){return c*x;};

            //define size of the array
            for (SizeType i = 0; i < nImg; i++)
            {
                dataList[i].resize(dim);
            }
            
            for (SizeType i = 0; i < dim; i++)
            {
                dataList[0].at(i) = fun(i, 1.0);
                dataList[1].at(i) = fun(i, 0.5);
                dataList[2].at(i) = fun(i, 1.5);
                dataList[3].at(i) = fun(i, 2);
            }

            canvas.setTerminal(gp, _png, "test_gnuplot-LinFunc");

            canvas.setTitle(gp, "Test Linear Function Plot");
            canvas.setXLabel(gp, "x");
            canvas.setYLabel(gp, "y");

            canvas.setColour(img1, _black);
            canvas.setColour(img2, _red);
            canvas.setColour(img3, _dark_orange);
            canvas.setColour(img4, _web_green);

            _Line2D line1, line2, line3, line4;
            line1.style = lines2D;
            line1.ls = 1;
            line1.lw = 2;
            line2.style = linespoints2D;
            line3.style = points2D;
            line3.ls = 3;
            line3.lw = 1.25;
            line4.style = impulses2D;
            canvas.setLineStyle(img1, line1);
            canvas.setLineStyle(img2, line2, 2, 1);
            canvas.setLineStyle(img3, line3);
            canvas.setLineStyle(img4, line4, 3, 1.5);

            canvas.setMultiplot(gp, true);

            _Range2D range;
            canvas.setRange2D(range, 0, data1.size(), 0, 100);
            //or direct assignment to range.[..] = ...
            imgList[0] = img1;
            imgList[1] = img2;
            imgList[2] = img3;
            imgList[3] = img4;

            //Ready to plot
            canvas.plotArray2D(gp, imgList, range, dataList);
        }//Linear Function

        void defaultMultiExp()
        {
            Gnuplot gp = Gnuplot("tee test_gnuplot-MultiExp.gnu | gnuplot -persist");   //commands dump and make graphics
            GnuplotCanvas canvas = GnuplotCanvas();

            Image2D img1, img2;
            SizeType nImg = 2;
            Array<Image2D> imgList(nImg);
            
            SizeType dim = 50;

            Array<FloatDP> data1(dim);
            Array<FloatDP> data2(dim);
            Array<Array<FloatDP>> dataList(nImg);

            std::function<FloatDP(bool, FloatDP, FloatDP)> fun = [](bool x, FloatDP coeff, FloatDP c)
            {   if (x)
                {
                    return exp(coeff*c);
                }
                else
                {
                    return exp(-coeff*c);
                }   
            };

            for (SizeType i = 0; i < nImg; i++)
            {
                dataList[i].resize(dim);
            }
            
            for (SizeType i = 0; i < dim; i++)
            {
                dataList[0].at(i) = fun(true, 2.0, i);
                dataList[1].at(i) = fun(false, 1, i);
            }

            canvas.setTerminal(gp, _png, "test_gnuplot-MultiExp");

            //canvas.setTitle(gp, "Test Exponential Function Plot");
            canvas.setXLabel(gp, "x");
            canvas.setYLabel(gp, "y");

            //set default colour == black
            canvas.setColour(img1);
            canvas.setColour(img2);
            //set default linestyle -> style = lines, linewidth = 1, linestyle = 1
            canvas.setLineStyle(img1);
            canvas.setLineStyle(img2);

            canvas.setMultiplotLayout(gp, 2, 1, "Multiplot Layout with simple Exponential function");

            //set automatic range and put imgs inside array
            imgList[0] = img1;
            imgList[1] = img2;

            canvas.plotArray2D(gp, imgList, dataList);
        }//Exponential Function

        void defaultSincFunc()
        {
            Gnuplot gp = Gnuplot("tee test_gnuplot-Sinc.gnu | gnuplot -persist");   //commands dump and make graphics
            GnuplotCanvas canvas = GnuplotCanvas();

            Image2D img;
            SizeType dim = 100;

            Array<FloatDP> data(dim, 0);

            std::function<FloatDP(FloatDP)> fun = [](FloatDP x)
            {   if (x != 0)
                {
                    return sin(x)/x;
                }
                else
                {
                    return FloatDP(1);
                }   
            };

            for (SizeType i = 0; i < dim; i++)
            {
                data[i] = fun(i);
  
            }
        
            canvas.setTerminal(gp, _png, "test_gnuplot-Sinc");
            canvas.setTitle(gp, "Test Sinc Function Plotting");
            //No Labeling

            //Range from -50 to 50
            _Range2D range;
            canvas.setRange2D(range, 0, dim, -0.5, 1);

            canvas.plotArray2D(gp, img, range, data);
        }//Sinc Function

        void stringAnimation()
        {
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
        
            Gnuplot gp = Gnuplot("tee test_gnuplot-StringEvol.gnu | gnuplot -persist");
            GnuplotCanvas canvas = GnuplotCanvas();

            Image2D image;
            _Range2D range;

            canvas.setTerminal(gp, _gif, "test_gnuplot-StringEvol");
            canvas.setMultiplot(gp, false);
            canvas.setTitle(gp, "String Evolution with triangular IC");
            canvas.setXLabel(gp, "x - Space");
            canvas.setYLabel(gp, "Amplitude");
            canvas.setRange2D(range, 0, Nx+1, -1, 1);

            canvas.plotTensor2D(gp, image, range, data);
        }//String Evolution over time

        void gauss3D()
        {
            Gnuplot gp = Gnuplot("tee test_gnuplot-Gauss3D | gnuplot -persist");
            GnuplotCanvas canvas = GnuplotCanvas();

            Image3D gauss;
            SizeType dim = 50;
            Tensor<3, FloatDP> data({dim, dim, 1}, 0);

            std::function<FloatDP(FloatDP, FloatDP)> fun = [&](FloatDP x, FloatDP y){return (exp(-0.005*pow(x-dim/2, 2)-0.005*pow(y-dim/2, 2)));};

            for (SizeType step = 0; step < 1; step++)
            {
                for (SizeType x = 0; x < dim; x++)
                {
                    for (SizeType y = 0; y < dim; y++)
                    {
                        data[{x, y, step}] = fun(x, y);
                    }   
                }
            }

            _Range3D range;

            canvas.setTerminal(gp, _png, "test_gnuplot-Gauss3D");
            canvas.setTitle(gp, "3D plot");
            canvas.setXYZLabel(gp, "x", "y", "z");

            canvas.set3DPalette(gp, gauss, -0.5, 1, 0.2, true);

            canvas.setMultiplot(gp, false);
            
            canvas.setRange3D(range, 0, dim, 0, dim, 0, 1);

            canvas.plotTensor3D(gp, gauss, range, data);
        }//Gauss 3D

        void gauss3DAnimation()
        {
            SizeType Nx = 10;   //Mesh 1° dim
            SizeType Ny = 10;   //Mesh 2° dim

            pde2D solver2D = pde2D();

            Parameter2D firstDim, secondDim;

            firstDim.length = 10;
            secondDim.length = 10;

            
            std::function<FloatDP(FloatDP,FloatDP)> IC_Gauss = [&](FloatDP x, FloatDP y){return (exp(-0.5*pow(x-firstDim.length/2, 2)-0.5*pow(y-secondDim.length/2, 2)));};

            std::function<FloatDP(FloatDP, FloatDP, FloatDP)> source = [&](FloatDP x, FloatDP y, FloatDP t){return 0;};

            auto data = solver2D.pde_2Dsolver(IC_Gauss, source, firstDim, secondDim, Nx+1, Ny+1);

            Gnuplot gp = Gnuplot("tee test_gnuplot-GaussAnimation | gnuplot -persist");
            GnuplotCanvas canvas;

            Image3D image;
            _Range3D range3D;
            _Line3D line3D;
            line3D.style = surface3D;

            canvas.setTerminal(gp, _gif, "test_gnuplot-GaussAnimation");
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
        }
};

int main(int argc, const char** argv) {

    TestGnuplot testGnuplot = TestGnuplot();
    
    testGnuplot.test();

    return 0;
}