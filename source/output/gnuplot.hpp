/***************************************************************************
 *            output/gnuplot.hpp
 *
 *  
 *
 ****************************************************************************/
/*
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



#include "../output/graphics.hpp"
#include "../config.hpp"
#include "../utility/string.hpp"
#include "algebra/tensor.hpp"
#include <string>
#include <ctype.h>

#include "gnuplot-iostream.h"

namespace Ariadne{

enum _Colours // %c - es: lc, fc
{
    _black,
    _dark_grey,
    _red,
    _light_red,
    _dark_red,
    _web_blue,
    _blue,
    _light_blue,
    _steelblue,
    _green,
    _dark_green,
    _light_green,
    _web_green,
    _dark_spring_green,
    _yellow,
    _dark_yellow,
    _magenta,
    _light_magenta,
    _dark_magenta,
    _cyan,
    _light_cyan,
    _dark_cyan,
    _dark_orange
};

const char *_colours[] = {"black", "dark-grey", "red", 
    "light-red", "dark-red", "web_blue", "blue", "light-blue",
    "steelblue", "green", "dark-green", "light-green", "web-green",
    "dark-spring-green", "yellow","dark-yellow", "magenta", 
    "light-magenta", "dark-magenta", "cyan", "light-cyan", 
    "dark-cyan", "dark-orange"};

enum _LineStyle2D // %s - es: ls, fs
{
    lines2D,
    linespoints2D, //lp
};

struct _Line2D
{
    _LineStyle2D style;
    int lw = 1;
    int ls = 1;
};

const char *_linestyle2D[] = {"lines", "linespoints"};// "points", "dots", "steps", "impulses"};

enum _LineStyle3D // %s - es: ls, fs
{
    lines3D,
    surface3D,
    pm3d,
    points3D,
};

struct _Line3D
{
    _LineStyle3D style;
    int lw = 1;
    int ls = 1;
};

const char *_linestyle3D[] = {"lines", "surface", "pm3d", "point",/* "impulse", "solid"*/};

enum _Format
{
    _png,
    _gif
};

const char *_format[] = {"png", "gif"};

struct _Range2D
{
    FloatMP Xmin = 0;
    FloatMP Xmax;
    FloatMP Ymin = 0;
    FloatMP Ymax;
};

struct _Range3D
{
    FloatMP Xmin = 0;
    FloatMP Xmax;
    FloatMP Ymin = 0;
    FloatMP Ymax;
    FloatMP Zmin = 0;
    FloatMP Zmax;
};

struct Image2D
{
    _Colours colour;
    _Line2D linestyle2D;
    _Range2D range2D;

};

struct Image3D
{
    _Colours colour;
    _Line3D linestyle3D;
    _Range3D range3D;
    
};

//Command
// show colorname
// test

class GnuplotCanvas
{
private:
    int sizeX;
    int sizeY;
    String _filename;
protected:  
    bool noCanvas;
    bool isMultiplot;
    bool is3DPalette;
public:
    ~GnuplotCanvas();
    // Constructors - Create the canvas
    //Create empty canvas
    GnuplotCanvas();
    //Create canvas with dimensions
    GnuplotCanvas(Image2D& image, int X, int Y);
    GnuplotCanvas(Image3D& image, int X, int Y);

    //Set Multiplot - Multiple plot on same screen
    void setMultiplot(Gnuplot& gp, bool s);
    //Plot 2D data from Tensor
    void plotTensor2D(Gnuplot& gp, Image2D& image, Tensor<2, FloatMP>& tensor);
    // Plot 3D data from Tensor
    void plotTensor3D(Gnuplot& gp, Image3D& image, Tensor<3, FloatMP>& tensor);
    // Plot from Array
    void plot2D(Gnuplot& gp, Image2D& image, Array<double> data);
    // 2D Plot from file
    void plot2D(Gnuplot& gp, Image2D& image, String filename);
    // 3D Plot with tensor - Evolution in Time
    void plot3D(Gnuplot& gp, Image3D& image, Array<Array<double>> data);
    // Set Terminal output 2D
    void setTerminal(Gnuplot& gp, Image2D& image, _Format format, String nameFile);
    // Set Terminal output 3D
    void setTerminal(Gnuplot& gp, Image3D& image, _Format format, String nameFile);
    // Set Title and Label
    void setXLabel(Gnuplot& gp, String xLabel);
    void setYLabel(Gnuplot& gp, String yLabel);
    void setZLabel(Gnuplot& gp, String zLabel);
    void setTitle(Gnuplot& gp, String title);
    void setXYZLabel(Gnuplot& gp, String xLabel, String yLabel, String zLabel);
    void setLabels(Gnuplot& gp, String xLabel, String yLabel, String zLabel, String title);
    // Set X, Y range
    void setRange2D(Image2D& image, FloatMP maxX, FloatMP maxY);
    // Set X, Y range
    void setRange2D(Image2D& image, FloatMP minX, FloatMP maxX, 
                FloatMP minY, FloatMP maxY);
    // Set X, Y, Z range
    void setRange3D(Image3D& image, FloatMP minX, FloatMP maxX, 
                FloatMP minY, FloatMP maxY,
                FloatMP minZ, FloatMP maxZ);
    void setRange3D(Image3D& image, FloatMP maxX, FloatMP maxY, FloatMP maxZ);
    // Set default 2D Linestyle
    void setLineStyle(Image2D& image);
    // Set default 3D linestyle
    void setLineStyle(Image3D& image);
    // Set Line style with parameter
    void setLineStyle(Image2D& image, _Line2D line, int ls, int lw);
    // Set Line style
    void setLineStyle(Image2D& image, _Line2D line);
    // Set line with parameter
    void setLineStyle(Image3D& image, _Line3D line, int ls, int lw);
    // set surface with no parameter
    void setLineStyle(Image3D& image, _Line3D line);
    // Set default 2D color
    void setColour(Image2D& image);
    // Set default 3D color
    void setColour(Image3D& image);
    // Set Colour with color
    void setColour(Image2D& image, _Colours color);
    // Set colour with color
    void setColour(Image3D& image, _Colours color);
    // Set X Log axis
    void setXLogAxis(Gnuplot& gp);
    // Set Y Log axis
    void setYLogAxis(Gnuplot& gp);
    // Set XY Log axis
    void setXYLogAxis(Gnuplot& gp);
    // Set XZ Log axis
    void setXZLogAxis(Gnuplot& gp);
    // Set YZ Log axis
    void setYZLogAxis(Gnuplot& gp);
    // Set XYZ Log axis
    void setXYZLogAxis(Gnuplot& gp);
    // Set Legend
    void setLegend(Gnuplot& gp);
    // Set View Projection of a 3D rapresentation
    void setMap(Gnuplot& gp);
    //Set 3D palette
    void set3DPalette(Gnuplot& gp, Image3D& image, FloatMP min, FloatMP max, bool s);
    //Unset colorbox
    void unsetColorbox(Gnuplot& gp);
    

};

GnuplotCanvas::~GnuplotCanvas()
{
}



} // namespace Ariadne