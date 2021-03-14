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
#include "numeric/float_bounds.hpp"
#include "geometry2d.hpp"
#include <string>
#include <ctype.h>

#ifdef HAVE_GNUPLOT_H

#include "gnuplot-iostream.hpp"

namespace Ariadne{

struct  _Range
{       
    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax; 
};

struct  _Labels
{
    String xLabel = "";
    String yLabel = "";
    String zLabel = "";
};

class GnuplotCanvas : public CanvasInterface
{
friend class Figure;
private:
    Gnuplot *gnuplot;
    List<Point2d> geom;
    Colour lc, fc;
    double lw;
    _Range rng;
    Nat dim;
    Point2d Cpoint;
    double dr;
    bool isdot;
    Nat sizeX;
    Nat sizeY;
    bool isMultiplot;
    bool is2DPalette;
    bool is3DPalette;
    _Labels labels;
    GnuplotFileType fileType;

private:
    Void plot2D(Array<double> data);

    Void plot2D(Array<Array<double>> data);

    Void plot3D(Array<Array<double>> data);

public:
    ~GnuplotCanvas();
    // Constructors - Create the canvas
    //Create canvas with dimensions
    GnuplotCanvas(String filename, GnuplotFileType fileType, Nat X = 800, Nat Y = 800);

    //CanvasInterface
    Void initialise(StringType x, StringType y, StringType z, double xl, double xu, double yl, double yu, double lz, double uz);
    Void initialise(StringType x, StringType y, double xl, double xu, double yl, double yu);
    Void write(const char* filename) const;
    Void finalise();
    Void move_to(double x, double y);
    Void line_to(double x, double y);
    Void circle(double x, double y, double r);
    Void dot(double x, double y);
    Void stroke();
    Void fill();
    Void set_dot_radius(double dr);
    Void set_line_width(double lw);
    Void set_line_colour(double r, double g, double b);
    Void set_fill_opacity(double o);
    Void set_fill_colour(double r, double g, double b);
    Vector2d scaling() const;
    Box2d bounds() const;

    Void set_3D_palette();

    Void plotData(Array<double> data);
    Void plotBounds(Array<Array<double>> bounds);
    Void plotTensor2DImage(Tensor<2, double> tensor);
    Void plotTensor3DImage(Tensor<3, double> tensor);
    Void plotXYProjection(Tensor<3, double> tensor);
    Void plotYZProjection(Tensor<3, double> tensor);
    Void plotXZProjection(Tensor<3, double> tensor);
    //CanvasInterface

    //Set Multiplot - Multiple plot on same screen
    void setMultiplot(bool s);
    //Set Multiplot Layout
    void setMultiplotLayout(int nRow, int nCol, String Title);
    // Set X Label
    void setXLabel(String xLabel);
    // Set Y Label
    void setYLabel(String yLabel);
    // Set Z Label
    void setZLabel(String zLabel);
    // Set Title
    void setTitle(String title);
    // Set Labels
    void setXYZLabel(String xLabel, String yLabel, String zLabel);
    // Set Labels and Title
    void setLabels(String xLabel, String yLabel, String zLabel, String title);
    // Set X, Y Range
    void setRange2D(double minX, double maxX, double minY, double maxY);

    void setRange3D(double minX, double maxX, double minY,  double maxY, double minZ, double maxZ);
    // Set X Log axis
    void setXLogAxis();
    // Set Y Log axis
    void setYLogAxis();
    // Set XY Log axis
    void setXYLogAxis();
    // Set XZ Log axis
    void setXZLogAxis();
    // Set YZ Log axis
    void setYZLogAxis();
    // Set XYZ Log axis
    void setXYZLogAxis();
    // Set Legend
    void setLegend();
    // Set View Projection of a 3D rapresentation
    void setMap();
    //Set 3D palette
    void set3DPalette(double min, double max, double step);
    //Set 2D palette
    //void set2DPalette(Image2D& image, double min, double max, double step);
    void set2DPalette(double min, double max, double step);
    //Unset colorbox
    void unsetColorbox();

};

GnuplotCanvas::~GnuplotCanvas()
{
    delete gnuplot;
}

} // namespace Ariadne



#endif