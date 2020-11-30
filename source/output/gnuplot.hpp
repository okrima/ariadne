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

#include "utility/gnuplot-iostream.h"

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
    points2D,
    dots2D,
    steps2D,
    impulses2D,
};

struct _Line2D
{
    _LineStyle2D style = lines2D;
    //line weight
    int lw = 1;
    //line size
    int ls = 1;
};

const char *_linestyle2D[] = {"lines", "linespoints", "points", "dots", "steps", "impulses"};

enum _LineStyle3D // %s - es: ls, fs
{
    lines3D,
    surface3D,
    pm3d,
    points3D,
};

struct _Line3D
{
    _LineStyle3D style = lines3D;
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
    FloatDP Xmin = 0;
    FloatDP Xmax;
    FloatDP Ymin = 0;
    FloatDP Ymax;
};

struct _Range3D
{
    FloatDP Xmin = 0;
    FloatDP Xmax;
    FloatDP Ymin = 0;
    FloatDP Ymax;
    FloatDP Zmin = 0;
    FloatDP Zmax;
};

struct Image2D
{
    _Colours colour = _black;
    _Line2D linestyle2D;
    //_Range2D range2D;

};

struct Image3D
{
    _Colours colour = _black;
    _Line3D linestyle3D;
    //_Range3D range3D;
    
};

class GnuplotCanvas : public CanvasInterface
{
friend class Figure;
private:
    Gnuplot *gnuplot;
    List<Point2d> geom;
    Colour lc, fc;
    double lw;
    //_Range2D range2D;
    Nat dim;
    Point2d Cpoint;
    double dr;
    bool isdot;
    Nat sizeX;
    Nat sizeY;
    bool noCanvas;
    bool isMultiplot;
    bool is2DPalette;
    bool is3DPalette;

private:
    // Plot from Array
    void plot2D(Gnuplot& gp, Image2D& image, Array<double> data);
    // Plot Bounds from Array
    void plot2D(Gnuplot& gp, Image2D& image, Array<Array<double>> dataBound);
    // Plot from Array with Range
    void plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<double> data);
    // Plot from Array Bound with Range
    void plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<Array<double>> dataBound);
    // 3D Plot with tensor Array - Evolution in Time
    void plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Array<Array<double>> data);
    // 3D Plot Bounds with tensor Array - Evolution in Time
    void plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Array<Array<Array<double>>> dataBound);
    //Print 2D projection from 3D
    void XZProjection(Gnuplot& gp, Image2D& image, _Range2D& range, String filename);


public:
    ~GnuplotCanvas();
    // Constructors - Create the canvas
    //Create empty canvas
    GnuplotCanvas();
    //Create canvas with dimensions
    GnuplotCanvas(String filename, Nat X, Nat Y, GnuplotFileType fileType);

    //CanvasInterface
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
    //CanvasInterface


    //Set Multiplot - Multiple plot on same screen
    void setMultiplot(bool s);
    //Set Multiplot Layout
    void setMultiplotLayout(int nRow, int nCol, String Title);
    //Plot 2D Array
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, Array<T>& array)
    {
        //Create Array<double> and feed into plot2D w/ array
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(array.size(), 0);
        for (SizeType x = 0; x < array.size(); x++)
        {
            data[x] = array[x].get_d();
        }
        plot2D(gp, image, data);
    }
    //Plot 2D Bounds - Array
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, Array<Bounds<T>>& array)
    {
        //Create Array<double> and feed into plot2D w/ array
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";

        Array<Array<double>> dataBound(2);   //dataBound[0] = Lower Value
                                                //dataBound[1] = Upper Value
        for (SizeType i = 0; i < 2; i++)    //Resize each bound
        {
            dataBound[i].resize(array.size());
        }
        
        for (SizeType x = 0; x < array.size(); x++)
        {
            dataBound[0].at(x) = array[x].lower().get_d();
            dataBound[1].at(x) = array[x].upper().get_d();
        }
        plot2D(gp, image, dataBound);
    }
    //Plot 2D Vector
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, Vector<T>& array)
    {
        //Create Array<double> and feed into plot2D w/ array data
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(array.size(), 0);
        for (SizeType x = 0; x < array.size(); x++)
        {
            data[x] = array[x].get_d();
        }
        plot2D(gp, image, data);
    }
    //Plot 2D Bounds - Vector
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, Vector<Bounds<T>>& array)
    {
        //Create Array<double> and feed into plot2D w/ array
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";

        Array<Array<double>> dataBound(2); //dataBound[0] = Lower Value
                                                //dataBound[1] = Upper Value
        for (SizeType i = 0; i < 2; i++)    //Resize the bound value
        {
            dataBound[i].resize(array.size());
        }
        for (SizeType x = 0; x < array.size(); x++)
        {
            dataBound[0].at(x) = array[x].lower().get_d();
            dataBound[1].at(x) = array[x].upper().get_d();
        }
        plot2D(gp, image, dataBound);
    }

    //Plot 2D Array
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<T>& array)
    {
        //Create Array<double> and feed into plot2D w/ array data
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(array.size(), 0);
        for (SizeType x = 0; x < array.size(); x++)
        {
            data[x] = array[x].get_d();
        }
        plot2D(gp, image, range2D, data);
    }
    //Plot 2D Bounds - Array
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<Bounds<T>>& array)
    {
        //Create Array<double> and feed into plot2D w/ array
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";

        Array<Array<double>> dataBound(2); //dataBound[0] = Lower Value
                                                //dataBound[1] = Upper Value
        for (SizeType i = 0; i < 2; i++)
        {
            dataBound[i].resize(array.size());
        }
        for (SizeType x = 0; x < array.size(); x++)
        {
            dataBound[0].at(x) = array[x].lower().get_d();
            dataBound[1].at(x) = array[x].upper().get_d();
        }
        plot2D(gp, image, dataBound);
    }

    //Plot 2D Vector
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Vector<T>& array)
    {
        //Create Array<double> and feed into plot2D w/ array data
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(array.size(), 0);
        for (SizeType x = 0; x < array.size(); x++)
        {
            data[x] = array[x].get_d();
        }
        plot2D(gp, image, range2D, data);
    }
    //Plot 2D Bounds - Vector
    template<typename T>
    void plotArray2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Vector<Bounds<T>>& array)
    {
        //Create Array<double> and feed into plot2D w/ array
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";

        Array<Array<double>> dataBound(2); //dataBound[0] = Lower Value
                                                //dataBound[1] = Upper Value
        for (SizeType i = 0; i < 2; i++)
        {
            dataBound[i].resize(array.size());
        }
        for (SizeType x = 0; x < array.size(); x++)
        {
            dataBound[0].at(x) = array[x].lower().get_d();
            dataBound[1].at(x) = array[x].upper().get_d();
        }
        plot2D(gp, image, dataBound);
    }

    //Plot 2D Array with list
    template<typename T>
    void plotArray2D(Gnuplot& gp, Array<Image2D> imgList, Array<Array<T>> dataList)
    {
        ARIADNE_ASSERT(dataList.size() > 0);
        ARIADNE_ASSERT(imgList.size() > 0);
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(dataList[0].size(), 0);
        for (SizeType i = 0; i < dataList.size(); i++)
        {
            for (SizeType j = 0; j < data.size(); j++)
            {
                data[j] = dataList[i].at(j).get_d();
            }
            
            plot2D(gp, imgList[i], data);
        }
    }
    //Plot 2D Bounds Array with list
    template<typename T>
    void plotArray2D(Gnuplot& gp, Array<Image2D> imgList, Array<Array<Bounds<T>>> dataList)
    {
        ARIADNE_ASSERT(dataList.size() > 0);
        ARIADNE_ASSERT(imgList.size() > 0);
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        for (SizeType i = 0; i < dataList.size(); i++)  //For each Function
        {
            Array<Array<double>> dataBound(2);       //Create the bound value
            for (SizeType k = 0; k < 2; k++)
            {
                dataBound[i].resize(dataList[0].size());//Resize each bound array
            }
            
            for (SizeType j = 0; j < dataBound[0].size(); j++)
            {
                dataBound[0].at(j) = dataList[i].at(j).lower().get_d();
                dataBound[1].at(j) = dataList[i].at(j).upper().get_d();
            }
            plot2D(gp, imgList[i], dataBound);
        }        
    }

    //Plot 2D Vector with list
    template<typename T>
    void plotArray2D(Gnuplot& gp, Vector<Image2D> imgList, Vector<Vector<T>> dataList)
    {
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(dataList[0].size(), 0);
        for (SizeType i = 0; i < dataList.size(); i++)
        {
            for (SizeType j = 0; j < data.size(); j++)
            {
                data[j] = dataList[i].at(j).get_d();
            }
            
            plot2D(gp, imgList[i], data);
        }
    }
    //Plot 2D Bounds Vector with list
    template<typename T>
    void plotArray2D(Gnuplot& gp, Vector<Image2D> imgList, Vector<Vector<Bounds<T>>> dataList)
    {
        ARIADNE_ASSERT(dataList.size() > 0);
        ARIADNE_ASSERT(imgList.size() > 0);
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        for (SizeType i = 0; i < dataList.size(); i++)  //For each Function
        {
            Array<Array<double>> dataBound(2);       //Create the bound value
            for (SizeType k = 0; k < 2; k++)
            {
                dataBound[i].resize(dataList[0].size());//Resize each bound array
            }
            
            for (SizeType j = 0; j < dataBound[0].size(); j++)
            {
                dataBound[0].at(j) = dataList[i].at(j).lower().get_d();
                dataBound[1].at(j) = dataList[i].at(j).upper().get_d();
            }
            plot2D(gp, imgList[i], dataBound);
        } 
    }

    //Plot 2D Array with list and range
    template<typename T>
    void plotArray2D(Gnuplot& gp, Array<Image2D> imgList, _Range2D& range2D, Array<Array<T>> dataList)
    {
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(dataList[0].size(), 0);
        for (SizeType i = 0; i < dataList.size(); i++)
        {
            for (SizeType j = 0; j < data.size(); j++)
            {
                data[j] = dataList[i].at(j).get_d();
            }
            
            plot2D(gp, imgList[i], range2D, data);
        }
    }
    //Plot 2D Bounds Array with list and range
    template<typename T>
    void plotArray2D(Gnuplot& gp, Array<Image2D> imgList, _Range2D& range2D, Array<Array<Bounds<T>>> dataList)
    {
        ARIADNE_ASSERT(dataList.size() > 0);
        ARIADNE_ASSERT(imgList.size() > 0);
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        for (SizeType i = 0; i < dataList.size(); i++)  //For each Function
        {
            Array<Array<double>> dataBound(2);       //Create the bound value
            for (SizeType k = 0; k < 2; k++)
            {
                dataBound[i].resize(dataList[0].size());//Resize each bound array
            }
            
            for (SizeType j = 0; j < dataBound[0].size(); j++)
            {
                dataBound[0].at(j) = dataList[i].at(j).lower().get_d();
                dataBound[1].at(j) = dataList[i].at(j).upper().get_d();
            }
            plot2D(gp, imgList[i], dataBound);
        } 
    }

    //Plot 2D Vector with list
    template<typename T>
    void plotArray2D(Gnuplot& gp, Vector<Image2D> imgList, _Range2D& range2D, Vector<Vector<T>> dataList)
    {
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        Array<double> data(dataList[0].size(), 0);
        for (SizeType i = 0; i < dataList.size(); i++)
        {
            for (SizeType j = 0; j < data.size(); j++)
            {
                data[j] = dataList[i].at(j).get_d();
            }
            
            plot2D(gp, imgList[i], range2D, data);
        }
    }
    //Plot 2D Bounds Vector with list and range
    template<typename T>
    void plotArray2D(Gnuplot& gp, Vector<Image2D> imgList, _Range2D& range2D, Vector<Vector<Bounds<T>>> dataList)
    {
        ARIADNE_ASSERT(dataList.size() > 0);
        ARIADNE_ASSERT(imgList.size() > 0);
        ARIADNE_ASSERT(imgList.size() == dataList.size());
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        for (SizeType i = 0; i < dataList.size(); i++)  //For each Function
        {
            Array<Array<double>> dataBound(2);       //Create the bound value
            for (SizeType k = 0; k < 2; k++)
            {
                dataBound[i].resize(dataList[0].size());//Resize each bound array
            }
            
            for (SizeType j = 0; j < dataBound[0].size(); j++)
            {
                dataBound[0].at(j) = dataList[i].at(j).lower().get_d();
                dataBound[1].at(j) = dataList[i].at(j).upper().get_d();
            }
            plot2D(gp, imgList[i], dataBound);
        } 
    }

    //Plot 2D data from Tensor - time evolution
    template<typename T>
    void plotTensor2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Tensor<2, T>& tensor)
    {
        //Create Array<double> and feed into plot2D w/ array data
        Array<double> data(tensor.size(0), 0);
        std::cout << "\n\nPlotting 2D with GNUPLOT\n\n";
        for (SizeType step = 0; step < tensor.size(1); step++)
        {
            for (SizeType x = 0; x < data.size(); x++)
            {
                data[x] = (tensor[{x, step}]).get_d();
            }
            plot2D(gp, image, range2D, data);
        }
    }
    //Plot 2D Bounds data from Tensor - time evolution
    template<typename T>
    void plotTensor2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Tensor<2, Bounds<T>>& tensor)
    {
        //Create Array<double> and feed into plot2D w/ array data
        Array<Array<double>> dataBound(2);
        for (SizeType step = 0; step < tensor.size(1) ; step++) //For each time step
        {
            for (SizeType i = 0; i < 2; i++)
            {
                dataBound[i].resize(tensor.size(0));
            }
            
            for (SizeType x = 0; x < tensor.size(0); x++)
            {
                dataBound[0] = tensor[{x, step}].lower().get_d();
                dataBound[1] = tensor[{x, step}].upper().get_d();
            }
            plot2D(gp, image, range2D, dataBound);
        }
        
    }

    // Plot 3D Array from Tensor - time evolution
    template<typename T>
    void plotTensor3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Tensor<3, T>& tensor)
    {
        //Create Array<Array<double>> and send into plot3D w/ array data
        SizeType dimX = tensor.size(0);
        SizeType dimY = tensor.size(1);
        SizeType dimTime = tensor.size(2);
        Array<Array<double>> data(dimX);

        std::cout << "\n\nPlotting 3D with GNUPLOT\n\n";

        for(SizeType step = 0; step < dimTime; step++)
        {
            for (SizeType i = 0; i < dimY; i++)
            {
                data[i].resize(dimX);
                for (SizeType j = 0; j < dimX; j++)
                {
                    data[i].at(j) = tensor[{j, i, step}].get_d();
                }    
            }
            plot3D(gp, image, range3D, data);
        }
    }
    // Plot 3D Bounds Array from Tensor - time evolution

    template<typename T>
    void plotTensor3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Tensor<3, Bounds<T>>& tensor)
    {
        //Create Array<Array<double>> and send into plot3D w/ array data
        SizeType dimX = tensor.size(0);
        SizeType dimY = tensor.size(1);
        SizeType dimTime = tensor.size(2);
        

        std::cout << "\n\nPlotting 3D with GNUPLOT\n\n";

        for (SizeType step = 0; step < dimTime; step++)
        {
            Array<Array<Array<double>>> dataBound(2);
            for (SizeType i = 0; i < 2; i++)    //Resize
            {
                dataBound[i].resize(dimX);
            }
            
            for (SizeType x = 0; x < dimX; x++)
            {
                dataBound[0].at(x).resize(dimY);
                dataBound[1].at(x).resize(dimY);
                for (SizeType y = 0; y < dimY; y++)
                {
                    dataBound[0].at(x).at(y) = tensor[{x, y, step}].lower().get_d();
                    dataBound[1].at(x).at(y) = tensor[{x, y, step}].upper().get_d();
                }  
            }
            plot3D(gp, image,range3D, dataBound); 
        }  
    }

    //Plot XZ projection
    template<typename T>
    void plotXZProjection(Gnuplot& gp, Image3D& image, _Range3D& range, Tensor<3, T>& tensor)
    {
        _Line2D line;
        line.style = lines2D;
        Image2D img;
        img.linestyle2D = line;
        _Range2D rng;
        rng.Xmin = range.Xmin;
        rng.Xmax = range.Xmax;
        rng.Ymin = range.Zmin;
        rng.Ymax = range.Zmax;
        set2DPalette(img, rng.Ymin, rng.Ymax, 0.2);
        Array<double> data(tensor.size(0), 0);
        if (!isMultiplot)
        { 
            for (SizeType step = 0; step < tensor.size(2); step++)
            { 
                setMultiplot(true);
                for (SizeType y = 0; y < tensor.size(1); y++)
                {
                    for (SizeType x = 0; x < tensor.size(0); x++)
                    {
                        data[x] = tensor[{x, y, step}].get_d();
                    }
                    plot2D(gp, img, rng, data);    
                }
            }
            setMultiplot(false);
        }
        else{
            ARIADNE_ERROR("Impossible to plot 3D image or projections with multiplot, please set multiplot = FALSE");
        }
    }

    //Plot XZ projection
    template<typename T>
    void plotYZProjection(Gnuplot& gp, Image3D& image, _Range3D& range, Tensor<3, T>& tensor)
    {
        _Line2D line;
        line.style = lines2D;
        Image2D img;
        img.linestyle2D = line;
        _Range2D rng;
        rng.Xmin = range.Ymin;
        rng.Xmax = range.Ymax;
        rng.Ymin = range.Zmin;
        rng.Ymax = range.Zmax;
        set2DPalette(img, rng.Ymin, rng.Ymax, 0.2);
        Array<double> data(tensor.size(1), 0);
        if (!isMultiplot)
        { 
            for (SizeType step = 0; step < tensor.size(2); step++)
            { 
                setMultiplot(true);
                for (SizeType x = 0; x < tensor.size(0); x++)
                {
                    for (SizeType y = 0; y < tensor.size(1); y++)
                    {
                        data[y] = tensor[{x, y, step}].get_d();
                    }
                    plot2D(gp, img, rng, data);
                    //setMultiplot(gp, true);     
                }
            }
            setMultiplot(false);
        }
        else{
            ARIADNE_ERROR("Impossible to plot 3D image or projections with multiplot, please set multiplot = FALSE");
        }
    }
    // Set Terminal output
    void setTerminal(_Format format, String nameFile);
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
    void setRange2D(FloatDP minX, FloatDP maxX, FloatDP minY, FloatDP maxY);
    // Set X, Y specular symmetric range
    void setRange2D(_Range2D& range2D, FloatDP maxX, FloatDP maxY);
    // Set X, Y range
    void setRange2D(_Range2D& range2D, FloatDP minX, FloatDP maxX, 
                FloatDP minY, FloatDP maxY);
    // Set X, Y, Z range
    void setRange3D(_Range3D& range3D, FloatDP minX, FloatDP maxX, 
                FloatDP minY, FloatDP maxY,
                FloatDP minZ, FloatDP maxZ);
    void setRange3D(_Range3D& range3D, FloatDP maxX, FloatDP maxY, FloatDP maxZ);
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
    void set3DPalette(Image3D& image, FloatDP min, FloatDP max, FloatDP step, bool s);
    //Set 2D palette
    void set2DPalette(Image2D& image, FloatDP min, FloatDP max, FloatDP step);
    //Unset colorbox
    void unsetColorbox();
    //Set plane projection
    void setXYprojection();

};

GnuplotCanvas::~GnuplotCanvas()
{
}

} // namespace Ariadne



#endif