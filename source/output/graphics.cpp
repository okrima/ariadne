/***************************************************************************
 *            output/graphics.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "../utility/standard.hpp"
#include "../config.hpp"

#include "../utility/macros.hpp"
#include "../utility/stlio.hpp"
#include "../numeric/numeric.hpp"
#include "../function/function.hpp"
#include "../geometry/point.hpp"
#include "../geometry/box.hpp"
#include "../output/geometry2d.hpp"
#include "../output/graphics.hpp"
#include "../output/cairo.hpp"

#include "../output/gnuplot.hpp"

namespace Ariadne {

static const Int DEFAULT_WIDTH = 800;
static const Int DEFAULT_HEIGHT = 800;

#ifdef HAVE_CAIRO_H

static const Int LEFT_MARGIN = 160;
static const Int BOTTOM_MARGIN = 40;
static const Int TOP_MARGIN = 10;
static const Int RIGHT_MARGIN = 10;

#endif

OutputStream& operator<<(OutputStream& os, const DrawableInterface& drawable) {
    if(auto writable=dynamic_cast<const WritableInterface*>(&drawable)) {
        os << *writable;
    } else {
        os << "Drawable";
    }
    return os;
}



Colour::Colour()
    : Colour("transparant", 1.0, 1.0, 1.0, 0.0) { }
Colour::Colour(Dbl rd, Dbl gr, Dbl bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(Dbl rd, Dbl gr, Dbl bl, Dbl op)
    : Colour("",rd,gr,bl,op) { }
Colour::Colour(const Char* nm, Dbl rd, Dbl gr, Dbl bl, Bool tr)
    : Colour("",rd,gr,bl,tr?0.0:1.0) { }
Colour::Colour(const char* nm, Dbl rd, Dbl gr, Dbl bl, Dbl op)
    : name(nm), red(rd), green(gr), blue(bl), opacity(op) { }
OutputStream& operator<<(OutputStream& os, const Colour& c) {
    return os << "Colour( name=" << c.name << ", r=" << c.red << ", g=" << c.green << ", b=" << c.blue << ", op=" << c.opacity << " )"; }


const Colour transparant=Colour();

const Colour white=Colour("white",1.0,1.0,1.0);
const Colour black=Colour("black",0.0,0.0,0.0);
const Colour red=Colour("red",1.0,0.0,0.0);
const Colour green=Colour("green",0.0,1.0,0.0);
const Colour blue=Colour("blue",0.0,0.0,1.0);
const Colour yellow=Colour("yellow",1.0,1.0,0.0);
const Colour cyan=Colour("cyan",0.0,1.0,1.0);
const Colour magenta=Colour("magenta",1.0,0.0,1.0);
const Colour ariadneorange=Colour("ariadneorange",1.0,0.75,0.5);

OutputStream& operator<<(OutputStream& os, GraphicsProperties const& gp) {
    return os << "GraphicsProperties(" << "dot_radius=" << gp.dot_radius
        << ", line_style=" << gp.line_style<<", line_width=" << gp.line_width << ", line_colour=" << gp.line_colour
        << ", fill_style=" << gp.fill_style << ", fill_colour=" << gp.fill_colour << ")"; }


inline StringType str(FloatDP x) {
    StringStream ss;
    ss << x;
    return ss.str();
}

Void draw(Figure& fig, const DrawableInterface& shape) {
    fig.draw(shape);
}

Void draw(Figure& fig, FloatDPApproximateBox const& box) {
    fig.draw(box);
}

struct GraphicsObject {
    GraphicsObject(const GraphicsProperties& gp, const DrawableInterface& sh)
        : properties(gp), shape_ptr(sh.clone()) { }
    GraphicsProperties properties;
    std::shared_ptr<const DrawableInterface> shape_ptr;
};



struct Figure::Data
{
    Data() : bounding_box(0), projection(2,0,1), properties() { }
    Data(ApproximateBoxType bbx, PlanarProjectionMap prj) : bounding_box(bbx), projection(prj), properties() { }
    ApproximateBoxType bounding_box;
    PlanarProjectionMap projection;
    GraphicsProperties properties;
    std::vector<GraphicsObject> objects;
};

Figure::~Figure()
{
    delete this->_data;
}


Figure::Figure()
    : Figure(ApproximateBoxType({{-1,+1},{-1,+1}}),PlanarProjectionMap(2,0,1))
{
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, const PlanarProjectionMap& proj)
    : _data(new Data(bbx,proj))
{
    ARIADNE_ASSERT_MSG(proj.argument_size() == bbx.dimension(), "Coordinate projection "<<proj<<" must take same number of arguments as the dimension of the bounding box "<<bbx);
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, DimensionType ix, DimensionType iy)
    : Figure(bbx,PlanarProjectionMap(bbx.dimension(),ix,iy))
{
}

Figure::Figure(const GraphicsBoundingBoxType& bbx, Pair<DimensionType,DimensionType> ixy)
    : Figure(bbx,PlanarProjectionMap(bbx.dimension(),ixy.first,ixy.second))
{
}

Figure& Figure::draw(const DrawableInterface& shape)
{
    this->_data->objects.push_back(GraphicsObject(this->_data->properties,shape)); return *this;
}


Figure& Figure::set_projection(DimensionType as, DimensionType ix, DimensionType iy)
{
    this->_data->projection=PlanarProjectionMap(as,ix,iy); return *this;
}

Figure& Figure::set_projection_map(const PlanarProjectionMap& p)
{
    this->_data->projection=p; return *this;
}

Figure& Figure::set_bounding_box(const ApproximateBoxType& bx)
{
    this->_data->bounding_box=bx; return *this;
}

PlanarProjectionMap Figure::get_projection_map() const
{
    return this->_data->projection;
}

ApproximateBoxType Figure::get_bounding_box() const
{
    return this->_data->bounding_box;
}

Figure& Figure::set_dot_radius(Dbl dr)
{
    this->_data->properties.dot_radius=dr; return *this;
}

Figure& Figure::set_line_style(Bool ls)
{
    this->_data->properties.line_style=ls; return *this;
}

Figure& Figure::set_line_width(Dbl lw)
{
    this->_data->properties.line_width=lw; return *this;
}

Figure& Figure::set_line_colour(Colour lc)
{
    this->_data->properties.line_colour=lc; return *this;
}

Figure& Figure::set_line_colour(Dbl r, Dbl g, Dbl b)
{
    this->set_line_colour(Colour(r,g,b)); return *this;
}

Figure& Figure::set_fill_style(Bool fs)
{
    this->_data->properties.fill_style=fs; return *this;
}

Figure& Figure::set_fill_opacity(Dbl fo)
{
    this->_data->properties.fill_colour.opacity=fo; return *this;
}

Figure& Figure::set_fill_colour(Colour fc)
{
    this->_data->properties.fill_colour=fc; return *this;
}

Figure& Figure::set_fill_colour(Dbl r, Dbl g, Dbl b)
{
    this->set_fill_colour(Colour(r,g,b,this->_data->properties.fill_colour.opacity)); return *this;
}

Bool Figure::get_line_style() const
{
    return this->_data->properties.line_style;
}

Dbl Figure::get_line_width() const
{
    return this->_data->properties.line_width;
}

Dbl Figure::get_dot_radius() const
{
    return this->_data->properties.dot_radius;
}

Colour Figure::get_line_colour() const
{
    return this->_data->properties.line_colour;
}


Bool Figure::get_fill_style() const
{
    return this->_data->properties.fill_style;
}

Dbl Figure::get_fill_opacity() const
{
    return this->_data->properties.fill_colour.opacity;
}

Colour Figure::get_fill_colour() const
{
    return this->_data->properties.fill_colour;
}

Figure& Figure::draw(ApproximateBoxType const& box)
{
    ApproximateBoxSetType box_set(box);
    DrawableInterface const& shape=box_set;
    this->draw(shape); return *this;
}

Figure& Figure::draw(RealBox const& box)
{
    return this->draw(ApproximateBoxType(box,DoublePrecision()));
}

Figure& Figure::clear() {
    this->_data->objects.clear(); return *this;
}

Void set_properties(CanvasInterface& canvas, const GraphicsProperties& properties) {
    const Colour& line_colour=properties.line_colour;
    const Colour& fill_colour=properties.fill_colour;
    if (properties.line_style==false) { canvas.set_line_width(0); }
    else { canvas.set_line_width(properties.line_width); }
    if (properties.fill_style==false) { canvas.set_fill_opacity(0); }
    else { canvas.set_fill_opacity(properties.fill_colour.opacity); }
    canvas.set_fill_colour(fill_colour.red, fill_colour.green, fill_colour.blue);
    canvas.set_line_colour(line_colour.red, line_colour.green, line_colour.blue);
}

Void Figure::_paint_all(CanvasInterface& canvas) const
{
    ApproximateBoxType bounding_box=this->_data->bounding_box;
    const PlanarProjectionMap projection=this->_data->projection;
    const std::vector<GraphicsObject>& objects=this->_data->objects;

    DimensionType dimension=projection.argument_size();

    // Don't attempt to compute a bounding box, as this relies on
    // a drawable object having one. Instead, the bounding box must be
    // specified explicitly
    if(bounding_box.dimension()==0) {
        bounding_box=ExactBoxType(dimension,ExactIntervalType(-1,1));
    }

    // Check projection and bounding box have same values.
    ARIADNE_ASSERT_MSG(bounding_box.dimension()==projection.argument_size(),"bounding_box="<<bounding_box<<", projection="<<projection);
    ARIADNE_ASSERT(bounding_box.dimension()>projection.x_coordinate());
    ARIADNE_ASSERT(bounding_box.dimension()>projection.y_coordinate());

    // Project the bounding box onto the canvas
    Dbl xl=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].lower());
    Dbl xu=numeric_cast<Dbl>(bounding_box[projection.x_coordinate()].upper());
    Dbl yl=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].lower());
    Dbl yu=numeric_cast<Dbl>(bounding_box[projection.y_coordinate()].upper());

    StringType tx=StringType("x")+str(projection.x_coordinate());
    StringType ty=StringType("x")+str(projection.y_coordinate());
    canvas.initialise(tx,ty,xl,xu,yl,yu);

    // Draw shapes
    for(SizeType i=0; i!=objects.size(); ++i) {
        const DrawableInterface& shape=objects[i].shape_ptr.operator*();
        if(shape.dimension()==0) { break; } // The dimension may be equal to two for certain empty sets.
        ARIADNE_ASSERT_MSG(dimension==shape.dimension(),
                           "Shape "<<shape<<", dimension="<<shape.dimension()<<", bounding_box="<<bounding_box);
        set_properties(canvas, objects[i].properties);
        shape.draw(canvas,this->_data->projection);
    }

    canvas.finalise();
}


Void
Figure::write(const char* cfilename) const
{
    this->write(cfilename, DEFAULT_WIDTH, DEFAULT_HEIGHT);
}


Void
Figure::write(const char* cfilename, Nat drawing_width, Nat drawing_height) const
{
    SharedPointer<CanvasInterface> canvas=make_canvas(drawing_width,drawing_height);

    this->_paint_all(*canvas);

    StringType filename(cfilename);
    if(filename.rfind(".") != StringType::npos) {
    } else {
        filename=filename+".png";
    }

    canvas->write(filename.c_str());
}



#ifdef HAVE_CAIRO_H

SharedPointer<CanvasInterface> make_canvas(Nat drawing_width, Nat drawing_height) {
    return std::make_shared<CairoCanvas>(ImageSize2d(drawing_width,drawing_height));
}

#else

SharedPointer<CanvasInterface> make_canvas(Nat drawing_width, Nat drawing_height) {
    ARIADNE_WARN_ONCE("No facilities for displaying graphics are available.");
    return std::make_shared<NullCanvas>();
}

#endif


#ifdef HAVE_CAIRO_H

CairoCanvas::~CairoCanvas()
{
    cairo_surface_destroy(cairo_get_target(cr));
    cairo_destroy(cr);
}

CairoCanvas::CairoCanvas(cairo_t *c)
    : cr(c), lw(1.0), dr(1.0), lc(0.0,0.0,0.0), fc(1.0,1.0,1.0, 1.0)
{
}

CairoCanvas::CairoCanvas(const ImageSize2d& size)
    : cr(0), lw(1.0), dr(1.0), lc(0.0,0.0,0.0), fc(1.0,1.0,1.0, 1.0)
{
    const Int canvas_width = static_cast<Int>(size.nx)+LEFT_MARGIN+RIGHT_MARGIN;
    const Int canvas_height = static_cast<Int>(size.ny)+BOTTOM_MARGIN+TOP_MARGIN;

    cairo_surface_t* surface = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, canvas_width, canvas_height);
    cr = cairo_create (surface);
}

CairoCanvas::CairoCanvas(const ImageSize2d& size, const Box2d& bounds)
    : CairoCanvas(size)
{
}

Vector2d CairoCanvas::scaling() const
{
    ImageSize2d sz = this->size_in_pixels();
    Box2d bb=this->bounds();
    return Vector2d((bb.xu-bb.xl)/sz.nx,(bb.yu-bb.yl)/sz.ny);
}

Box2d CairoCanvas::bounds() const
{
    double xl=LEFT_MARGIN;
    double yu=TOP_MARGIN;
    double xu=cairo_image_surface_get_width(cairo_get_target(cr))-RIGHT_MARGIN;
    double yl=cairo_image_surface_get_height(cairo_get_target(cr))-BOTTOM_MARGIN;
    cairo_device_to_user(cr,&xl,&yu);
    cairo_device_to_user(cr,&xu,&yl);
    return Box2d(xl,xu,yl,yu);
}

Void CairoCanvas::stroke()
{
    cairo_save(cr);

    // Set user and device space identical so that the line width is interpreted as pixels
    cairo_identity_matrix(cr);

    cairo_set_source_rgb(cr, lc.red,lc.green,lc.blue);
    cairo_set_line_width(cr, lw);
    cairo_stroke (this->cr);

    cairo_restore(cr);
}

Void CairoCanvas::fill() {
    cairo_set_source_rgba(cr,fc.red,fc.green,fc.blue,fc.opacity);
    cairo_fill_preserve (cr);
    this->stroke();
}

ImageSize2d CairoCanvas::size_in_pixels() const {
    return ImageSize2d(cairo_image_surface_get_width(cairo_get_target(cr))-(LEFT_MARGIN+RIGHT_MARGIN),
                       cairo_image_surface_get_height(cairo_get_target(cr))-(BOTTOM_MARGIN+TOP_MARGIN));
}

Void CairoCanvas::move_to(double x, double y) { cairo_move_to (cr, x, y); }
Void CairoCanvas::line_to(double x, double y) { cairo_line_to (cr, x, y); }
Void CairoCanvas::circle(double x, double y, double r) { cairo_arc (cr, x, y, r, 0, 2*M_PI); }
Void CairoCanvas::dot(double x, double y) { cairo_arc (cr, x, y, dr/1000, 0, 2*M_PI); }
Void CairoCanvas::set_dot_radius(double radius) { this->dr=radius; }
Void CairoCanvas::set_line_width(double width) { this->lw=width; }
Void CairoCanvas::set_line_colour(double r, double g, double b) { lc.red=r; lc.green=g; lc.blue=b; }
Void CairoCanvas::set_fill_opacity(double o) { fc.opacity=o; }
Void CairoCanvas::set_fill_colour(double r, double g, double b) { fc.red=r; fc.green=g; fc.blue=b; }


// TODO: Use generic canvas routines; move cairo-specific functionality
// into CairoCanvas class.
Void CairoCanvas::initialise(StringType text_x, StringType text_y, double xl, double xu, double yl, double yu)
{


    CairoCanvas& cairo_canvas=*this;
    cairo_t *crp=cairo_canvas.cr;

    const ImageSize2d drawing_size = cairo_canvas.size_in_pixels();
    const Int drawing_width = static_cast<Int>(drawing_size.nx);
    const Int drawing_height = static_cast<Int>(drawing_size.ny);

    //const Int canvas_width = cairo_image_surface_get_width(cairo_get_target(cr));
    //const Int canvas_height = cairo_image_surface_get_height(cairo_get_target(cr));

    const Int left_margin = LEFT_MARGIN;
    //const Int right_margin = RIGHT_MARGIN;
    //const Int bottom_margin = BOTTOM_MARGIN;
    const Int top_margin = TOP_MARGIN;

    // clear background
    cairo_set_source_rgb (crp, 1,1,1);
    cairo_paint (crp);

    // Set text font
    cairo_select_font_face (crp, "roman",CAIRO_FONT_SLANT_NORMAL,CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size (crp, 30);

    // Set text colour
    cairo_set_source_rgb (crp, 0., 0., 0.);

    // Get axis label text
    StringType text_xl=str(xl);
    StringType text_xu=str(xu);
    StringType text_yl=str(yl);
    StringType text_yu=str(yu);

    // Write axis labels
    cairo_text_extents_t te;
    cairo_text_extents (crp, text_xl.c_str(), &te);
    cairo_move_to(crp, left_margin-2, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_xl.c_str());
    cairo_text_extents (crp, text_xu.c_str(), &te);
    cairo_move_to(crp, left_margin+drawing_width-te.width-4, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_xu.c_str());
    cairo_text_extents (crp, text_x.c_str(), &te);
    cairo_move_to(crp, left_margin+drawing_width/2-te.width/2-3, top_margin+drawing_height+4+te.height);
    cairo_show_text (crp, text_x.c_str());

    cairo_text_extents (crp, text_yl.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+drawing_height+2);
    cairo_show_text (crp, text_yl.c_str());
    cairo_text_extents (crp, text_yu.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+te.height+2);
    cairo_show_text (crp, text_yu.c_str());
    cairo_text_extents (crp, text_y.c_str(), &te);
    cairo_move_to(crp, left_margin-te.width-6, top_margin+drawing_height/2+te.height+2);
    cairo_show_text (crp, text_y.c_str());


    // Save unclipped state and canvas coordinates
    cairo_save (crp);

    // Set clipping region
    cairo_move_to (crp, left_margin, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin);
    cairo_line_to (crp, left_margin, top_margin);
    cairo_line_to (crp, left_margin, top_margin+drawing_height);

    // Fill clipping region with a very light colour to indicate where figure
    // should be drawn
    cairo_set_source_rgb(crp, 0.97,0.97,0.97);
    cairo_fill_preserve (crp);
    cairo_clip (crp);
    cairo_new_path (crp);

    //std::cerr<<"cw="<<canvas_width<<" lm="<<left_margin<<" dw="<<drawing_width<<" rm="<<right_margin<<" xl="<<xl<<" xu="<<xu<<"\n";
    //std::cerr<<"ch="<<canvas_height<<"tm="<<top_margin<<" dw="<<drawing_height<<" bm="<<bottom_margin<<" yl="<<yl<<" yu="<<yu<<"\n";

    // compute device to user coordinate transformation
    double ctr0=left_margin;
    double ctr1=top_margin;
    double sc0=(drawing_width)/(xu-xl);
    double sc1=-(drawing_height)/(yu-yl);
    double utr0=(-xl);
    double utr1=(-yu);

    // Scale to user coordinates
    cairo_translate(crp, ctr0, ctr1);
    cairo_scale (crp, sc0,sc1);
    cairo_translate(crp, utr0, utr1);
}

Void CairoCanvas::write(const char* filename) const {
    cairo_surface_t* surface = cairo_get_target (cr);
    cairo_surface_write_to_png (surface, filename);
}

Void CairoCanvas::finalise()
{
    cairo_t *crp=this->cr;

    // Restore canvas coordinates and unclipped state
    cairo_restore (crp);

    const ImageSize2d drawing_size = this->size_in_pixels();
    const Int drawing_width = static_cast<Int>(drawing_size.nx);
    const Int drawing_height = static_cast<Int>(drawing_size.ny);

    const Int left_margin = LEFT_MARGIN;
    const Int top_margin = TOP_MARGIN;

    cairo_set_line_width (crp, 2.0);
    cairo_set_source_rgb (crp, 0.0, 0.0, 0.0);
    cairo_move_to (crp, left_margin, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin+drawing_height);
    cairo_line_to (crp, left_margin+drawing_width, top_margin);
    cairo_line_to (crp, left_margin, top_margin);
    cairo_line_to (crp, left_margin, top_margin+drawing_height);
    cairo_stroke (crp);
}

#endif

#ifdef HAVE_GNUPLOT_H

//Create the canvas
GnuplotCanvas::GnuplotCanvas()
{
    noCanvas = true;
    isMultiplot = false;
    is3DPalette = false;

}

GnuplotCanvas::GnuplotCanvas(Image2D& image, int X, int Y)
{
    ARIADNE_ASSERT(sizeX >= 0);
    ARIADNE_ASSERT(sizeY >= 0);
    noCanvas = false;
    isMultiplot = false;
    sizeX = X;
    sizeY = Y;
}

GnuplotCanvas::GnuplotCanvas(Image3D& image, int X, int Y)
{
    ARIADNE_ASSERT(sizeX >= 0);
    ARIADNE_ASSERT(sizeY >= 0);
    noCanvas = false;
    isMultiplot = false;
    is3DPalette = false;
    sizeX = X;
    sizeY = Y;
}

void GnuplotCanvas::setMultiplot(Gnuplot& gp, bool s)
{
    if (s)
    {
        isMultiplot = true;
        gp << "set multiplot\n";
    }  
}

void GnuplotCanvas::setMultiplotLayout(Gnuplot& gp, int nRow, int nCol, String title)
{
    gp << "set multiplot layout " << to_string(nRow) << "," << to_string(nCol) << " title \"" << title << "\"\n";
}

void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, Array<double> data)
{
    // START PLOT SINTAX
    gp << "plot ";
    gp << "[] ";

    // Gnuplot wait for input
    gp << "'-' ";
    // set linestyle
    gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
        " lw " << to_string(image.linestyle2D.lw);
    // set colour
    gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    //Send data through pipe gp
    gp.send1d(data);
}
/*
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, Vector<double> data)
{
    // START PLOT SINTAX
    gp << "plot ";
    gp << "[] ";

    // Gnuplot wait for input
    gp << "'-' ";
    // set linestyle
    gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
        " lw " << to_string(image.linestyle2D.lw);
    // set colour
    gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    //Send data through pipe gp
    gp.send1d(data);
}
*/
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, Array<Array<double>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)    //For each value - Upper and Lower
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);      
        
    }
    else
    {
        for (SizeType i = 0; i < 2; i++)    //Lower = 0 and Upper = 1 value
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
        }
    }
}
/*
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, Vector<Vector<double>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);         
    }
    else
    {
        for (SizeType i = 0; i < 2; i++)    //Lower = 0 and Upper = 1 value
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
        }
    }   
}
*/
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<double> data)
{
    // START PLOT SINTAX
    gp << "plot ";
    // set Range
    if (!to_string(range2D.Xmax).empty()){
        gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
    }else
    {
        gp << "[] ";
    }
    if (!to_string(range2D.Ymax).empty()){
        gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
    }
    else
    {
        gp << "[] ";
    }

    // Gnuplot wait for input
    gp << "'-' ";
    // set linestyle
    gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
        " lw " << to_string(image.linestyle2D.lw);
    // set colour
    gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    //Send data through pipe gp
    gp.send1d(data);
}

void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Array<Array<double>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)
        {
            // START PLOT SINTAX
            gp << "plot ";
            if (!to_string(range2D.Xmax).empty()){
                gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
            }else
            {
                gp << "[] ";
            }
            if (!to_string(range2D.Ymax).empty()){
                gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
            }

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);         
    }
    else
    {
        for (SizeType i = 0; i < 2; i++)    //Lower = 0 and Upper = 1 value
        {
            // START PLOT SINTAX
            gp << "plot ";
            if (!to_string(range2D.Xmax).empty()){
                gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
            }else
            {
                gp << "[] ";
            }
            if (!to_string(range2D.Ymax).empty()){
                gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
            }

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
        }
    }      
}
/*
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Vector<double> data)
{
    // START PLOT SINTAX
    gp << "plot ";
    // set Range
    if (!to_string(range2D.Xmax).empty()){
        gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
    }else
    {
        gp << "[] ";
    }
    if (!to_string(range2D.Ymax).empty()){
        gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
    }
    else
    {
        gp << "[] ";
    }

    // Gnuplot wait for input
    gp << "'-' ";
    // set linestyle
    gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
        " lw " << to_string(image.linestyle2D.lw);
    // set colour
    gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    //Send data through pipe gp
    gp.send1d(data);
}
*/
/*
void GnuplotCanvas::plot2D(Gnuplot& gp, Image2D& image, _Range2D& range2D, Vector<Vector<double>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)
        {
            // START PLOT SINTAX
            gp << "plot ";
            if (!to_string(range2D.Xmax).empty()){
                gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
            }else
            {
                gp << "[] ";
            }
            if (!to_string(range2D.Ymax).empty()){
                gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
            }

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);         
    }
    else
    {
        for (SizeType i = 0; i < 2; i++)    //Lower = 0 and Upper = 1 value
        {
            // START PLOT SINTAX
            gp << "plot ";
            if (!to_string(range2D.Xmax).empty()){
                gp << "[" << to_string(range2D.Xmin) << ":" << to_string(range2D.Xmax-1) << "] ";
            }else
            {
                gp << "[] ";
            }
            if (!to_string(range2D.Ymax).empty()){
                gp << "[" << to_string(range2D.Ymin) << ":" << to_string(range2D.Ymax) << "] ";
            }

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle2D.style] << " ls " << to_string(image.linestyle2D.ls) <<
                " lw " << to_string(image.linestyle2D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send1d(dataBound[i]);
        }
    }      
}

*/
void GnuplotCanvas::plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Array<Array<double>> data)
{  
    // START PLOT SINTAX
    gp << "splot ";
    // get Range
    if (!to_string(range3D.Xmax).empty()){
        gp << "[" << to_string(range3D.Xmin) << ":" << to_string(range3D.Xmax) << "] ";
    }else
    {
        gp << "[ ] ";
    }
    if (!to_string(range3D.Ymax).empty()){
        gp << "[" << to_string(range3D.Ymin) << ":" << to_string(range3D.Ymax) << "] ";
    }
    else
    {
        gp << "[ ] ";
    }
    if (!to_string(range3D.Zmax).empty()){
        gp << "[" << to_string(range3D.Zmin) << ":" << to_string(range3D.Zmax) << "] ";
    }
    else
    {
        gp << "[ ] ";
    }
    // Gnuplot wait for input 
    gp << "'-' ";
    // get linestyle
    gp << "with " << _linestyle3D[image.linestyle3D.style] << " ";

    if (image.linestyle3D.style != surface3D && image.linestyle3D.style != pm3d)
    {
        gp << " ls " << to_string(image.linestyle3D.ls) << " lw " << to_string(image.linestyle3D.lw);
    }
    // get colour
    if (!is3DPalette)
    {
        gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    }
    else
    {
        gp<< "\n";
    }
    //Send data through pipe gp
    gp.send2d(data);
}
/*
void GnuplotCanvas::plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Vector<Vector<double>> data)
{
    // START PLOT SINTAX
    gp << "splot ";
    // get Range
    if (!to_string(range3D.Xmax).empty()){
        gp << "[" << to_string(range3D.Xmin) << ":" << to_string(range3D.Xmax) << "] ";
    }else
    {
        gp << "[ ] ";
    }
    if (!to_string(range3D.Ymax).empty()){
        gp << "[" << to_string(range3D.Ymin) << ":" << to_string(range3D.Ymax) << "] ";
    }
    else
    {
        gp << "[ ] ";
    }
    if (!to_string(range3D.Zmax).empty()){
        gp << "[" << to_string(range3D.Zmin) << ":" << to_string(range3D.Zmax) << "] ";
    }
    else
    {
        gp << "[ ] ";
    }
    // Gnuplot wait for input 
    gp << "'-' ";
    // get linestyle
    gp << "with " << _linestyle3D[image.linestyle3D.style] << " ";

    if (image.linestyle3D.style != surface3D && image.linestyle3D.style != pm3d)
    {
        gp << " ls " << to_string(image.linestyle3D.ls) << " lw " << to_string(image.linestyle3D.lw);
    }
    // get colour
    if (!is3DPalette)
    {
        gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
    }
    else
    {
        gp<< "\n";
    }
    //Send data through pipe gp
    gp.send2d(data);
} 
*/
/*
void GnuplotCanvas::plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Vector<Vector<Vector<double>>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)    //Lower and Upper Value
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle3D.style] << " ls " << to_string(image.linestyle3D.ls) <<
                " lw " << to_string(image.linestyle3D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send2d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);         
    }
    else{
        ARIADNE_ERROR("Impossible to plot 3D image with multiplot, please set multiplot = FALSE");
    }
} 
*/
void GnuplotCanvas::plot3D(Gnuplot& gp, Image3D& image, _Range3D& range3D, Array<Array<Array<double>>> dataBound)
{
    if (!isMultiplot)
    {
        setMultiplot(gp, false);
        for (SizeType i = 0; i < 2; i++)    //Lower and Upper Value
        {
            // START PLOT SINTAX
            gp << "plot ";
            gp << "[] ";

            // Gnuplot wait for input
            gp << "'-' ";
            // set linestyle
            gp << "with " << _linestyle2D[image.linestyle3D.style] << " ls " << to_string(image.linestyle3D.ls) <<
                " lw " << to_string(image.linestyle3D.lw);
            // set colour
            gp << " lc rgb \"" << _colours[image.colour] << "\"\n";
            //Send data through pipe gp
            gp.send2d(dataBound[i]);
            setMultiplot(gp, true);
        }
        setMultiplot(gp, false);         
    }
    else{
        ARIADNE_ERROR("Impossible to plot 3D image with multiplot, please set multiplot = FALSE");
    }
}


void GnuplotCanvas::setTerminal(Gnuplot& gp, /*Image2D& image, */_Format format, String nameFile)
{
    if (format == _gif)
        {
            gp << "set term " << _format[format] << " animate\n";
        }
        else
        {
            gp << "set term " << _format[format] << " \n";

        }   

    if (noCanvas == true){ // If no dimensions

    }
    else    
    {
        gp << "size " << to_string(sizeX) << " " <<
            to_string(sizeY) << "\n";

    }
    // create file
    gp << "set output \"" << nameFile << "." << _format[format] << "\"\n";
    

    //_filename = nameFile + "." + _format[format];
    //setColour(image);
    //setLineStyle(image);
}

void GnuplotCanvas::setXLabel(Gnuplot& gp, String xLabel)
{
    gp << "set xlabel '" << xLabel << "'\n";
}

void GnuplotCanvas::setYLabel(Gnuplot& gp, String yLabel)
{
    gp << "set ylabel '" << yLabel << "'\n";
}

void GnuplotCanvas::setZLabel(Gnuplot& gp, String zLabel)
{
    gp << "set zlabel '" << zLabel << "'\n";
}

void GnuplotCanvas::setTitle(Gnuplot& gp, String title)
{
    gp << "set title '" << title << "'\n";
}

void GnuplotCanvas::setXYZLabel(Gnuplot& gp, String xLabel, String yLabel, String zLabel = "")
{
    gp << "set xlabel '" << xLabel << "'\n";
    gp << "set ylabel '" << yLabel << "'\n";
    
    if (zLabel != "")
    {
        gp << "set zlabel '" << zLabel << "'\n";
    }

}

void GnuplotCanvas::setLabels(Gnuplot& gp, String xLabel, String yLabel, String zLabel, String title)
{
    gp << "set xlabel '" << xLabel << "'\n";
    gp << "set ylabel '" << yLabel << "'\n";
    gp << "set zlabel '" << zLabel << "'\n";
    gp << "set title '" << title << "'\n";
}

void GnuplotCanvas::setRange2D(_Range2D& range2D, FloatDP maxX, FloatDP maxY)
{
    range2D.Xmin = 0;
    range2D.Xmax = maxX;
    range2D.Ymin = 0;
    range2D.Ymax = maxY;
}

void GnuplotCanvas::setRange2D(_Range2D& range2D, FloatDP minX, FloatDP maxX, 
                FloatDP minY, FloatDP maxY)
{
    range2D.Xmin = minX;
    range2D.Xmax = maxX;
    range2D.Ymin = minY;
    range2D.Ymax = maxY;
}

void GnuplotCanvas::setRange3D(_Range3D& range3D, FloatDP minX, FloatDP maxX, 
                FloatDP minY, FloatDP maxY,
                FloatDP minZ, FloatDP maxZ)
{
    range3D.Xmin = minX;
    range3D.Xmax = maxX;
    range3D.Ymin = minY;
    range3D.Ymax = maxY;
    range3D.Zmin = minZ;
    range3D.Zmax = maxZ;
}

void GnuplotCanvas::setRange3D(_Range3D& range3D, FloatDP maxX, FloatDP maxY, FloatDP maxZ)
{
    range3D.Xmin = 0;
    range3D.Xmax = maxX;
    range3D.Ymin = 0;
    range3D.Ymax = maxY;
    range3D.Zmin = -maxZ;
    range3D.Zmax = maxZ;
}

void GnuplotCanvas::setLineStyle(Image2D& image)
{
    image.linestyle2D.style = lines2D;
    image.linestyle2D.ls = 1;
    image.linestyle2D.lw = 1;
}

void GnuplotCanvas::setLineStyle(Image3D& image)
{
    image.linestyle3D.style = lines3D;
    image.linestyle3D.lw = 1;
    image.linestyle3D.ls = 1;
}

void GnuplotCanvas::setLineStyle(Image2D& image, _Line2D line)
{
    image.linestyle2D.style = line.style;
    image.linestyle2D.ls = 1;
    image.linestyle2D.lw = 1;
}

void GnuplotCanvas::setLineStyle(Image2D& image, _Line2D line, int ls, int lw)
{
    image.linestyle2D.style = line.style;
    if (lw < 0)
    {
        lw = 0;
    }
    if (lw > 7)
    {
        lw = 7;
    }
    image.linestyle2D.lw = lw;
    if (ls < 0)
    {
        ls = 0;
    }
    if (ls > 7)
    {
        ls = 7;
    }
    image.linestyle2D.ls = ls;
}

void GnuplotCanvas::setLineStyle(Image3D& image, _Line3D line, int ls, int lw)
{
    image.linestyle3D.style = line.style;
    if (lw < 0)
    {
        lw = 0;
    }
    if (lw > 7)
    {
        lw = 7;
    }
    image.linestyle3D.lw = lw;
    if (ls < 0)
    {
        ls = 0;
    }
    if (ls > 7)
    {
        ls = 7;
    }
    image.linestyle3D.ls = ls;
}
//surface
void GnuplotCanvas::setLineStyle(Image3D& image, _Line3D line)
{
    image.linestyle3D.style = line.style;
    image.linestyle3D.lw = 1;
    image.linestyle3D.ls = 1;
}
// Set default colour
void GnuplotCanvas::setColour(Image2D& image)
{
    image.colour = _black;
}
void GnuplotCanvas::setColour(Image3D& image)
{
    image.colour = _black;
}
// Set colour
void GnuplotCanvas::setColour(Image2D& image, _Colours color)
{
    image.colour = color;
}

// Set Colour
void GnuplotCanvas::setColour(Image3D& image, _Colours color)
{
    image.colour = color;
}
void GnuplotCanvas::setXLogAxis(Gnuplot& gp)
{
    gp << "set logscale x\n";
}

void GnuplotCanvas::setYLogAxis(Gnuplot& gp)
{
    gp << "set logscale y\n";
}

void GnuplotCanvas::setXYLogAxis(Gnuplot& gp)
{
    gp << "set logscale xy\n";
}

void GnuplotCanvas::setXZLogAxis(Gnuplot& gp)
{
    gp << "set logscale xz\n";
}

void GnuplotCanvas::setYZLogAxis(Gnuplot& gp)
{
    gp << "set logscale yz\n";
}

void GnuplotCanvas::setXYZLogAxis(Gnuplot& gp)
{
    gp << "set logscale xyz\n";
}

void GnuplotCanvas::setLegend(Gnuplot& gp)
{
    gp << "set key default\n";
}
void GnuplotCanvas::setMap(Gnuplot& gp)
{
    gp << "set pm3d map\n";
}
void GnuplotCanvas::set3DPalette(Gnuplot& gp, Image3D& image, FloatDP min, FloatDP max, FloatDP step, bool s)
{
    if (s)
    {
        is3DPalette = true;
        gp << "set cbrange [" << to_string(min) << ":" << to_string(max) << "]\n";
        gp << "set cbtics " << to_string(step) << "\n";
        gp << "set palette defined\n";
        image.linestyle3D.style = pm3d;
    }  
}
void GnuplotCanvas::unsetColorbox(Gnuplot& gp)
{
    gp << "unset colorbox\n";
}

#endif // DEBUG

} // namespace Ariadne


