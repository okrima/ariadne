/***************************************************************************
 *            solver.cc
 *
 *  Copyright  2006-9  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "functional.h"
#include "config.h"

#include "interval.h"
#include "function_model.h"

#include "solver.h"

#include "logging.h"
#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "taylor_model.h"
#include "function.h"
#include "function_model.h"

namespace Ariadne {

namespace {

Vector<Interval> ranges(const Vector<ValidatedTaylorModel>& f) {
    Vector<Interval> r(f.size()); for(uint i=0; i!=f.size(); ++i) { r[i]=f[i].range(); } return r;
}

Vector<ValidatedTaylorModel>& clobber(Vector<ValidatedTaylorModel>& h) {
    for(uint i=0; i!=h.size(); ++i) { h[i].set_error(0.0); } return h; }

// Compute the Jacobian over an arbitrary domain
Matrix<ValidatedNumberType>
jacobian2(const Vector<ValidatedTaylorModel>& f, const Vector<ValidatedNumberType>& x)
{
    Vector< Differential<ValidatedNumberType> > dx(x.size());
    for(uint i=0; i!=x.size()-f.size(); ++i) {
        dx[i]=Differential<ValidatedNumberType>::constant(f.size(),1u,x[i]); }
    for(uint i=0; i!=f.size(); ++i) {
        uint j=i+(x.size()-f.size());
        dx[j]=Differential<ValidatedNumberType>::variable(f.size(),1u,x[j],i); }
    Vector< Differential<ValidatedNumberType> > df(f.size());
    for(uint i=0; i!=f.size(); ++i) {
        df[i]=evaluate(f[i].expansion(),dx);
    }
    Matrix<ValidatedNumberType> J=jacobian(df);
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<ExactFloatType>
jacobian2_value(const Vector<ValidatedTaylorModel>& f)
{
    const uint rs=f.size();
    const uint fas=f.zero_element().argument_size();
    const uint has=fas-rs;
    Matrix<ExactFloatType> J(rs,rs);
    MultiIndex a(fas);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=rs; ++j) {
            a[has+j]=1; const ExactFloatType x=f[i][a]; J[i][j]=x; a[has+j]=0;
        }
    }
    return J;
}

// Compute the Jacobian over the unit domain
Matrix<ValidatedNumberType>
jacobian2_range(const Vector<ValidatedTaylorModel>& f)
{
    uint rs=f.size();
    uint fas=f.zero_element().argument_size();
    uint has=fas-rs;
    Matrix<ValidatedNumberType> J(rs,rs);
    for(uint i=0; i!=rs; ++i) {
        for(ValidatedTaylorModel::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
            for(uint k=0; k!=rs; ++k) {
                const uint c=iter->key()[has+k];
                if(c>0) {
                    const ExactFloatType& x=iter->data();
                    if(iter->key().degree()==1) { J[i][k]+=x; }
                    else { J[i][k]+=ValidatedNumberType(-1,1)*x*c; }
                    //std::cerr<<"  J="<<J<<" i="<<i<<" a="<<iter->key()<<" k="<<k<<" c="<<c<<" x="<<x<<std::endl;
                }
            }
        }
    }
    return J;
}



} // namespace


static const bool ALLOW_PARTIAL_FUNCTION = true;

FunctionModelFactoryInterface<ValidatedTag>* make_taylor_function_factory();

ValidatedVectorFunctionModel operator*(const Matrix<ExactFloatType>& A,const ValidatedVectorFunctionModel& v) {
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorFunctionModel r(A.row_size(),v[0].create_zero());
    for(uint i=0; i!=r.size(); ++i) {
        ValidatedScalarFunctionModel t=r[i];
        for(uint j=0; j!=v.size(); ++j) {
            t+=ExactFloat(A[i][j])*v[j];
        }
        r[i]=t;
    }
    return r;
}

ValidatedVectorFunctionModel operator*(const Matrix<ValidatedNumberType>& A,const ValidatedVectorFunctionModel& v) {
    ARIADNE_ASSERT(v.size()!=0);
    ValidatedVectorFunctionModel r(A.row_size(),v[0].create_zero());
    for(uint i=0; i!=r.size(); ++i) {
        ValidatedScalarFunctionModel t=r[i];
        for(uint j=0; j!=v.size(); ++j) {
            t+=A[i][j]*v[j];
        }
        r[i]=t;
    }
    return r;
}

ErrorFloatType sup_error(const ValidatedVectorFunctionModel& x) {
    ErrorFloatType r=0.0;
    for(uint i=0; i!=x.size(); ++i) { r=std::max(r,x[i].error()); }
    return r;
}







SolverBase::SolverBase(double max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps), _function_factory_ptr(make_taylor_function_factory())
{
}

void
SolverBase::set_function_factory(const FunctionModelFactoryInterface<ValidatedTag>& factory)
{
    this->_function_factory_ptr=std::shared_ptr< FunctionModelFactoryInterface<ValidatedTag> >(factory.clone());
}

const FunctionModelFactoryInterface<ValidatedTag>&
SolverBase::function_factory() const
{
    return *this->_function_factory_ptr;
}


Vector<ValidatedNumberType>
SolverBase::solve(const ValidatedVectorFunction& f,
                  const Vector<ValidatedNumberType>& ix) const
{
    return this->solve(f,static_cast<Vector<Interval>>(ix));
}

Vector<ValidatedNumberType>
SolverBase::solve(const ValidatedVectorFunction& f,
                  const Box& bx) const
{
    Set< Vector<ValidatedNumberType> > r = this->solve_all(f,bx);
    if(r.size()==0u) { ARIADNE_THROW(NoSolutionException,"SolverBase::solve","no solution in solve("<<f<<","<<bx<<")"); }
    if(r.size()!=1u) { ARIADNE_THROW(SolverException,"SolverBase::solve","non-unique solution in solve("<<f<<","<<bx<<")"); }
    return *r.begin();
}

void
solve_all(Set< Vector<ValidatedNumberType> >& r,
          const SolverInterface& s,
          const ValidatedVectorFunction& f,
          const Box& ix);

Set< Vector<ValidatedNumberType> >
SolverBase::solve_all(const ValidatedVectorFunction& f,
                      const Box& ix) const
{
    ARIADNE_LOG(5,"SolverBase::solve_all(f,ix): f="<<f<<", ix="<<ix<<"\n");

    // Create result set
    Set< Vector<ValidatedNumberType> > r;

    Vector<ValidatedFloatType> x=make_singleton(ix);

    // Test for no solution
    const Vector<ValidatedNumberType> z(ix.size());
    if(inconsistent(f.evaluate(x),z)) {
        return r;
    }

    bool invertible_jacobian=true;
    //Vector<ValidatedNumberType> nx=2*ix-make_exact(ix);
    Vector<ValidatedNumberType> nx=x;
    try {
        Matrix<ValidatedNumberType> Jinv=inverse(f.jacobian(nx));
    }
    catch(const SingularMatrixException& e) {
        invertible_jacobian=false;
    }

    bool need_to_split=true;

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            Vector<ValidatedNumberType> y=this->solve(f,nx);
            bool is_new=true;
            for(Set<Vector<ValidatedNumberType> >::const_iterator iter=r.begin(); iter!=r.end(); ++iter) {
                if(consistent(y,*iter)) {
                    is_new=false;
                    break;
                }
            }
            if(is_new) {
                r.insert(y);
            }
            need_to_split=false;
        }
        catch(const NoSolutionException& e) {
            //ARIADNE_WARN("NoSolutionException exception: "<<e.what()<<"\n");
            need_to_split=false;
        }
        catch(const SolverException& e) {
            ARIADNE_WARN("SolverException exception: "<<e.what()<<"\n");
            need_to_split=true;
            // No solution found, try splitting
        }
    } else {
        need_to_split=true;
    }

    if(need_to_split) {
        // If sup_error is too small, assume solution is not verified
        if(sup_error(make_singleton(ix))<this->maximum_error()) {
            if(!invertible_jacobian) {
                ARIADNE_WARN("Cannot verify solution in "<<ix<<" with f="<<f(make_singleton(ix))<<"); "
                             <<"Jacobian "<<f.jacobian(nx)<<" is not invertible; "
                             <<"approximate inverse="<<inverse(midpoint(f.jacobian(nx)))<<"\n");
            } else {
                ARIADNE_WARN("Cannot verify or falsify solution in "<<ix<<"; f("<<ix<<")="<<f(make_singleton(ix))<<".\n");
            }
            return r;
        }

        //std::cerr<<"  Splitting "<<ix<<"\n";
        std::pair< Vector<Interval>, Vector<Interval> > splt=split(ix);
        r.adjoin(this->solve_all(f,splt.first));
        r.adjoin(this->solve_all(f,splt.second));
    }

    return r;
}





Vector<ValidatedNumberType>
SolverBase::zero(const ValidatedVectorFunction& f,
                 const Box& bx) const
{
    const double& e=this->maximum_error();
    uint n=this->maximum_number_of_steps();
    ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
    Vector<ValidatedNumberType> r=make_singleton(bx);
    Vector<ValidatedNumberType> nr(r.size());
    bool has_solution=false;
    while(n>0) {
        nr=this->step(f,r);
        ARIADNE_LOG(5,"  nr="<<nr<<"\n");

        if(!has_solution && refines(nr,r)) {
            has_solution=true;
        }

        if(has_solution && sup_error(nr) < e) {
            return nr;
        }

        if(!consistent(nr,r)) {
            ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; "<<nr<<" is inconsistent with "<<r);
        }
        r=refinement(nr,r);
        n=n-1;
    }
    if(!consistent(f.evaluate(r),Vector<ValidatedNumberType>(f.result_size()))) {
        ARIADNE_THROW(NoSolutionException,"SolverBase::zero","No result found in "<<bx<<"; f("<<r<<") is inconsistent with zero");
    } else {
        UpperFloatType widen=ExactFloatType(eps())*sup_error(r);
        r+=Vector<ValidatedNumberType>(r.size(),ValidatedNumberType(-widen,+widen));
        nr=this->step(f,r);
        if(refines(nr,r)) {
            return nr;
        }
        ARIADNE_THROW(SolverException,"SolverBase::zero","No result verified in "<<bx<<"; maximum number of steps reached with approximation "<<r<<" which cannot be robustly checked");
    }
}



Vector<ValidatedNumberType>
SolverBase::fixed_point(const ValidatedVectorFunction& f, const Box& bx) const
{
    ValidatedVectorFunction id=ValidatedVectorFunction::identity(f.argument_size());
    return this->solve(f-id,bx);
}


ValidatedVectorFunctionModel
SolverBase::implicit(const ValidatedVectorFunction& f,
                      const Box& ip,
                      const Box& ix) const
{
    ARIADNE_LOG(4,"SolverBase::implicit(ValidatedVectorFunction f, IntervalVector ip, IntervalVector ix)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_ASSERT(f.result_size()==ix.size());
    ARIADNE_ASSERT(f.argument_size()==ip.size()+ix.size());

    const uint n=ix.size();
    const double err=this->maximum_error();

    ValidatedVectorFunctionModel id(this->function_factory().create_identity(ip));
    ValidatedVectorFunctionModel h(this->function_factory().create_constants(ip,make_singleton(ix)));
    ValidatedVectorFunctionModel nh(this->function_factory().create_zeros(n,ip));
    ValidatedVectorFunctionModel fnh(this->function_factory().create_zeros(n,ip));

    uint steps_remaining=this->maximum_number_of_steps();
    uint number_unrefined=n;
    Array<bool> refinement(n,false);

    while(steps_remaining>0) {
        ARIADNE_LOG(5,"\n");
        ARIADNE_LOG(5,"step="<<this->maximum_number_of_steps()-steps_remaining<<"\n");
        nh=this->implicit_step(f,id,h);
        fnh=compose(f,join(id,nh));
        ARIADNE_LOG(5,"  nh="<<nh<<"\n");
        ARIADNE_LOG(5,"  fnh="<<fnh<<"\n");

        if(ALLOW_PARTIAL_FUNCTION) {
            for(uint i=0; i!=n; ++i) {
                if(!refinement[i]) {
                    ARIADNE_LOG(6,"refines(nh["<<i<<"],x["<<i<<"])="<<refines(nh[i],h[i])<<"\n");
                    if(refines(nh[i],h[i])) {
                        refinement[i]=true;
                        --number_unrefined;
                    } else if(disjoint(ValidatedScalarFunctionModel(nh[i]).range(),ix[i])) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; function "<<nh<<" is disjoint from "<<h<<" for at least one point.");
                    }
                }
            }

            ARIADNE_LOG(6,"refinement="<<refinement<<"\n");
            ARIADNE_LOG(6,"nh.range()="<<nh.range()<<"\n");
            ARIADNE_LOG(6,"sup_error(nh)="<<sup_error(nh)<<"\n");
            ARIADNE_LOG(6,"sup_error(fnh)="<<sup_error(fnh)<<"\n");

            h=nh;

        } else { // !ALLOW_PARTIAL_FUNCTION
            for(uint i=0; i!=n; ++i) {
                if(!refinement[i]) {
                    if(refines(nh[i],h[i])) {
                        refinement[i]=true;
                        --number_unrefined;
                    } else if(disjoint(nh[i],h[i])) {
                        ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","No result found in "<<ix<<"; "<<nh<<" is disjoint from "<<h);
                    }
                }
            }

            h=intersection(nh,h);
        }

        steps_remaining=steps_remaining-1;
        if( (number_unrefined==0) && ( (steps_remaining==0) || ((sup_error(nh)<err) && (sup_error(fnh)<err)) ) ) {
            return h;
        }
    }

    ARIADNE_THROW(NoSolutionException,"SolverBase::implicit","Could not prove existence of a solution in "<<ix<<".");


}

ValidatedScalarFunctionModel
SolverBase::implicit(const ValidatedScalarFunction& f,
                     const Box& ip,
                     const Interval& ix) const
{
    ARIADNE_LOG(4,"SolverBase::implicit(ValidatedScalarFunction f, IntervalVector ip, Interval ix)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ValidatedVectorFunctionModel res=this->implicit(ValidatedVectorFunction(List<ValidatedScalarFunction>(1u,f)),ip,Box(1u,ix));
    return res[0];
}

ValidatedVectorFunctionModel
SolverBase::continuation(const ValidatedVectorFunction& f,
                         const Vector<ApproximateNumberType>& p,
                         const Box& ix,
                         const Box& ip) const
{
    ARIADNE_NOT_IMPLEMENTED;
}




Vector<ValidatedNumberType>
IntervalNewtonSolver::step(const ValidatedVectorFunction& f,
                           const Vector<ValidatedNumberType>& x) const
{
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumberType> im(m);
    Vector<ValidatedNumberType> w=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<ValidatedNumberType> A=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<ValidatedNumberType> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<ValidatedNumberType> dx=Ainv*w;
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumberType> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    return nx;
}

Vector<ValidatedNumberType>
KrawczykSolver::step(const ValidatedVectorFunction& f,
                     const Vector<ValidatedNumberType>& x) const
{
    Matrix<ValidatedNumberType> I=Matrix<ValidatedNumberType>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumberType> im(m);
    Vector<ValidatedNumberType> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<ValidatedNumberType> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<ValidatedNumberType> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<ValidatedNumberType> dx=M*fm-(I-M*J)*(x-m);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumberType> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<ValidatedNumberType> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


Vector<ValidatedNumberType>
FactoredKrawczykSolver::step(const ValidatedVectorFunction& f,
                             const Vector<ValidatedNumberType>& x) const
{
    Matrix<ValidatedNumberType> I=Matrix<ValidatedNumberType>::identity(x.size());
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    Vector<ExactFloat> m(make_exact(x));
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<ValidatedNumberType> im(m);
    Vector<ValidatedNumberType> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<ValidatedNumberType> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<ValidatedNumberType> mJ(midpoint(J));
    Matrix<ValidatedNumberType> M=inverse(mJ);
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<ValidatedNumberType> dx=M*(fm+(J-mJ)*(x-m));
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<ValidatedNumberType> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<ValidatedNumberType> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");
    return nr;
}


ValidatedVectorFunctionModel
IntervalNewtonSolver::implicit_step(const ValidatedVectorFunction& f,
                                    const ValidatedVectorFunctionModel& id,
                                    const ValidatedVectorFunctionModel& h) const
{
    const uint m=id.size();
    const uint n=h.size();

    ARIADNE_LOG(6,"IntervalNewtonSolver::implicit_step(ValidatedVectorFunction f, ValidatedVectorFunctionModel id, ValidatedVectorFunctionModel h)\n");
    ARIADNE_LOG(7,"f="<<f<<"\n");
    ARIADNE_LOG(7,"h="<<h<<"\n");

    ValidatedVectorFunctionModel mh=h; mh.clobber();
    ARIADNE_LOG(7,"midpoint(h)="<<mh<<"\n");

    Matrix<ValidatedScalarFunction> D2f(n,n);
    for(uint i=0; i!=n; ++i) {
        for(uint j=0; j!=n; ++j) {
            D2f[i][j]=f[i].derivative(m+j);
        }
    }
    ARIADNE_LOG(7,"D2f="<<D2f<<"\n");

    ValidatedScalarFunctionModel z=h[0]*Real(0);
    ValidatedVectorFunctionModel idh=join(id,h);

    Matrix<ValidatedScalarFunctionModel> J(n,n,z);
    for(uint i=0; i!=n; ++i) {
        for(uint j=0; j!=n; ++j) {
            J[i][j]=compose(D2f[i][j],idh);
        }
    }
    ARIADNE_LOG(7,"J="<<J<<"\n");

    Matrix<Interval> rngJ(n,n);
    for(uint i=0; i!=n; ++i) {
        for(uint j=0; j!=n; ++j) {
            Interval D2fij=Interval(unchecked_evaluate(D2f[i][j],make_singleton(join(id.range(),h.range()))));
            rngJ[i][j]=intersection(J[i][j].range(),D2fij);
        }
    }
    ARIADNE_LOG(7,"rngJ="<<rngJ<<"\n");

    ValidatedVectorFunctionModel fidmh=compose(f,join(id,mh));
    ARIADNE_LOG(7,"compose(f,join(id,midpoint(h)))="<<fidmh<<"\n");

    ValidatedVectorFunctionModel dh(n,z);
    if(n==1) {
        if(contains(rngJ[0][0],ExactFloatType(0.0))) {
            ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","D2f(P,X)="<<rngJ[0][0]<<" which contains zero.");
        }
        if(contains(J[0][0].range(),ExactFloatType(0.0))) {
            dh[0]=fidmh[0]/make_singleton(rngJ[0][0]);
        } else {
            dh[0]=fidmh[0]/J[0][0];
        }
    } else {
        dh=inverse(make_singleton(rngJ))*fidmh;
    }
    ARIADNE_LOG(7,"dh="<<dh<<"\n");

    ValidatedVectorFunctionModel nh=mh-dh;
    ARIADNE_LOG(7,"nh="<<nh<<"\n");
    return nh;
}


ValidatedVectorFunctionModel
KrawczykSolver::implicit_step(const ValidatedVectorFunction& f,
                              const ValidatedVectorFunctionModel& p,
                              const ValidatedVectorFunctionModel& x) const
{
    const uint np=p.size();
    const uint nx=x.size();
    Matrix<ValidatedNumberType> I=Matrix<ValidatedNumberType>::identity(nx);
    ARIADNE_LOG(4,"  Contracting x="<<x<<"\n");
    ARIADNE_LOG(4,"    p="<<p<<"\n");
    ARIADNE_LOG(4,"    f="<<f<<"\n");
    //ARIADNE_LOG(5,"  e="<<sup_error(x)<<"  x="<<x<<"\n");
    ValidatedVectorFunctionModel mx(x);
    for(uint i=0; i!=mx.size(); ++i) { mx[i].set_error(0.0); }
    ARIADNE_LOG(5,"    mx="<<mx<<"\n");
    Vector<ErrorFloatType> ex(nx);
    for(uint i=0; i!=nx; ++i) { ex[i]=x[i].error(); }
    Vector<ValidatedNumberType> eix=make_bounds(ex);
    ARIADNE_LOG(5,"    ex="<<ex<<"\n");
    ValidatedVectorFunctionModel fm=compose(f,join(p,mx));
    ARIADNE_LOG(5,"    f(p,mx)="<<fm<<"\n");
    Vector<ValidatedNumberType> rp(np);
    for(uint i=0; i!=np; ++i) { rp[i]=make_singleton(p[i].range()); }
    Vector<ValidatedNumberType> rx(nx);
    for(uint i=0; i!=nx; ++i) { rx[i]=make_singleton(x[i].range()); }
    Matrix<ValidatedNumberType> J=project(f.jacobian(join(rp,rx)),range(0,nx),range(np,np+nx));
    ARIADNE_LOG(5,"    D2f(r)=J="<<J<<"\n");
    Matrix<ValidatedNumberType> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"    inverse(D2f(m))=M="<<M<<"\n");
    ARIADNE_LOG(5,"    M*f(p,mx)="<<M*fm<<"\n");
    ARIADNE_LOG(5,"    (I-M*J)="<<(I-M*J)<<"\n");
    ARIADNE_LOG(5,"    (I-M*J) * (ex*ValidatedNumberType(-1,+1))="<<(I-M*J)<<"*"<<eix<<"="<<(I-M*J) * eix<<"\n");
    ValidatedVectorFunctionModel dx= M*fm - (I-M*J) * eix;
    ARIADNE_LOG(5,"    dx="<<dx<<"\n");
    ValidatedVectorFunctionModel nwx= mx - dx;
    ARIADNE_LOG(5,"    nwx="<<nwx<<"\n");
    return nwx;
}


/*
ValidatedScalarFunctionModel
IntervalNewtonSolver::implicit(const ValidatedScalarFunction& f,
                               const Box& ip,
                               const Interval& ix) const
{
    return this->SolverBase::implicit(f,ip,ix);
    ARIADNE_LOG(4,"IntervalNewtonSolver::implicit(ValidatedScalarFunction f, IntervalVector P, Interval X)\n");
    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_LOG(5,"P="<<ip<<"\n");
    ARIADNE_LOG(5,"X="<<std::setprecision(17)<<ix<<"\n");
    ARIADNE_ASSERT_MSG(f.argument_size()==ip.size()+1u,"f="<<f<<", P="<<ip<<", X="<<ix<<"\n");
    const uint n=ip.size();
    ValidatedScalarFunction df=f.derivative(n);
    ARIADNE_LOG(5,"df="<<df<<"\n");

    // Simple check to see if there is definitely no solution
    Interval fipix=evaluate(f,join(ip,ix));
    ARIADNE_LOG(5,"f(P,X)="<<fipix<<"\n");
    ARIADNE_LOG(5,"f(p,X)="<<f(join(IntervalVector(midpoint(ip)),ix))<<"\n");
    if(fipix.lower()>0.0 || fipix.upper()<0.0) {
        ARIADNE_THROW(NoSolutionException,"IntervalNewtonSolver","f(P,X)="<<fipix<<" which does not contain zero.");
    }

    // Simple test of nonsingularity of Jacobian
    Interval dfipix=evaluate(df,join(ip,ix));
    ARIADNE_LOG(5,"df(P,X)="<<dfipix<<"\n");
    if(dfipix.lower()<=0.0 && dfipix.upper()>=0.0) {
        ARIADNE_THROW(SingularJacobianException,"IntervalNewtonSolver","df(P,X)="<<dfipix<<" which contains zero.");
    }


    // Set up auxiliary functions
    ValidatedScalarFunctionModel h=this->function_factory().create_constant(ip,ix);
    ValidatedVectorFunctionModel id=this->function_factory().create_identity(ip);
    ValidatedScalarFunctionModel dh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel mh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel nh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel dfidh=this->function_factory().create_zero(ip);
    ValidatedScalarFunctionModel fidmh=this->function_factory().create_zero(ip);

    ARIADNE_LOG(5,"f="<<f<<"\n");
    ARIADNE_LOG(5,"df="<<df<<"\n");
    ARIADNE_LOG(5,"P="<<ip<<"\n");
    ARIADNE_LOG(5,"X="<<ix<<"\n");
    ARIADNE_LOG(5,"h="<<h<<"\n");
    ARIADNE_LOG(5,"id="<<id<<"\n");

    // Compute solution
    bool refinement=false;
    for(uint i=0; i!=this->maximum_number_of_steps(); ++i) {
        mh=h; mh.clobber();
        ARIADNE_LOG(7,"mh="<<mh<<"\n");
        ValidatedNumberType dfiph=evaluate(df,join(ip,intersection(ix,h.range())));
        ARIADNE_LOG(7,"df(P,h(P))="<<dfiph<<"\n");
        dfidh=compose(df,join(id,h));
        ARIADNE_LOG(7,"df(id,h)="<<dfidh<<"\n");
        fidmh=compose(f,join(id,mh));
        ARIADNE_LOG(7,"f(p,mh(p))="<<fidmh<<"\n");
        // Prefer to divide by df(id,h), but may divide by constant df(P,h(P)) if latter is tighter
        if(refines(dfidh.range(),dfiph)) {
            dh=fidmh/dfidh;
        } else {
            dh=fidmh/dfiph;
        }
        ARIADNE_LOG(6,"dh="<<dh<<"\n");
        ARIADNE_LOG(6,"dh.range()="<<dh.range()<<"\n");
        ARIADNE_LOG(8,"norm(dh)="<<norm(dh)<<", h.error()="<<h.error()<<"\n");
        ErrorFloatType herr=h.error();
        ErrorFloatType dhnrm=norm(dh);
        // Prefer to check refinement using norm of dh, since this is faster
        // if(refines(nh,h)) { refinement=true; }
        if(dhnrm<=herr) { refinement=true; }
        nh=mh-dh;
        ARIADNE_LOG(6,"nh="<<nh<<"\n");
        ARIADNE_LOG(6,"nh.range()="<<nh.range()<<"\n");
        if(disjoint(nh.range(),ix)) {
            ARIADNE_THROW(NoSolutionException,"IntervalNewtonSolver",
                          "No solution of f(a,x)=0 with f="<<f<<", a in "<<ip<<", x in "<<ix<<".");
        }
        h=nh;
        if(refinement && dhnrm<this->maximum_error()/2) { break; }
    }
    return h;
    if(!refinement) {
        ARIADNE_LOG(6,"refinement="<<refinement<<"\n");
        ARIADNE_THROW(UnknownSolutionException,"IntervalNewtonSolver",
                      "Cannot find solution of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<".");
    }
    if(!refines(h.range(),ix)) {
        ARIADNE_LOG(6,"(h.range(),ix)="<<h.range()<<" "<<ix<<"\n");
        ARIADNE_THROW(UnknownSolutionException,"IntervalNewtonSolver",
                      "Range "<<h.range()<<" of solution h of f(x,h(x))=0 with f="<<f<<", x="<<ip<<", h in "<<ix<<" lies outside permitted interval.");
    }

}
*/


void IntervalNewtonSolver::write(std::ostream& os) const
{
    os << "IntervalNewtonSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

void KrawczykSolver::write(std::ostream& os) const
{
    os << "KrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}

void FactoredKrawczykSolver::write(std::ostream& os) const
{
    os << "FactoredKrawczykSolver"
       << "( maximum_error=" << this->maximum_error()
       << ", maximum_number_of_steps=" << this->maximum_number_of_steps()
       << " )";
}


} // namespace Ariadne
