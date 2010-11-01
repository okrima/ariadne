/***************************************************************************
 *            procedure.h
 *
 *  Copyright 2010  Pieter Collins
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


/*! \file procedure.h
 *  \brief Procedure to compute a real function
 */

#ifndef ARIADNE_PROCEDURE_H
#define ARIADNE_PROCEDURE_H

#include <iostream>

#include "container.h"
#include "vector.h"

#include "operators.h"
#include "formula.h"
#include "expansion.h"

namespace Ariadne {

typedef uint Nat;
template<class X> class Formula;

struct ProcedureInstruction
{
    explicit ProcedureInstruction(Operator o, uint a) : op(o), arg(a) { }
    explicit ProcedureInstruction(Operator o, uint a1, uint a2) : op(o), arg1(a1), arg2(a2) { }
    explicit ProcedureInstruction(Operator o, uint a, int n) : op(o), arg(a), np(n) { }
    Operator op;
    union {
        struct { uint arg; int np; };
        struct { uint arg1; uint arg2; };
    };
};

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class X>
class Procedure {
  public:
    explicit Procedure();
    explicit Procedure(const Formula<X>& f);
    explicit Procedure(const Expansion<X>& f);
  public:
    List<X> _constants;
    List<ProcedureInstruction> _instructions;
  public:
    void new_instruction(Operator o, uint a) { _instructions.append(ProcedureInstruction(o,a)); }
    void new_instruction(Operator o, uint a, int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
    void new_instruction(Operator o, uint a1, uint a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
  public:
};

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class X>
class Vector< Procedure<X> > {
  public:
    explicit Vector(const Vector< Formula<X> >& f);
    explicit Vector(const Procedure<X>& p);
  public:
    List<X> _constants;
    List<ProcedureInstruction> _instructions;
    Vector<Nat> _results;
  public:
    Nat result_size() const { return _results.size(); }
    void new_instruction(Operator o, Nat a) { _instructions.append(ProcedureInstruction(o,a)); }
    void new_instruction(Operator o, Nat a, int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
    void new_instruction(Operator o, Nat a1, Nat a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
    void set_return(Nat i, Nat a) { _results[i]=a; }
};

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> void _execute(List<T>& v, const List<ProcedureInstruction>& p, const List<X>& c, const Vector<T>& x)
{
    T z=x[0]*0;
    for(uint i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        switch(instruction.op) {
            case CNST: v.append(z+c[instruction.arg]); break;
            case IND:  v.append(x[instruction.arg]); break;
            case ADD:  v.append(v[instruction.arg1]+v[instruction.arg2]); break;
            case SUB:  v.append(v[instruction.arg1]-v[instruction.arg2]); break;
            case MUL:  v.append(v[instruction.arg1]*v[instruction.arg2]); break;
            case DIV:  v.append(v[instruction.arg1]/v[instruction.arg2]); break;
            case POW:  v.append(pow(v[instruction.arg],instruction.np)); break;
            case ABS:  v.append(abs(v[instruction.arg])); break;
            case POS:  v.append(pos(v[instruction.arg])); break;
            case NEG:  v.append(neg(v[instruction.arg])); break;
            case REC:  v.append(rec(v[instruction.arg])); break;
            case SQR:  v.append(sqr(v[instruction.arg])); break;
            case SQRT: v.append(sqrt(v[instruction.arg])); break;
            case EXP:  v.append(exp(v[instruction.arg])); break;
            case LOG:  v.append(log(v[instruction.arg])); break;
            case SIN:  v.append(sin(v[instruction.arg])); break;
            case COS:  v.append(cos(v[instruction.arg])); break;
            case TAN:  v.append(tan(v[instruction.arg])); break;
            case ATAN:  v.append(atan(v[instruction.arg])); break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
    }
}

template<class T> void _propagate(Vector<T>& x, List<T>& v, const List<ProcedureInstruction>& p)
{
    static const Float inf=Ariadne::inf<Float>();
    ARIADNE_ASSERT(v.size()==p.size());
    for(uint r=p.size()-1u; r!=-1u; --r) {
        uint a=p[r].arg; uint a1=p[r].arg1; uint a2=p[r].arg2;
        switch(p[r].op) {
            case CNST: break;
            case IND: restrict(x[a],v[r]); break;
            case ADD: restrict(v[a1],v[r]-v[a2]); restrict(v[a2],v[r]-v[a1]); break;
            case SUB: restrict(v[a1],v[a2]+v[r]); restrict(v[a2],v[a1]-v[r]); break;
            case MUL: restrict(v[a1],v[r]/v[a2]); restrict(v[a2],v[r]/v[a1]); break;
            case DIV: restrict(v[a1],v[a2]*v[r]); restrict(v[a2],v[a1]/v[r]); break;
            case MAX: restrict(v[a1],max(v[r],v[a2])); restrict(v[a2],max(v[r],v[a1])); break;
            case POS: restrict(v[a],v[r]); break;
            case NEG: restrict(v[a],neg(v[r])); break;
            case REC: restrict(v[a],rec(v[r])); break;
            case SQR: restrict(v[a],sqrt(v[r])); break;
            case POW: restrict(v[a],exp(log(v[r])/p[r].np)); break;
            case SQRT: restrict(v[a],sqr(v[r])); break;
            case EXP: restrict(v[a],log(v[r])); break;
            case LOG: restrict(v[a],exp(v[r])); break;
            case SIN: restrict(v[a],asin(v[r])); break;
            case COS: restrict(v[a],acos(v[r])); break;
            case TAN: restrict(v[a],atan(v[r])); break;
            case ATAN: restrict(v[a],tan(v[r])); break;
            case EQ: restrict(v[a1],v[r]); restrict(v[a2],v[r]); break;
            case LEQ: restrict(v[a1],Interval(-inf,v[a2].upper())); restrict(v[a1],Interval(v[a2].lower(),+inf)); break;
            default: ARIADNE_THROW(std::runtime_error,"_propagate(Vector<T>,List<T>,List<ProcedureInstruction>)","Unhandled operator "<<p[r].op);
        }
    }
}

template<class X> uint _convert(List<ProcedureInstruction>& p, List<X>& c, const FormulaNode<X>* f, Map<const FormulaNode<X>*,uint>& ind) {
    typedef ProcedureInstruction PI;
    if(ind.has_key(f)) { return ind[f]; }
    switch(f->op) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case CNST: p.append(PI(CNST,c.size())); c.append(*f->val); break;
        case IND: p.append(PI(IND,f->ind)); break;
        case ADD: case SUB: case MUL: case DIV:
            p.append(PI(f->op,_convert(p,c,f->arg1,ind),_convert(p,c,f->arg2,ind))); break;
        case POW:
            p.append(PI(f->op,_convert(p,c,f->arg,ind),f->np)); break;
        case NEG: case REC: case SQR: case SQRT:
        case EXP: case LOG: case SIN: case COS: case TAN: case ATAN:
            p.append(PI(f->op,_convert(p,c,f->arg,ind))); break;
        default: assert(false);
    }
    const uint num=p.size()-1;
    ind.insert(f,num);
    return num;
}

template<class X>
void _write(std::ostream& os, const List<ProcedureInstruction>& p, const List<X>& c)
{
    for(uint i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        os << "v[" << i << "]=";
        switch(instruction.op) {
            case CNST:
                os << c[instruction.arg]; break;
            case IND:
                os << "x[" << instruction.arg << "]"; break;
            case ADD: case SUB: case MUL: case DIV:
                os << "v[" << instruction.arg1 << "]" << symbol(instruction.op) << "v[" << instruction.arg2 << "]"; break;
            case POW:
                os<<"pow(v[" << instruction.arg << "],"<<instruction.np<<")"; break;
            case POS: case NEG:
                os << symbol(instruction.op) << "v[" << instruction.arg << "]"; break;
            case ABS: case REC: case SQR: case SQRT:
            case EXP: case LOG: case SIN: case COS: case TAN:
            case ASIN: case ACOS: case ATAN:
                os << instruction.op << "(v[" << instruction.arg << "])"; break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
        os << "; ";
    }
}


template<class X>
Procedure<X>::Procedure()
{
}

template<class X>
Procedure<X>::Procedure(const Formula<X>& f)
{
    Map<const FormulaNode<X>*, uint> ind;
    _convert(this->_instructions,this->_constants,f._root.operator->(), ind);
}

template<class X>
Procedure<X>::Procedure(const Expansion<X>& e)
{
    Vector< Formula<X> > id=Formula<X>::identity(e.argument_size());
    Formula<X> f=horner_evaluate(e,id);
    Map<const FormulaNode<X>*, uint> ind;
    _convert(this->_instructions,this->_constants,f._root.operator->(), ind);
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> inline T evaluate(const Procedure<X>& f, const Vector<T>& x)
{
    List<T> e;
    _execute(e,f._instructions,f._constants,x);
    return e.back();
}

template<class X> Procedure<X>& operator+=(Procedure<X>& f, const X& c) {
    f._constants.append(c);
    f.new_instruction(CNST,f._constants.size()-1);
    f.new_instruction(ADD,f._instructions.size()-1,f._instructions.size()-2);
    return f;
}

template<class X>
std::ostream& operator<<(std::ostream& os, const Procedure<X> f) {
    os<<"Procedure( ";
    _write(os,f._instructions,f._constants);
    os << "r=v[" << f._instructions.size()-1u << "] )";
    return os;
}


template<class X>
Vector< Procedure<X> >::Vector(const Vector< Formula<X> >& f)
    : _results(f.size(),0u)
{
    Map<const FormulaNode<X>*, uint> ind;
    for(uint i=0; i!=f.size(); ++i) {
        _convert(this->_instructions,this->_constants,f[i]._root.operator->(), ind);
    }
    for(Nat i=0; i!=f.size(); ++i) {
        this->_results[i]=ind[f[i]._root.operator->()];
    }
}

template<class X>
Vector< Procedure<X> >::Vector(const Procedure<X>& f)
    : _instructions(f._instructions)
    , _constants(f._constants)
    , _results(1u,f._instructions.size()-1u)
{
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> inline Vector<T> evaluate(const Vector< Procedure<X> >& f, const Vector<T>& x)
{
    List<T> v;
    _execute(v,f._instructions,f._constants,x);
    Vector<T> r(f.result_size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=v[f._results[i]];
    }
    return r;
}

template<class X>
std::ostream& operator<<(std::ostream& os, const Vector< Procedure<X> > f) {
    os<<"Procedure( ";
    _write(os,f._instructions,f._constants);
    os << "r=v" << f._results << " )";
    return os;
}


} // namespace Ariadne


#include "numeric.h"

namespace Ariadne {

inline
void restrict(Interval& r, const Interval& x) {
    r.set_lower(max(r.lower(),x.lower()));
    r.set_upper(min(r.upper(),x.upper()));
};

template<class X>
void simple_hull_reduce(Vector<Interval>& dom, const Procedure<X>& f, Interval codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<X>& c=f._constants;
    List<Interval> v; v.reserve(p.size());

    _execute(v,p,c,dom);
    restrict(v.back(),codom);
    _propagate(dom,v,p);
}

template<class X>
void simple_hull_reduce(Vector<Interval>& dom, const Vector< Procedure<X> >& f, Vector<Interval> codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<X>& c=f._constants;
    List<Interval> v; v.reserve(p.size());

    ARIADNE_ASSERT(codom.size()==f._results.size());

    _execute(v,p,c,dom);
    for(uint i=0; i!=codom.size(); ++i) {
        restrict(v[f._results[i]],codom[i]);
    }
    _propagate(dom,v,p);
}

} // namespace Ariadne


#endif // ARIADNE_PROCEDURE_H
