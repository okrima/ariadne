/***************************************************************************
 *            variables.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file variables.h
 *  \brief Internal variables
 */

#ifndef ARIADNE_VARIABLES_H
#define ARIADNE_VARIABLES_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "macros.h"
#include "pointer.h"
#include "container.h"
#include "tribool.h"

namespace Ariadne {


template<class T> class Set;

typedef bool Boolean;
typedef tribool Tribool;
typedef std::string String;
class Integer;
class Real;

class Identifier : public String
{
  public:
    Identifier() : String() { }
    Identifier(const char* cstr) : String(cstr) { }
    Identifier(const String& str) : String(str) { }
};

class UntypedVariable;
template<class T> class ExtendedVariable;
template<class T> class Variable;
template<class T> class LetVariable;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;

template<class T> class Constant;
template<class T> class Expression;
template<class LHS,class RHS> class Assignment;

// Simplifying typedefs
typedef Variable<String> StringVariable;
typedef Constant<String> StringConstant;
typedef PrimedVariable<String> PrimedStringVariable;
typedef Variable<Integer> IntegerVariable;
typedef PrimedVariable<Integer> PrimedIntegerVariable;
typedef ExtendedVariable<Real> ExtendedRealVariable;
typedef LetVariable<Real> LetRealVariable;
typedef DottedVariable<Real> DottedRealVariable;
typedef PrimedVariable<Real> PrimedRealVariable;
typedef Variable<Real> RealVariable;
typedef Constant<Real> RealConstant;

//! \ingroup ExpressionModule
//! A named constant of type \a T.
template<class T> class Constant
{
  public:
    template<class X> explicit Constant(const String& str, const X& value)
        : _name_ptr(new Identifier(str)), _value_ptr(new T(value)) { }
    //explicit Constant(const String& str, const T& value)
    //    : _name_ptr(new String(str)), _value_ptr(new T(value)) { }
    const Identifier& name() const { return *_name_ptr; }
    const T& value() const { return *_value_ptr; }
    bool operator==(const Constant<T>& other) const {
        if(this->name()==other.name()) { assert(this->value()==other.value()); return true; } else { return false; } }
  private:
    shared_ptr<Identifier> _name_ptr;
    shared_ptr<T> _value_ptr;
};

template<> class Constant<String>
{
  public:
    explicit Constant(const String& value) : _value(value) { }
    const Identifier& name() const { return static_cast<const Identifier&>(_value); }
    const String& value() const { return _value; }
    operator const String& () const { return _value; }
    bool operator==(const Constant<String>& other) const {
        return (this->value()==other.value()); }
  private:
    String _value;
};

enum VariableType { type_bool, type_tribool, type_enumerated, type_string, type_integer, type_real };
enum VariableCategory { simple, dotted, primed };

//! \ingroup ExpressionModule
//! A named variable of unknown type.
class UntypedVariable {
  public:
    const Identifier& name() const { return *_name_ptr; }
    const VariableType& type() const { return this->_type; }
    bool operator==(const UntypedVariable& other) const {
        return (this->name()==other.name()) && (this->_category==other._category); }
    bool operator!=(const UntypedVariable& other) const { return !(*this==other); }
    bool operator<(const UntypedVariable& other) const {
        return this->name()<other.name() || (this->name()==other.name() && this->_type < other._type); }
    virtual std::ostream& write(std::ostream&) const;
  public:
    static std::string name(const VariableType& tp) {
        switch(tp) {
            case type_bool: return "Bool";
            case type_tribool: return "Tribool";
            case type_enumerated: return "Enumerated";
            case type_string: return "String";
            case type_integer: return "Integer";
            case type_real: return "Real";
        }
        return "Unknown";
    }
  protected:
    explicit UntypedVariable(const String& nm, VariableType tp, VariableCategory cat=simple)
        : _name_ptr(new Identifier(nm)), _type(tp), _category(cat) { }
  private:
    shared_ptr<Identifier> _name_ptr;
    VariableType _type;
    VariableCategory _category;
};

template<class T> inline VariableType variable_type() { ARIADNE_FAIL_MSG("Unknown variable type"); }
template<> inline VariableType variable_type<Boolean>() { return type_bool; }
template<> inline VariableType variable_type<Tribool>() { return type_tribool; }
template<> inline VariableType variable_type<String>() { return type_string; }
template<> inline VariableType variable_type<Integer>() { return type_integer; }
template<> inline VariableType variable_type<Real>() { return type_real; }

inline std::ostream& UntypedVariable::write(std::ostream& os) const {
    switch(this->_category) {
        case simple: os << this->name(); break;
        case dotted: os << "dot("<<this->name()<<")"; break;
        case primed: os << "next("<<this->name()<<")"; break;
    }
    //os << ":" << name(this->_type);
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const UntypedVariable& var) {
    return var.write(os); }



//! \ingroup ExpressionModule
//! A named variable of type \a T, possibly decorated by a "dot" or "prime"
//! representing a time derivative or updated value.
template<class T> class ExtendedVariable
    : public UntypedVariable
{
  public:
    Assignment< ExtendedVariable<T>,Expression<T> > operator=(const Expression<T>& e) const;
  protected:
    explicit ExtendedVariable(const String& nm, VariableCategory cat=simple)
        : UntypedVariable(nm, variable_type<T>(), cat) { }
};

template<class R, class D, class T=void> struct EnableIfRealDouble { };
template<class T> struct EnableIfRealDouble<Real,double,T> { typedef T Type; };

//! \ingroup ExpressionModule
//! A named variable of type \a T.
template<class T> class Variable
    : public ExtendedVariable<T>
{
    typedef Assignment< Variable<T>, Expression<T> > AssignmentType;
  public:
    explicit Variable(const String& nm) : ExtendedVariable<T>(nm) { }
    Variable<T> const& base() const { return *this; }
    inline AssignmentType operator=(const T& e) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& e) const;
    inline AssignmentType operator=(const Expression<T>& e) const;
    template<class D> inline typename EnableIfRealDouble<T,D,AssignmentType>::Type operator=(D e) const;
};



template<class T> LetVariable<T> let(const Variable<T>&);

template<class T> class LetVariable
    : public ExtendedVariable<T>
{
    typedef Assignment< Variable<T>, Expression<T> > AssignmentType;
  public:
    friend LetVariable<T> let<>(const Variable<T>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& val) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& var) const;
    inline AssignmentType operator=(const Expression<T>& expr) const;
    template<class D> inline typename EnableIfRealDouble<T,D,AssignmentType>::Type operator=(D e) const;
  private:
    explicit LetVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),simple) { }
};

template<class T> inline LetVariable<T> let(const Variable<T>& var) {
    return LetVariable<T>(var); }



DottedVariable<Real> dot(const Variable<Real>&);

template<class T> class DottedVariable;

template<class T> class DottedVariable
    : public ExtendedVariable<T>
{
    typedef Assignment< DottedVariable<T>, Expression<T> > AssignmentType;
  public:
    friend DottedVariable<Real> dot(const Variable<Real>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& e) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& e) const;
    inline AssignmentType operator=(const Expression<T>& e) const;
    template<class D> inline typename EnableIfRealDouble<T,D,AssignmentType>::Type operator=(D e) const;
  private:
    explicit DottedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),dotted) { }
};

inline DottedVariable<Real> dot(const Variable<Real>& var) {
    return DottedVariable<Real>(var); }


template<class T> PrimedVariable<T> next(const Variable<T>&);

template<class T> class PrimedVariable
    : public ExtendedVariable<T>
{
    typedef Assignment< PrimedVariable<T>, Expression<T> > AssignmentType;
  public:
    friend PrimedVariable<T> next<>(const Variable<T>&);
    Variable<T> base() const { return Variable<T>(this->name()); }
    inline AssignmentType operator=(const T& val) const;
    inline AssignmentType operator=(const Constant<T>& cnst) const;
    inline AssignmentType operator=(const Variable<T>& var) const;
    inline AssignmentType operator=(const Expression<T>& expr) const;
    template<class D> inline typename EnableIfRealDouble<T,D,AssignmentType>::Type operator=(D e) const;
  private:
    explicit PrimedVariable(const Variable<T>& var) : ExtendedVariable<T>(var.name(),primed) { }
};

template<class T> inline PrimedVariable<T> next(const Variable<T>& var) {
    return PrimedVariable<T>(var); }

template<class T> inline List< Variable<T> > operator,(const Variable<T>& v1, const Variable<T>& v2) {
    List< Variable<T> > r; r.append(v1); r.append(v2); return r; }


} // namespace Ariadne

#endif /* ARIADNE_VARIABLES_H */
