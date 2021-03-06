// Copyright (C) 2010 - 2015 Zuse Institute Berlin
//
// first added : 2010-11-29 td
// last changed:
//
#ifndef __EXPR_NODE_H
#define __EXPR_NODE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "myTypes.h"
#include "Expression.h"


namespace PARKIN
{
    class ExprNode
    {
            friend class Expression;

        protected:
            typedef Expression::Param           Param;
            // typedef Expression::ExprNodeType    ExprNodeType;

            ExprNode() : _count(1) {}
            virtual ~ExprNode() {}

            virtual char        type()                   const = 0;
            virtual Real        eval(Param&)             const = 0;
            virtual Real        eval(double*[])          const = 0;
            virtual Expression  df  (std::string const&) const = 0;
            virtual void        prt (std::ostream&)      const = 0;
            virtual bool        eq  (Real)               const = 0;
            virtual void        off (std::string const&, long const) = 0;
            virtual void        off (Expression::Names const&,
                                      std::vector<long> const&) { }

        private:
            long _count;
    };

    ///
    ///
    ///

    class Literal : public ExprNode
    {
        public:
            Literal(Real v) : _val(v) { }
            // virtual ~Literal() { }

            virtual char        type()                         const;
            virtual Real        eval(Param& /*x*/)             const;
            virtual Real        eval(double* /*y*/[])          const;
            virtual Expression  df  (std::string const& /*n*/) const;
            virtual void        prt (std::ostream& s)          const;
            virtual bool        eq  (Real v)                   const;
            virtual void        off (std::string const& /*n*/, long const /*a*/);

        private:
            Real _val;
    };

    class Variable : public ExprNode
    {
        public:
            Variable(std::string n) : _name(n), _off(0) { }
            // 24.07.12 td: used to be:  Variable(std::string n) : _name(n), _off(-1) { }
            // virtual ~Variable() { }

            virtual char        type()                     const;
            virtual Real        eval(Param& x)             const;
            virtual Real        eval(double* y[])          const;
            virtual Expression  df  (std::string const& n) const;
            virtual void        prt (std::ostream& s)      const;
            virtual bool        eq  (Real v)               const;
            virtual void        off (std::string const& n, long const a);

        private:
            std::string _name;
            long        _off;
    };

    ///

    class UnaryExpr : public ExprNode
    {
            // static std::map<std::string, long> _opMap;
            static std::map<ExprNodeType, std::string> _opMap;

        public:
            UnaryExpr(ExprNodeType const& op, Expression const& e) :
                _op(op), _expr(e)
            { createOpMap(); }
            // virtual ~UnaryExpr() { }

            virtual char        type()                     const;
            virtual Real        eval(Param& x)             const;
            virtual Real        eval(double* y[])          const;
            virtual Expression  df  (std::string const& n) const;
            virtual void        prt (std::ostream& s)      const;
            virtual bool        eq  (Real v)               const;
            virtual void        off (std::string const& n, long const a);

        private:
//            enum { UMINUS, EXP, LOG, SIN, COS, TAN };
            static void createOpMap()
            {
                _opMap[MINUS]   = "-";
                _opMap[ABS]     = "abs";
                _opMap[ARCCOS]  = "acos";
                _opMap[ARCOSH]  = "acosh";
                _opMap[ARCSIN]  = "asin";
                _opMap[ARSINH]  = "asinh";
                _opMap[ARCTAN]  = "atan";
                _opMap[ARTANH]  = "atanh";
                _opMap[CEIL]    = "ceil";
                _opMap[COS]     = "cos";
                _opMap[COSH]    = "cosh";
                _opMap[DELTA]   = "delta";
                _opMap[EXP]     = "exp";
                _opMap[FLOOR]   = "floor";
                _opMap[LN]      = "log";
                _opMap[LOG]     = "log10";
                _opMap[SIGN]    = "sign";
                _opMap[SIN]     = "sin";
                _opMap[SINH]    = "sinh";
                _opMap[TAN]     = "tan";
                _opMap[TANH]    = "tanh";
            }
            Real opEval(Real const& v) const;

            ExprNodeType    _op;
            Expression      _expr;
    };

    //

    class BinaryExpr : public ExprNode
    {
        public:
            BinaryExpr(ExprNodeType const& op, Expression const& e1, Expression const& e2) :
                _op(op), _lexpr(e1), _rexpr(e2)
            { }
            // virtual ~BinaryExpr() { }

            virtual char        type()                     const;
            virtual Real        eval(Param& x)             const;
            virtual Real        eval(double* y[])          const;
            virtual Expression  df  (std::string const& n) const;
            virtual void        prt (std::ostream& s)      const;
            virtual bool        eq  (Real v)               const;
            virtual void        off (std::string const& n, long const a);

        private:
            ExprNodeType    _op;
            Expression      _lexpr;
            Expression      _rexpr;
    };

    //

    class TernaryExpr : public ExprNode
    {
            static std::map<ExprNodeType, std::string> _opMap;

        public:
            TernaryExpr(ExprNodeType const& op,
                        Expression const& e1,
                        Expression const& e2,
                        Expression const& e3) :
                _op(op), _lexpr(e1), _mexpr(e2), _rexpr(e3)
            { createOpMap(); }
            // virtual ~TernaryExpr() { }

            virtual char        type()                     const;
            virtual Real        eval(Param& x)             const;
            virtual Real        eval(double* y[])          const;
            virtual Expression  df  (std::string const& n) const;
            virtual void        prt (std::ostream& s)      const;
            virtual bool        eq  (Real v)               const;
            virtual void        off (std::string const& n, long const a);

        private:
            static void createOpMap()
            {
                _opMap[TPLUS]     = "+++";
                _opMap[HILLplus]  = "Hplus";
                _opMap[HILLminus] = "Hminus";
            }

            ExprNodeType    _op;
            Expression      _lexpr;
            Expression      _mexpr;
            Expression      _rexpr;
    };

    //

    class UserExpr : public ExprNode
    {
        public:
            UserExpr(std::string const& name, Expression const& e) :
                _name(name), _expr(e)
            { }
            // virtual ~UserExpr() { }

            virtual char        type()                     const;
            virtual Real        eval(Param& x)             const;
            virtual Real        eval(double* y[])          const;
            virtual Expression  df  (std::string const& n) const;
            virtual void        prt (std::ostream& s)      const;
            virtual bool        eq  (Real v)               const;
            virtual void        off (std::string const& n, long const a);

        private:
            std::string     _name;
            Expression      _expr;
    };

    ///

    class Polynomial : public ExprNode
    {
        public:
            Polynomial(Expression::Names const&                 names,
                       Expression::Coeff const&                 coeff,
                       std::vector<Expression::Multix> const&  multIdx) :
                    _N(names.size()), _names(names), _offs(),
                    _Jcoeff(coeff.size()), _coeff(coeff), _multIdx(multIdx)
            { _offs.resize(_N); }
            // virtual ~Polynomial() {}

            virtual char        type()                      const;
            virtual Real        eval(Param& x)              const;
            virtual Real        eval(double* y[])           const;
            virtual Expression  df  (std::string const& n)  const;
            virtual void        prt (std::ostream& s)       const;
            virtual bool        eq  (Real v)                const;
            virtual void        off (std::string const&        n,
                                     long const                a);
            virtual void        off (Expression::Names const&  nVec,
                                     std::vector<long> const&  aVec);

        private:
            long                   _N;  // _N == _names.size() == _multIdx[j].size()
            Expression::Names       _names;
            std::vector<long>     _offs;

            long                               _Jcoeff;
            Expression::Coeff                   _coeff;
            std::vector<Expression::Multix>   _multIdx;
    };

    ///

    Real sign(Real x);
}

#endif // __EXPR_NODE_H
