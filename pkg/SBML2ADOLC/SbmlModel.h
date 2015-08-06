/**
 * Copyright (C) 2014 - 2015 Zuse Institute Berlin
 *
 * @file    SbmlModel.h
 *
 * @brief   Data structure to capture the important bits of an 
 *          SBML Document for eventually fast numerics
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    30.06.2015
 *
 */
#ifndef __SBML_MODEL_H
#define __SBML_MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sbml/SBMLTypes.h>

#include "myTypes.h"
// #include "SbmlJacobian.h"

//


using namespace std;
LIBSBML_CPP_NAMESPACE_USE


//

/*
typedef  map< unsigned, string >              ListIndex;
typedef  map< string, unsigned >              IndexList;
typedef  map< string, double >                ValueList;
typedef  map< string, string >                StringList;

// typedef  map< string, double >                CompartList;
// typedef  map< string, double >                SpeciesList;
// typedef  map< string, double >                ParameterList;
typedef  map< string, vector<string> >        FunctionDefList;
// typedef  map< string, string >                AssignRuleList;
// typedef  map< string, string >                ReactionList;
// typedef  map< string, string >                RateOdeList;
// // typedef  map< string, string >                TriggerList;
// // typedef  map< string, vector<string> >        EventMathList;
typedef  map< string, vector<string> >        EventAssignList;
typedef  map< string, vector<int> >           StoichiometryMat;
*/

//


class SbmlModel
{

    public:
        /// SbmlModel()  { }
        SbmlModel(Model const* m, string const& fname);
        SbmlModel(Model const* m, string const& fname, bool adolc);

        ~SbmlModel();

        SbmlModel(SbmlModel const& other);

        ///

        string const& getName() const;

        void setModelIDs(Model const* m);

        void setCompart(Model const* m);
        void setCompart(Compartment const* c, unsigned n = 0);

        void setSpecies(Model const* m);
        void setSpecies(Species const* s, unsigned n = 0);

        void setParameter(Model const* m);
        void setParameter(Reaction const* r);
        void setParameter(Parameter const* p, string const& prefix = "global");

        void setFunctionDef(Model const* m);
        void setFunctionDef(FunctionDefinition const* fd, unsigned n = 0);

        void setAssignRule(Model const* m);
        void setAssignRule(Rule const* r, unsigned n = 0);

        void setReaction(Model const* m);
        void setReaction(Reaction const* r, unsigned n = 0);

        void setEvent(Model const* m);
        void setEvent(Event const* e, unsigned n = 0);
        // void setEventAssign(EventAssignment const* ea, unsigned n = 0);
        
        void setRateOde(Model const* m);

        void translateEventAssignVector(vector<string>& vec);
        void translateFunctionDefVector(vector<string>& vec);
        void translateFormulaString(string& str, string const& rId = ""); 
        SbmlModel const& toAdolC(string const& fname = 
                                           "YDOT_LIMEXcpp_TEMPLATE.txt");

    private:
        SbmlModel const& operator= (SbmlModel const) { }
        void initModel(Model const* m);

        void doEventBlock(ifstream& infile);
        void doFunctionBlock(ifstream& infile);
        void doPiecewiseBlock(ifstream& infile);
        void doMultiLine(string const& str);
    
    private:
        // SbmlJacobian*    _jac;
        ostringstream    _outs;
        bool             _adolc;
        string           _sbmlfile;
        string           _modelname;
        ValueList        _comp, _spec, _parm;
        // ValueList        _bcsp, _cons;
        FunctionDefList  _func;
        StringList       _rule, _reac, _rate;
        StringList       _trig, _dmmy;
        EventAssignList  _emth;
        IndexList        _icomp, _ispec, _iparm, _ifunc;
        IndexList        _irule, _ireac;    // , _irate;
        IndexList        _itrig; // , _iemth;
        ListIndex        _vmodel;
        ListIndex        _vcomp, _vspec, _vparm, _vfunc;
        ListIndex        _vrule, _vreac;    // , _vrate;
        ListIndex        _vtrig; // , _vemth;
	IndexList        _piece;

    public:
        friend ostream& operator<< (ostream& os, SbmlModel const& model);

};


//


ostream& operator<< (ostream& os, SbmlModel const& model);

ostream& operator<< (ostream& os, pair<ValueList,ListIndex> const& plist);
ostream& operator<< (ostream& os, pair<ValueList,IndexList> const& plist);
ostream& operator<< (ostream& os, ValueList const& list);

ostream& operator<< (ostream& os, pair<StringList,ListIndex> const& plist);
ostream& operator<< (ostream& os, pair<StringList,IndexList> const& plist);
ostream& operator<< (ostream& os, StringList const& list);

ostream& operator<< (ostream& os, pair<FunctionDefList,ListIndex> const& plist);
ostream& operator<< (ostream& os, pair<FunctionDefList,IndexList> const& plist);
ostream& operator<< (ostream& os, FunctionDefList const& list);

ostream& operator<< (ostream& os, vector<string> const& vec);


//


vector<string> split(string const& str, string const& delim = " ,()+*/");
                              // '-' not(!!) in delim as it could be unary...
void replaceAll(string& str, string const& from, string const& to);
size_t replaceSingle(string& str, string const& from, string const& to, 
                     size_t start_pos = 0);


#endif // __SBML_MODEL_H
