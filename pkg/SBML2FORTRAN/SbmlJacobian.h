/**
 * Copyright (C) 2015 Zuse Institute Berlin
 *
 * @file    SbmlJacobian.h
 *
 * @brief   Data structure to capture the important bits of an 
 *          SBML Document for eventually even faster numerics
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    11.03.2015
 *
 */
#ifndef __SBML_JACOBIAN_H
#define __SBML_JACOBIAN_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sbml/SBMLTypes.h>

#include "myTypes.h"

//


using namespace std;
LIBSBML_CPP_NAMESPACE_USE


//

typedef  pair< unsigned, unsigned >                     UnsignedPair;
typedef  pair< string, string >                         StringPair;

typedef  map< UnsignedPair, StringPair >                ArrayIndex;
typedef  map< string, map< string, string > >           StringArray;

/*
typedef  map< string, vector<string> >        FunctionDefList;
typedef  map< string, vector<string> >        EventAssignList;
typedef  map< string, vector<int> >           StoichiometryMat;
*/

//


class SbmlJacobian
{
    public:
        SbmlJacobian();
	SbmlJacobian(IndexList const& ispec,
	             IndexList const& iparm,
	             IndexList const& ifunc,
	             IndexList const& irule);
        /// SbmlJacobian(Model const* m);

        ~SbmlJacobian();

        void               set_druledy(StringList&      rule, 
                                       IndexList const& irule,
                                       IndexList const& ispec);
        StringArray const& get_druledy() const;
        ArrayIndex const& get_vdruledy() const;

        void               set_druledp(StringList&      rule,
                                       IndexList const& irule,
                                       IndexList const& iparm);
        StringArray const& get_druledp() const;
        ArrayIndex const& get_vdruledp() const;

        //

        void               set_dreacdy(StringList&      reac,
                                       IndexList const& ireac,
                                       IndexList const& ispec);
        StringArray const& get_dreacdy() const;
        ArrayIndex const& get_vdreacdy() const;

        void               set_dreacdp(StringList&      reac,
                                       IndexList const& ireac,
                                       IndexList const& iparm);
        StringArray const& get_dreacdp() const;
        ArrayIndex const& get_vdreacdp() const;

        //

        /// void           set_dfdy(StringList const& rate);
        void               set_dfdy(StringList&      rate,
                                    IndexList const& ispec);
        StringArray const& get_dfdy() const;
        ArrayIndex const& get_vdfdy() const;

        void               set_dfdp(StringList&      rate,
                                    IndexList const& ispec,
                                    IndexList const& iparm);
        StringArray const& get_dfdp() const;
        ArrayIndex const& get_vdfdp() const;

        SbmlJacobian const& toFortran(string const& fname =
                                          "JAC_LIMEX_TEMPLATE.txt");
    private:
	void initJacobian();
	void analyseFormula(string const& str, string const rId = "");
        string get_pId(string const& tok, 
                       string const& rId,
                       IndexList const& iparm) const;

    private:
	IndexList    _ispec, _iparm, _ifunc, _irule;

	StringArray  _druledy;
        StringArray  _druledp;

        StringArray  _dreacdy;
        StringArray  _dreacdp;

        StringArray  _dfdy;
        StringArray  _dfdp;

        ArrayIndex  _vdruledy;
        ArrayIndex  _vdruledp;

        ArrayIndex  _vdreacdy;
        ArrayIndex  _vdreacdp;

        ArrayIndex  _vdfdy;
        ArrayIndex  _vdfdp;

    public:
        friend ostream& operator<< (ostream& os, SbmlJacobian const& modeljac);

};

//

ostream& operator<< (ostream& os, SbmlJacobian const& modeljac);

ostream& operator<< (ostream& os, 
                     pair<StringArray,ArrayIndex> const& parray);

#endif // __SBML_JACOBIAN_H
