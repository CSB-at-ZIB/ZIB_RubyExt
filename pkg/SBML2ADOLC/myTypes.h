/**
 * Copyright (C) 2015 Zuse Institute Berlin
 *
 * @file    myTypes.h
 *
 * @brief   Common typedefs for SBML translator necessary
 *          (so far) for the classed SbmlModel and SbmlJacobian
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    11.03.2015
 *
 */
#ifndef __MY_TYPES_H
#define __MY_TYPES_H

#include <string>
#include <vector>
#include <map>

//

using namespace std;

//

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

//

namespace PARKIN
{
    typedef  double  Real;
}

#endif // __MY_TYPES_H
