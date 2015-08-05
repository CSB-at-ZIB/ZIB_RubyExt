/**
 * Copyright (C) 2014 - 2015 Zuse Institute Berlin
 *
 * @file    SbmlModel.cpp
 *
 * @brief   Implementation of the data structure, 
 *          as declared in SbmlModel.h
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    20.05.2014
 *
 */

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include "SbmlModel.h"

/// using namespace std;

//=========================================================================

SbmlModel::SbmlModel(Model const* m, string const& fname) :
    // _jac(new SbmlJacobian()),
    _fortran(true), 
    _sbmlfile(fname),
    _modelname( m->getId() )
{
   // cerr << "c'tor SbmlModel " << this << endl;
 
   initModel(m);
}  

SbmlModel::SbmlModel(Model const* m, string const& fname, bool fortran) :
    // _jac(new SbmlJacobian()),
    _fortran(fortran),
    _sbmlfile(fname),
    _modelname( m->getId() )
{
   // cerr << "c'tor SbmlModel " << this << endl;

   initModel(m);
}

SbmlModel::~SbmlModel()
{
   // delete _jac;

   // cerr << "d'tor SbmlModel " << this << endl;
}

SbmlModel::SbmlModel(SbmlModel const& other)
{
  // _jac = new SbmlJacobian();

  // *_jac = *(other._jac);
  _outs.str( string() );
  _outs.clear();

  _fortran = other._fortran;
  _sbmlfile = other._sbmlfile;
  _modelname = other._modelname;

  _comp = other._comp; 
  _spec = other._spec; 
  _parm = other._parm;
  _func = other._func;
  _rule = other._rule; 
  _reac = other._reac; 
  _rate = other._rate;
  _trig = other._trig;
  _emth = other._emth;

  _icomp = other._icomp; 
  _ispec = other._ispec; 
  _iparm = other._iparm;
  _ifunc = other._ifunc;
  _irule = other._irule; 
  _ireac = other._ireac; 
  _itrig = other._itrig;

  _vcomp = other._vcomp; 
  _vspec = other._vspec; 
  _vparm = other._vparm;
  _vfunc = other._vfunc;
  _vrule = other._vrule; 
  _vreac = other._vreac; 
  _vtrig = other._vtrig;

  _piece = other._piece;
}

//=========================================================================

void 
SbmlModel::initModel(Model const* m)
{
   _piece.clear();

   setModelIDs(m);
   setCompart(m);
   setSpecies(m);
   setParameter(m);
   setFunctionDef(m);
   setAssignRule(m);
   setReaction(m);
   setEvent(m);
   setRateOde(m);

   /*
   if ( _jac != 0 )
   {
      _jac -> set_druledy(_rule,_irule,_ispec);
      _jac -> set_druledp(_rule,_irule,_iparm);

      _jac -> set_dreacdy(_reac,_ireac,_ispec);
      _jac -> set_dreacdp(_reac,_ireac,_iparm);

      _jac -> set_dfdy(_rate,_ispec);
      _jac -> set_dfdp(_rate,_ispec,_iparm);
   }
   */
}

//=========================================================================

string const&
SbmlModel::getName() const
{
  return _modelname;
}

//=========================================================================

void
SbmlModel::setModelIDs(Model const* m)
{
  if (m != 0)
  {
    char   str1[64],str2[64],str3[64],str4[64],str5[64];
    time_t now = time(0);
    string mId = m->getId();
    string fn = _sbmlfile;

    if ( mId.length() > 42 )
    {
       mId[43] = '\0';
    }

    if ( fn.length() > 42 )
    {
      fn[43] = '\0';
    }

    sprintf( str1, "%s", mId.c_str() );
    sprintf( str2, "%s", ctime(&now) );
    sprintf( str3, "%012lu", (unsigned long) now );
    sprintf( str4, "%s", fn.c_str() );
    sprintf( str5, "%s", "no vareq");

    str2[strlen(str2)-1] = '\0';

    _vmodel[1] = str1;
    _vmodel[2] = str2;
    _vmodel[3] = str3;
    _vmodel[4] = str4;
    _vmodel[5] = str5;
  }
}

//=========================================================================

void
SbmlModel::setCompart(Model const* m)
{
   if (m != 0) 
   {
      Compartment const* c;
      unsigned int       nCompart;

      nCompart = m->getNumCompartments();

      for (unsigned n = 0; n < nCompart; ++n)
      {
         c = m->getCompartment(n);
         setCompart(c,n+1);
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setCompart(Compartment const* c, unsigned n)
{
   if (c != 0)
   {
      string cId = c->getId();

      _icomp[cId] = n; 
      _vcomp[n] = cId;

      _comp[cId] = c->getVolume();
   }
}

//=========================================================================

void
SbmlModel::setSpecies(Model const* m)
{
   if (m != 0)
   {
      Species const* s;
      unsigned int   nSpecies;

      nSpecies = m->getNumSpecies();

      for (unsigned n = 0; n < nSpecies; ++n)
      {
         s = m->getSpecies(n);
         setSpecies(s,n+1);
      } 
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setSpecies(Species const* s, unsigned n)
{
   if (s != 0)
   {
      double val = 0.0;
      string sId = s->getId();

      if ( s->getBoundaryCondition() )
      {
         // _vspec[n] = sId;
         _rate[sId] = "0.0";
         return;
      }

      // if ( s->getHasOnlySubstanceUnits() || 
      //      s->getSubstanceUnits().size() > 0 )
      if ( s->isSetInitialAmount() )
      {
         val = s->getInitialAmount();
      }
      if ( s->isSetInitialConcentration() )
      {
         val = s->getInitialConcentration();
      }

      _ispec[sId] = n;
      _vspec[n] = sId;

      _spec[sId] = val;
   }
}

//=========================================================================

void
SbmlModel::setParameter(Model const* m)
{
   if (m != 0)
   {
      Parameter const* p;
      unsigned int     nParameter;

      nParameter = m->getNumParameters();

      for (unsigned n = 0; n < nParameter; ++n)
      {
         p = m->getParameter(n);
         setParameter(p);
      }

      Reaction const* r;
      unsigned int    nReaction;

      nReaction = m->getNumReactions();

      for (unsigned n = 0; n < nReaction; ++n)
      {
         r = m->getReaction(n);
         setParameter(r);
      }
   } 
}

//-------------------------------------------------------------------------

void
SbmlModel::setParameter(Reaction const* r)
{
   if (r != 0)
   {
      if ( r->isSetKineticLaw() )
      {
         Parameter const*  p;
         KineticLaw const* kl;
         unsigned int      nReactParameter;

         kl = r->getKineticLaw();
         nReactParameter = kl->getNumParameters();
         
         for (unsigned n = 0; n < nReactParameter; ++n)
         {
            p = kl->getParameter(n);
            setParameter( p, r->getId() );
         }

         unsigned int nLocalParameter;

         nLocalParameter = kl->getNumLocalParameters();

         for (unsigned n = 0; n < nLocalParameter; ++n)
         {
            p = kl->getLocalParameter(n);
            setParameter( p, r->getId()+"_local" );
         }
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setParameter(Parameter const* p, string const& prefix)
{
   if (p != 0)
   {
      unsigned n = _iparm.size() + 1;
      string   pId = prefix + "_" + p->getId();

      if (_iparm.count(pId) <= 0)
      {
         _iparm[pId] = n;
         _vparm[n] = pId;
      }

      _parm[pId] = p->getValue();
   }
}

//=========================================================================

void
SbmlModel::setFunctionDef(Model const* m)
{
   if (m != 0)
   {
      FunctionDefinition const* fd;
      unsigned int              nFunction;

      nFunction = m->getNumFunctionDefinitions();

      for (unsigned n = 0; n < nFunction; ++n)
      {
         fd = m->getFunctionDefinition(n);
         setFunctionDef(fd,n+1);
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setFunctionDef(FunctionDefinition const* fd, unsigned n)
{
   if (fd != 0)
   {
      char*          formula;
      ASTNode const* math;
      unsigned int   nChildren;
      vector<string> vec;

      vec.clear();

      math = fd->getMath();
      nChildren = math->getNumChildren();

      if ( nChildren > 1 )
      {
         for (unsigned j = 0; j < nChildren-1; ++j)
         {
            vec.push_back( (math->getChild(j))->getName() );
         }
      }

      if ( nChildren > 0 )
      {
         formula = SBML_formulaToString(math->getChild(nChildren-1));

         vec.push_back( formula );

         free(formula);
      }

      string fId = fd->getId();

      _ifunc[fId] = n;
      _vfunc[n] = fId;

      if (_fortran)
      {
         translateFunctionDefVector(vec);
      }

      _func[fId] = vec;
   }
}

//=========================================================================

void
SbmlModel::setAssignRule(Model const* m)
{
   if (m != 0)
   {
      Rule const*  r;
      unsigned int nRule;

      nRule = m->getNumRules();

      for (unsigned n = 0; n < nRule; ++n)
      {
         r = m->getRule(n);
         setAssignRule(r,n+1);
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setAssignRule(Rule const* r, unsigned n)
{
   if (r != 0)
   {
      char*  formula;
      string rId;

      rId = r->getVariable();

      if ( r->isSetMath() )
      {
         formula = SBML_formulaToString( r->getMath() );

         if (rId.length() <= 0)
         {
            char str[64];

            sprintf(str, "zerorule%u", (unsigned)(_rule.size()+1));
            rId = str;
         }

         _irule[rId] = n;
         _vrule[n] = rId;

         _rule[rId] = formula;

         free(formula);
      }

      if (_fortran)
      {
         translateFormulaString(_rule[rId]);
      }
   }
}

//=========================================================================

void
SbmlModel::setReaction(Model const* m)
{
   if (m != 0)
   {
      Reaction const* r;
      unsigned int    nReaction;

      nReaction = m->getNumReactions();

      for (unsigned n = 0; n < nReaction; ++n)
      {
         r = m->getReaction(n);
         setReaction(r,n+1);
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setReaction(Reaction const* r, unsigned n)
{
   if (r != 0)
   {
      string rId;

      rId = r->getId();

      _ireac[rId] = n;
      _vreac[n] = rId;

      _reac[rId] = "0.0";

      if ( r->isSetKineticLaw() )
      {
         KineticLaw const* kl;
         char*             formula;

         kl = r->getKineticLaw();

         if ( kl->isSetMath() )
         {
            formula = SBML_formulaToString( kl->getMath() );

            _reac[rId] = formula;

            free(formula);
         }
      }

      if (_fortran)
      {
         translateFormulaString(_reac[rId], rId);
      }
   }
}

//=========================================================================

void
SbmlModel::setEvent(Model const* m)
{
   if (m != 0)
   {
      Event const*  e;
      unsigned int  nEvent;

      nEvent = m->getNumEvents();

      for (unsigned n = 0; n < nEvent; ++n)
      {
         e = m->getEvent(n);
         setEvent(e,n+1);
      }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::setEvent(Event const* e, unsigned n)
{
   if (e != 0)
   {
      string          eId;
      unsigned int    nEventAssignment;
      char*           formula;
      vector<string>  vec;

      eId = e->getId();

      if (eId.length() <= 0)
      {
         char str[64];
   
         sprintf(str, "event%u", (unsigned)(_trig.size()+1));

         eId = str;
      }

      if (e->isSetDelay() )
      {
         formula = SBML_formulaToString( (e->getDelay())->getMath() );

         if ( strcmp(formula,"0") != 0 )
         { 
            free(formula);
            return;
         }

         free(formula);
      }

      _itrig[eId] = n;
      _vtrig[n] = eId;

      _trig[eId] = ".false.";

      if ( e->isSetTrigger() )
      {
         formula = SBML_formulaToString( (e->getTrigger())->getMath() );

         _trig[eId] = formula;

         free(formula);
      }

      //

      vec.clear();

      nEventAssignment = e->getNumEventAssignments();

      for (unsigned j = 0; j < nEventAssignment; ++j)
      {
         EventAssignment const* ea;

         ea = e->getEventAssignment(j);
         // setEventAssign(ea,j+1);

         if ( ea->isSetMath() )
         {
            vec.push_back( ea->getVariable() );

            formula = SBML_formulaToString( ea->getMath() );

            vec.push_back( formula );

            free(formula);
         }
      }

      if (_fortran)
      {
         translateFormulaString(_trig[eId]);

         translateEventAssignVector(vec);
      }

      _emth[eId] = vec;
   }
}

//=========================================================================

void
SbmlModel::setRateOde(Model const* m)
{
   if (m != 0)
   {
      StoichiometryMat stm;
      unsigned int     nReactant;
      unsigned int     nReaction;
      string           sId;
      Reaction const*  r;

      stm.clear();

      nReaction = m->getNumReactions();

      for (unsigned n = 0; n < nReaction; ++n)
      {
         r = m->getReaction(n);
         nReactant = r->getNumReactants();

         for (unsigned j = 0; j < nReactant; ++j)
         {
             sId = (r->getReactant(j))->getSpecies();

             if ( (m->getSpecies(sId))->getBoundaryCondition() )
             {
                continue;
             }

             if ( stm.count(sId) <= 0 )
             {
                stm[sId] = vector<int>(nReaction, 0);
             }

             stm[sId][n] -= (r->getReactant(j))->getStoichiometry(); 
         }

         unsigned nProduct;
         nProduct = r->getNumProducts();

         for (unsigned j = 0; j < nProduct; ++j)
         {
             sId = (r->getProduct(j))->getSpecies();

             if ( (m->getSpecies(sId))->getBoundaryCondition() )
             {
                continue;
             }

             if ( stm.count(sId) <= 0 )
             {
                stm[sId] = vector<int>(nReaction, 0);
             }

             stm[sId][n] += (r->getProduct(j))->getStoichiometry(); 
         } 
      }

      for (StoichiometryMat::const_iterator it = stm.begin(); 
                                            it != stm.end(); 
                                          ++it)
      {
          char   str[4096];
          string cId, eqn;

          cId = (m->getSpecies(it->first))->getCompartment();
          eqn.clear();

          if (!_fortran)
          {
             for (unsigned n = 0; n < nReaction; ++n)
             {
                int  j = (it->second)[n];

                if (j != 0)
                {
                    if (j == 1)
                    {
                       sprintf(str, " + %s", 
                               (m->getReaction(n))->getId().c_str() );
                    }
                    else if (j == -1)
                    {
                       sprintf(str, " - %s", 
                               (m->getReaction(n))->getId().c_str() );
                    }
                    else
                    {
                       char sign = (j < 0) ? '-' : '+';
   
                       sprintf(str, " %c (%d) * %s", 
                               sign, abs(j), 
                               (m->getReaction(n))->getId().c_str() );
                    }

                    eqn += str;
                }
             }

             sprintf(str, "(%s ) / %s", eqn.c_str(), cId.c_str() );
          }
          else
          {
             for (unsigned n = 0; n < nReaction; ++n)
             {
                int  j = (it->second)[n];

                if (j != 0)
                {
                    if (j == 1)
                    {
                       sprintf(str, " + rea(%d)", 
                               _ireac[(m->getReaction(n))->getId()] );
                    }
                    else if (j == -1)
                    {
                       sprintf(str, " - rea(%d)", 
                               _ireac[(m->getReaction(n))->getId()] );
                    }
                    else
                    {
                       char sign = (j < 0) ? '-' : '+';
   
                       sprintf(str, " %c (%d) * rea(%d)", 
                               sign, abs(j), 
                               _ireac[(m->getReaction(n))->getId()] );
                    }

                    eqn += str;
                }
             }

             sprintf(str, "(%s ) / com(%d)", eqn.c_str(), _icomp[cId] );
          }

          _rate[it->first] = str;
      }

      //

      for (StringList::const_iterator it = _rule.begin();
                                      it != _rule.end();
                                    ++it)
      {
          if (_fortran)
          {
             char str[4096];

             if ( _ispec.count(it->first) > 0 )
             {
                 sprintf(str, " + rul(%d) - spe(%d)", 
                         _irule[it->first], _ispec[it->first]);

                 _rate[_vspec[_ispec[it->first]]] = str;
             }
          }
      }
   }
}

//=========================================================================

void
SbmlModel::translateEventAssignVector(vector<string>& vec)
{
  if ( !vec.empty() )
  {
     for (unsigned n = 0; n < vec.size(); ++n)
     {
        translateFormulaString(vec[n]);
     }
  }
}

//-------------------------------------------------------------------------

void
SbmlModel::translateFunctionDefVector(vector<string>& vec)
{
   if ( !vec.empty() )
   {
      size_t               pos = 0;
      string               fun = vec.back();
      vector<string>       tok = split(fun);
      map<string,unsigned> vargs;

      vargs.clear();
      
      for (unsigned n = 0; n < vec.size()-1; ++n)
      {
         vargs[vec[n]] = n+1;
      }

      for (unsigned n = 0; n < tok.size(); ++n)
      {
         char   toF[64];
         string w = tok[n];

         char*  pEnd = &w[0];
         double val = strtod(w.c_str(), &pEnd);

         if ( *pEnd == '\0' )
         {
            if (val < 0.0)
            {
              sprintf(toF, "(%E)", val);
            }
            else
            {
              sprintf(toF, "%E", val);
            }
            pEnd = strchr(toF, 'E');
            *pEnd = 'D';
            pos = replaceSingle(fun, w, toF, pos);
         }
         else if ( vargs.count(w) > 0 )
         {
             sprintf(toF, "x%d", vargs[w]);
             pos = replaceSingle(fun, w, toF, pos);
         }
      }

      for (unsigned n = 0; n < vec.size()-1; ++n)
      {
         char toF[64];
         sprintf(toF, "x%d", n+1);
         vec[n] = toF;
      }

      vec.pop_back();
      vec.push_back( fun );
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::translateFormulaString(string& str, string const& rId)
{
   size_t         pos = 0;
   vector<string> tok = split(str);

   for (unsigned n = 0; n < tok.size(); ++n)
   {
      char   toF[64];
      string w = tok[n];

      char*  pEnd = &w[0];
      double val = strtod(w.c_str(), &pEnd);

      if ( *pEnd == '\0' )
      {
         if (val < 0.0)
         {
           sprintf(toF, "(%E)", val);
         }
         else
         {
           sprintf(toF, "%E", val);
         }
         pEnd = strchr(toF, 'E');
         *pEnd = 'D';
         pos = replaceSingle(str, w, toF, pos);
      }
      else
      { 
         if ( w[0] == '-' )  // get rid of unary minus!
         {
            w.erase(w.begin(),w.begin()+1);
            ++pos;
         }

         string ww = "global_" + w;
         string rw = rId + "_" + w;
         string lw = rId + "_local_" + w;

         if ( _irule.count(w) > 0 )
         {
            sprintf(toF, "rul(%d)", _irule[w]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _ifunc.count(w) > 0 )
         {
            sprintf(toF, "fun%d", _ifunc[w]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _iparm.count(lw) > 0 )  //first, check purely local parameters
         {
            sprintf(toF, "par(%d)", _iparm[lw]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _iparm.count(rw) > 0 )  // then, parameters based on reactions
         {
            sprintf(toF, "par(%d)", _iparm[rw]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _iparm.count(ww) > 0 )  // finally, global parameters
         {
            sprintf(toF, "par(%d)", _iparm[ww]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _ispec.count(w) > 0 )
         {
            sprintf(toF, "spe(%d)", _ispec[w]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( _icomp.count(w) > 0 )
         {
            sprintf(toF, "com(%d)", _icomp[w]);
            pos = replaceSingle(str, w, toF, pos);
         }
         else if ( w.compare("time") == 0 )
         {
            sprintf(toF, "t");
            pos = replaceSingle(str, w, toF, pos);
         }
      }
   }

   pos = str.find("piecewise");

   while( pos != string::npos )
   {
      size_t cur = pos+8;
      int    level = 0;
      int    nargs = 0;
      char   ch, toF[64];

      while ( 
              ( cur < str.size() ) && 
              ( (ch = str[++cur]) != ')' || level != 1 )
            )
      {
         if ( ch == '(' )
         {
            ++level;
         }
         else if ( ch == ')' )
         {
            --level;
            if ( level < 1 ) return;
         }
         else if ( str[cur] == ',' && level == 1 )
         {
            ++nargs;
         }
      }

      sprintf(toF, "piecewise%d", nargs);
      _piece[toF] = unsigned(nargs);
      pos = replaceSingle(str, "piecewise", toF, pos);
      pos = str.find("piecewise",pos);
   }
}

//-------------------------------------------------------------------------

SbmlModel const& 
SbmlModel::toFortran(string const& fname)
{
   ifstream infile;
   string   line;

   infile.open( fname.c_str() );

   _outs.str( string() );
   _outs.clear();

   if ( infile.is_open() )
   {
      while ( getline( infile, line ) )
      {
         char str[128];

         if ( line.size() <= 0 )
         {
            continue;
         }

         if ( line[0] == 'c' )
         {
            if ( line.find("c#block") != string::npos )
            {
               if ( line.find("pieces") != string::npos )
               {
                  doPiecewiseBlock( infile );
               }
               else if ( line.find("functions") != string::npos )
               {
                  doFunctionBlock( infile );
               }
               else if ( line.find("events") != string::npos )
               {
                  doEventBlock( infile );
               }
               line = "c ";
            }
            else if ( line.find("c#comment") != string::npos )
            {
               ListIndex::const_iterator lBeg;
               ListIndex::const_iterator lEnd;

               if ( line.find("modelids") != string::npos )
               {
                  lBeg = _vmodel.begin();
                  lEnd = _vmodel.end();
               }
               if ( line.find("compartments") != string::npos )
               {
                  lBeg = _vcomp.begin();
                  lEnd = _vcomp.end();
               }
               if ( line.find("species") != string::npos )
               {
                  lBeg = _vspec.begin();
                  lEnd = _vspec.end();
               }
               if ( line.find("parameters") != string::npos )
               {
                  lBeg = _vparm.begin();
                  lEnd = _vparm.end();
               }

               getline( infile, line );

               _outs << "c     ---" << endl;

               for (ListIndex::const_iterator it = lBeg;
                                              it != lEnd;
                                            ++it)
               {
                   sprintf( str, line.c_str(), 
                            it->first, (it->second).c_str() );

                   _outs << str << endl;
               }

               line = "c     ---";
            }
            else if ( line.compare("c#modelids") == 0 )
            {
               getline( infile, line );
               if ( !_vmodel.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vmodel.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#compartments") == 0 )
            {
               getline( infile, line );
               if ( !_vcomp.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vcomp.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#species") == 0 )
            {
               getline( infile, line );
               if ( !_vspec.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vspec.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#parameters") == 0 )
            {
               getline( infile, line );
               if ( !_vparm.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vparm.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#rules") == 0 )
            {
               getline( infile, line );
               if ( !_vrule.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vrule.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#reactions") == 0 )
            {
               getline( infile, line );
               if ( !_vreac.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vreac.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.compare("c#events") == 0 )
            {
               getline( infile, line );
               if ( !_vtrig.empty() )
               {
                  sprintf( str, line.c_str(), 
                           (_vtrig.rbegin())->first );
               }
               else
               {
                  sprintf( str, line.c_str(), 0 );
               } 
               line = str;
            }
            else if ( line.find("c#for_all") != string::npos )
            {
               ListIndex::const_iterator lBeg;
               ListIndex::const_iterator lEnd;
               ValueList&                val = _comp;
               StringList&               formula = _rule;
               char                      lstkind = 'v';

               if ( line.find("compartments") != string::npos )
               {
                  lBeg = _vcomp.begin();
                  lEnd = _vcomp.end();
                  val = _comp;
                  lstkind = 'v';
               }
               else if ( line.find("species") != string::npos )
               {
                  lBeg = _vspec.begin();
                  lEnd = _vspec.end();
                  val = _spec;
                  lstkind = 'v';
               }
               else if ( line.find("parameters") != string::npos )
               {
                  lBeg = _vparm.begin();
                  lEnd = _vparm.end();
                  val = _parm;
                  lstkind = 'v';
               }
               else if ( line.find("functions") != string::npos )
               {
                  lBeg = _vfunc.begin();
                  lEnd = _vfunc.end();
                  lstkind = 'f';
               }
               else if ( line.find("pieces") != string::npos )
               {
                  //lBeg = _pieces.begin();
                  // lEnd = _pieces.end();
                  lstkind = 'p';
               }
               else if ( line.find("algebraic") != string::npos )
               {
                  lBeg = _vspec.begin();
                  lEnd = _vspec.end();
                  formula = _rule;
                  lstkind = 'a';
               }
               else if ( line.find("assignments") != string::npos )
               {
                  lBeg = _vrule.begin();
                  lEnd = _vrule.end();
                  formula = _rule;
                  lstkind = 's';
               }
               else if ( line.find("reactions") != string::npos )
               {
                  lBeg = _vreac.begin();
                  lEnd = _vreac.end();
                  formula = _reac;
                  lstkind = 's';
               }
               // else if ( line.find("events") != string::npos )
               // {
               //    lBeg = _vtrig.begin();
               //    lEnd = _vtrig.end();
               //    formula = _trig;
               //    lstkind = 's';
               // }
               else if ( line.find("rates") != string::npos )
               {
                  lBeg = _vspec.begin();
                  lEnd = _vspec.end();
                  formula = _rate;
                  lstkind = 's';
               }
 
               getline( infile, line );

               string old = line;

               switch( lstkind )
               {
                  case 'p':  // piecewise decls :  'dble prec piecewise*'
                     for (IndexList::const_iterator it = _piece.begin(); 
                                                    it != _piece.end(); 
                                                  ++it)
                     {
                        sprintf( str, line.c_str(), it->second );
                        _outs << str << endl;
                     }
                     break;

                  case 'f':  // function decls :  'dble prec fun*'
                     for (ListIndex::const_iterator it = lBeg; 
                                                    it != lEnd; 
                                                  ++it)
                     {
                        sprintf( str, line.c_str(), it->first );
                        _outs << str << endl;
                     }
                     break;

                  case 'a':  // algebraic :  'b(*) = zero'
                     for (ListIndex::const_iterator it = lBeg; 
                                                    it != lEnd; 
                                                  ++it)
                     {
                        if ( formula.count(it->second) > 0 )
                        {
                            sprintf( str, line.c_str(), it->first );
                            _outs << str << endl;
                        }
                     }
                     break;
                     
                  case 'v':  // values :  'field(*) = dble prec value'
                     for (ListIndex::const_iterator it = lBeg; 
                                                    it != lEnd; 
                                                  ++it)
                     {
                         sprintf( str, line.c_str(), 
                                  it->first, val[it->second] );
                         char* pEnd = strchr(str, 'E');
                         *pEnd = 'D';
                         _outs << str << endl;
                     }
                     break;
                  
                  case 's':  // strings : 'field(*) = string formula'
                     for (ListIndex::const_iterator it = lBeg; 
                                                    it != lEnd; 
                                                  ++it)
                     {
                         char buf[8192]; 

                         sprintf( buf, line.c_str(), 
                                  it->first, formula[it->second].c_str() );

                         doMultiLine( buf );
/*
                         unsigned     len = 0;
                         stringstream sstr, ostr;
                         string       eqn;

                         ostr.flush();
                         sstr.flush();
                         sstr << formula[it->second];
                         len = 0;
                         line = old;

                         while (sstr >> skipws >> eqn)
                         {
                            len += eqn.size();

                            if ( ++len < 45 )
                            {
                                ostr << " " << eqn;
                            }
                            else
                            {
                                if ( line[5] != '&' )
                                {
                                    sprintf( str, line.c_str(), 
                                             it->first, ostr.str().c_str() );

                                    line = "     &          %s";
                                }
                                else
                                {
                                     sprintf( str, line.c_str(), 
                                              ostr.str().c_str() );
                                }

                                _outs << str << endl;

                                ostr.flush();
                                ostr.str( string() );
                                ostr.clear();
                                ostr << " " << eqn;
                                len = eqn.size()+1;
                            }
                         }

                         if ( len < 45 )
                         {
                            if ( line[5] != '&' )
                            {
                                sprintf( str, line.c_str(), 
                                         it->first, ostr.str().c_str() );
                            }
                            else
                            {
                                 sprintf( str, line.c_str(), 
                                          ostr.str().c_str() );
                            }

                            _outs << str << endl;
                         }
*/  
                     }
                     break;

                  default:
                     break;
               }

               line = "c ";
            }
         }

         _outs << line << endl; 
      }

      infile.close();
   }
   else
   {
      _outs << "Fatal error: .f template file not found!" << endl;
   }

   return *this;
}

//-------------------------------------------------------------------------

void
SbmlModel::doEventBlock(ifstream& infile)
{
   string         line;
   vector<string> block;

   block.clear();

   while ( getline( infile, line ) )
   {
      if ( line.compare("c#end_block") == 0 )
      {
         break;
      }

      block.push_back( line );
   }

   if ( line.compare("c#end_block") != 0 )
   {
      return;
   }

   for (ListIndex::const_iterator it = _vtrig.begin();
                                  it != _vtrig.end();
                                ++it)
   {
      char            str[4096];
      unsigned        eno = it->first;
      vector<string>  vec = _emth[it->second];

      for (unsigned n = 0; n < block.size(); ++n)
      {
         line = block[n];

         if ( line.compare("c#trigs") == 0 )
         {
            line = block[++n];

            sprintf( str, line.c_str(),
                     eno, _trig[it->second].c_str());

            doMultiLine( str );

            line = "c ";
         }
         else if ( line.compare("c#if_trig") == 0 )
         {
            line = block[++n];

            sprintf( str, line.c_str(), eno );

            line = str;
         }
         else if ( line.find("c#for_all") != string::npos )
         {
            line = block[++n];

            for (unsigned j = 0; j < vec.size(); j+=2)
            {
                sprintf( str, line.c_str(), 
                         vec[j].c_str(), vec[j+1].c_str() );

                doMultiLine( str );
            }

            line = "c ";
         }

         _outs << line << endl;
      }
   }
}

//-------------------------------------------------------------------------

void 
SbmlModel::doFunctionBlock(ifstream& infile)
{
   string         line;
   vector<string> block;

   block.clear();

   while ( getline( infile, line ) )
   {
      if ( line.compare("c#end_block") == 0 )
      {
         break;
      }

      block.push_back( line );
   }

   if ( line.compare("c#end_block") != 0 )
   {
       return;
   }

   for (ListIndex::const_iterator it = _vfunc.begin();
                                  it != _vfunc.end();
                                ++it)
   {
       char            str[4096];
       unsigned        fno = it->first;
       vector<string>  vec = _func[it->second];

       for (unsigned n = 0; n < block.size(); ++n)
       {
           line = block[n];

           if ( line.compare("c#args") == 0 )
           {
               string arg;

               arg.clear();

               if ( vec.size() > 1 )
               {
                   arg = "( " + vec[0];
                   for (unsigned j = 1; j < vec.size()-1; ++j)
                   { 
                       arg += ", " + vec[j];
                   }
                   arg += " )";
               }

               line = block[++n];

               sprintf( str, line.c_str(),
                        fno, arg.c_str());

               doMultiLine( str );

               line = "c ";
           }
           else if ( line.find("c#for_all") != string::npos )
           {
               line = block[++n];

               for (unsigned j = 0; j < vec.size()-1; ++j)
               {
                   sprintf( str, line.c_str(), vec[j].c_str() );

                   _outs << str << endl;
               }

               line = "c ";
           }
           else if ( line.compare("c#impls") == 0 )
           {
               line = block[++n];

               sprintf( str, line.c_str(), 
                        fno, vec.back().c_str() );

               doMultiLine( str );

               line = "c ";
           }
           else if ( line.compare("c#funcs") == 0 )
           {
               line = block[++n];

               sprintf( str, line.c_str(), fno );

               line = str;
           }

           _outs << line << endl;
       }
   }
}

//-------------------------------------------------------------------------

void
SbmlModel::doPiecewiseBlock(ifstream& infile)
{
   string         line;
   vector<string> block;

   block.clear();

   while ( getline( infile, line ) )
   {
      if ( line.compare("c#end_block") == 0 )
      {
         break;
      }

      block.push_back( line );
   }

   if ( line.compare("c#end_block") != 0 )
   {
      return;
   }
 
   for (IndexList::const_iterator it = _piece.begin();
                                  it != _piece.end();
                                ++it)
   {
      char      str[4096];
      unsigned  narg = _piece[it->first];

      for (unsigned n = 0; n < block.size(); ++n)
      {
         line = block[n];

         if ( line.compare("c#args") == 0 )
         {
               string arg;

               arg.clear();

               if ( narg%2 == 0 )
               {
                   arg = "( ";
                   for (unsigned j = 0; j < narg/2; ++j)
                   {
                       sprintf(str, "val%d, cond%d, ", j+1, j+1);
                       arg += str;
                   }
                   arg += "dflt )";
               }

               line = block[++n];

               sprintf( str, line.c_str(),
                        narg, arg.c_str());

               doMultiLine( str );

               line = "c ";
         }
         else if ( line.find("c#for_all") != string::npos )
         {
               char arg[32];

               if ( line.find("dbledecls") != string::npos )
               {
                    line = block[++n];

                    if ( narg%2 == 0 )
                    {
                       for (unsigned j = 0; j < narg/2; ++j)
                       {
                           sprintf(arg, "val%d", j+1);
                           sprintf(str, line.c_str(), arg);

                           _outs << str << endl;
                       }
                    }

                    sprintf(str, line.c_str(), "dflt");

                    _outs << str << endl;
               }
               else if ( line.find("logicdecls") != string::npos )
               {
                    line = block[++n];

                    if ( narg%2 == 0 ) 
                    {
                       for (unsigned j = 0; j < narg/2; ++j)
                       {
                           sprintf(arg, "cond%d", j+1);
                           sprintf(str, line.c_str(), arg);

                           _outs << str << endl;
                       }
                    }
               }
               else if ( line.find("remaining_cases") != string::npos )
               {
                    if ( narg%2 == 0 )
                    {
                       for (unsigned j = 1; j < narg/2; ++j)
                       {
                           char arg[32];

                           if ( j > 1 ) 
                           {
                              _outs << "c " << endl;
                           }
                           line = block[n+1];
                           sprintf(arg, "cond%d", j+1);
                           sprintf(str, line.c_str(), arg);
                           _outs << str << endl;
              
                           _outs << block[n+2] << endl;

                           line = block[n+3];
                           sprintf(arg, "val%d", j+1);
                           sprintf(str, line.c_str(), narg, arg);
                           _outs << str << endl;
                       } 
                    }
            
                    n += 3;
               }
               else
               {
                  if ( line.find("3lines") != string::npos ) n+=2;

                  ++n;
               }

               line = "c ";
         }
         else if ( line.compare("c#3lines_first_case") == 0 )
         {
            line = block[n+1];

            sprintf( str, line.c_str(), 
                     ((narg > 1) ? "cond1" : ".false.") );

            _outs << str << endl;
            _outs << block[n+2] << endl;

            line = block[n+3];

            sprintf( str, line.c_str(), 
                     narg, ((narg > 1) ? "val1" : "0.0D0") );
 
            _outs << str << endl;

            n += 3;

            line = "c ";
         }
         else if ( line.compare("c#3lines_default") == 0 )
         {
              _outs << block[n+1] << endl;
              _outs << block[n+2] << endl;

              line = block[n+3];

              sprintf(str, line.c_str(), narg, "dflt");

              line = str;

              n += 3;
         }
         else if ( line.compare("c#funcs") == 0 )
         {
               line = block[++n];

               sprintf( str, line.c_str(), narg );

               line = str;
         }

         _outs << line << endl;
      }
   }
}

//-------------------------------------------------------------------------

void 
SbmlModel::doMultiLine(string const& str)
{
   bool         firstline = true;
   unsigned     len = 0;
   stringstream sstr, ostr;
   string       bits;
   
   ostr.flush();
   sstr.flush();
   
   sstr << str;

   len = 0;
   firstline = true;

   while (sstr >> skipws >> bits)
   {
      len += bits.size();

      if ( ++len < 64 )
      {
          ostr << " " << bits;
      }
      else
      {
          len = 9;

          if ( firstline )
          {  //        12345
             _outs << "     ";
             firstline = false;
          }
          else
          {  //        123456789012345
             _outs << "     &         ";
          }

          _outs << ostr.str() << endl;

          ostr.flush();
          ostr.str( string() );
          ostr.clear();
          
          ostr << " " << bits;

          len += bits.size()+1;
      }
   }

   if ( len < 64 )
   {
       if ( firstline )
       {  //        12345
          _outs << "     ";
       }
       else
       {  //        123456789012345
          _outs << "     &         ";
       }

       _outs << ostr.str() << endl;
   }
}

//=========================================================================

//=========================================================================

ostream& operator<< (ostream& os, SbmlModel const& model)
{
   time_t now = time(0);

   if (model._fortran)
   {
      os << "c-----" << endl;
      os << "c SBML Model : " << setw(55) << left 
                              << model.getName() 
                              << endl;
      os << "c              " << setw( min(55,(int)model.getName().size()) ) 
                              << setfill('~') << "~" << endl; 
      os << "c       Date : " << ctime(&now);
      os << "c              " << endl;
      os << "c     Author : " << "automated transcription by 'sbml2fortran'"
                              << endl;
      os << "c              " << endl;
      os << "c Copyright (C) Zuse Institute Berlin, CSB Group" 
         << endl;
      os << "c-----" << endl;
      os << "c" << endl;
      os << model._outs.str();
   }
   else
   {
      os << pair<ValueList,ListIndex>( model._comp, model._vcomp );
      os << endl;
      os << pair<ValueList,ListIndex>( model._spec, model._vspec ); 
      os << endl;
      os << pair<ValueList,ListIndex>( model._parm, model._vparm ); 
      os << endl;

      os << pair<FunctionDefList,ListIndex>( model._func, model._vfunc );
      os << endl;

      os << pair<StringList,ListIndex>( model._rule, model._vrule );
      os << endl;
      os << pair<StringList,ListIndex>( model._trig, model._vtrig );
      os << endl;
      os << pair<StringList,ListIndex>( model._reac, model._vreac ); 
      os << endl;

      // os << model._rate;
      os << pair<StringList,ListIndex>( model._rate, model._vspec );
      os << endl;
 
      //

      /*
      os << pair<StringArray,ArrayIndex>( model._jac->get_druledy(),
                                          model._jac->get_vdruledy()
                                        );
      os << endl;
      os << pair<StringArray,ArrayIndex>( model._jac->get_druledp(),
                                          model._jac->get_vdruledp()
                                        );
      os << endl;
      os << pair<StringArray,ArrayIndex>( model._jac->get_dreacdy(),
                                          model._jac->get_vdreacdy()
                                        );
      os << endl;
      os << pair<StringArray,ArrayIndex>( model._jac->get_dreacdp(),
                                          model._jac->get_vdreacdp()
                                        );
      */
/*
      os << endl;
      // since model const, wrong here: model._jac->set_dfdy(model._rate);
      os << pair<StringArray,ArrayIndex>( model._jac->get_dfdy(),
                                          model._jac->get_vdfdy()
                                        );
      os << endl;
      os << pair<StringArray,ArrayIndex>( model._jac->get_dfdp(), 
                                          model._jac->get_vdfdp()
                                        );
*/
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<ValueList,ListIndex> const& plist)
{
   ValueList list = plist.first;
   ListIndex indx = plist.second;

   for (ListIndex::const_iterator it = indx.begin();
                                  it != indx.end();
                                ++it)
   {
      os << it->first << " : " << it->second << " : ";
      os << /* scientific << */ list[it->second] << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<StringList,ListIndex> const& plist)
{
   StringList list = plist.first;
   ListIndex  indx = plist.second;

   for (ListIndex::const_iterator it = indx.begin();
                                  it != indx.end();
                                ++it)
   {
      os << it->first << " : " << it->second << " : ";
      os << list[it->second] << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<FunctionDefList,ListIndex> const& plist)
{
   FunctionDefList list = plist.first;
   ListIndex       indx = plist.second;

   for (ListIndex::const_iterator it = indx.begin();
                                  it != indx.end();
                                ++it)
   {
      os << it->first << " : " << it->second << " : ";
      os << list[it->second] << endl;
   }

   return os;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<ValueList,IndexList> const& plist)
{
   ValueList list = plist.first;
   IndexList indx = plist.second;

   for (ValueList::const_iterator it = list.begin();
                                   it != list.end();
                                 ++it)
   {
      os << it->first << " : " << indx[it->first] << " : ";
      os << /* scientific << */ it->second << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, ValueList const& list)
{
   for (ValueList::const_iterator it = list.begin();
                                   it != list.end();
                                 ++it)
   {
      os << it->first << " : ";
      os << /* scientific << */ it->second << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<StringList,IndexList> const& plist)
{
   StringList list = plist.first;
   IndexList  indx = plist.second;

   for (StringList::const_iterator it = list.begin();
                                   it != list.end();
                                 ++it)
   {
      os << it->first << " : "
         << indx[it->first] << " : " 
         << it->second << endl;

      // os << it->first << " : "
      //    << indx[it->first] << " : "
      //    << split(it->second) << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, StringList const& list)
{
   for (StringList::const_iterator it = list.begin();
                                   it != list.end();
                                 ++it)
   {
      os << it->first << " : " << it->second << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, 
                     pair<FunctionDefList,IndexList> const& plist)
{
   FunctionDefList list = plist.first;
   IndexList       indx = plist.second;

   for (FunctionDefList::const_iterator itf = list.begin();
                                        itf != list.end();
                                      ++itf)
   {
      vector<string> vec = itf->second;

      os << itf->first;

      // for (unsigned j = 0; j < vec.size()-1; ++j)
      // {
      //    os << " , " << vec[j];
      // } 

      os << " : " << indx[itf->first] << " : ";

      // if ( !vec.empty() )
      // {
      //    os << vec.back();
      // }

      os << vec;

      // os << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, FunctionDefList const& list)
{
   for (FunctionDefList::const_iterator itf = list.begin();
                                        itf != list.end();
                                      ++itf)
   {
      vector<string> vec = itf->second;

      os << itf->first;

      // for (unsigned j = 0; j < vec.size()-1; ++j)
      // {
      //    os << " , " << vec[j];
      // } 

      os << " : ";

      // if ( !vec.empty() )
      // {
      //    os << vec.back();
      // }

      os << vec;

      // os << endl;
   }

   return os;
}

//-------------------------------------------------------------------------

ostream& operator<< (ostream& os, vector<string> const& vec)
{
   for (vector<string>::const_iterator it = vec.begin();
                                       it != vec.end();
                                     ++it)
   {
      os << "|" << *it; 
   }

   os << "|";

  return os;
}

//=========================================================================

vector<string> 
split(string const& str, string const& delim)
{
    vector<string> result;

    if ( str.empty() )
    {
        return result;
    }

    string::const_iterator sBeg = str.begin();
    string::const_iterator sEnd = str.end();
    string::const_iterator it = sBeg;

    do 
    {
        while (delim.find(*it) == string::npos && it != sEnd)
        {
            ++it;         // find the position of the first delimiter in str
        }

        string token = string(sBeg, it);                   // grab the token

        if ( !token.empty() ) 
        {
                        // empty token only when str starts with a delimiter
            result.push_back(token); // push the token into a vector<string>
        }

        while (delim.find(*it) != string::npos && it != sEnd)
        {
            ++it;            // ignore the additional consecutive delimiters
        }

        sBeg = it;                           // process the remaining tokens
    } 
    while (it != sEnd);

    return result;
}

//-------------------------------------------------------------------------

void replaceAll(string& str, string const& from, string const& to)
{
    if ( from.empty() )
    {
        return;
    }

    size_t start_pos = 0;

    while( (start_pos = str.find(from, start_pos)) != string::npos )
    {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();    // In case 'to' contains 'from',
                                     // like replacing 'x' with 'yx'
    }
}

//-------------------------------------------------------------------------

size_t replaceSingle(string& str, 
                     string const& from, string const& to, 
                     size_t start_pos)
{
    if ( from.empty() )
    {
        return start_pos;
    }

    if( (start_pos = str.find(from, start_pos)) != string::npos )
    {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }

    return start_pos;
}

//=========================================================================

