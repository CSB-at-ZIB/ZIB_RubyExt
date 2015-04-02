/**
 * Copyright (C) 2015 Zuse Institute Berlin
 *
 * @file    SbmlJacobian.cpp
 *
 * @brief   Implementation of the data structure, 
 *          as declared in SbmlJacobian.h
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    11.03.2015
 *
 */

#include <cstdlib>
#include "Parser.h"
#include "Expression.h"
#include "SbmlJacobian.h"
#include "SbmlModel.h"

/// using namespace std;

//=========================================================================

SbmlJacobian::SbmlJacobian() 
{
   // cerr << "c'tor SbmlJacobian " << this << endl;

   _druledy.clear();
   _druledp.clear();

   _dreacdy.clear();
   _dreacdp.clear();

   _dfdy.clear();
   _dfdp.clear();

   _vdruledy.clear();
   _vdruledp.clear();

   _vdreacdy.clear();
   _vdreacdp.clear();

   _vdfdy.clear();
   _vdfdp.clear();
}

/*
SbmlJacobian::SbmlJacobian(Model const* m) 
{
}  
*/

//-------------------------------------------------------------------------

SbmlJacobian::~SbmlJacobian()
{
   // cerr << "d'tor SbmlJacobian " << this << endl;
}

//=========================================================================

void
SbmlJacobian::set_druledy(StringList&      rule, 
                          IndexList const& irule,
                          IndexList const& ispec)
{
   PARKIN::Parser p;
   ostringstream  outs;
   int            j = 0;
   int            k = 0;

   _druledy.clear();
   _vdruledy.clear();

   for (IndexList::const_iterator it1 = irule.begin();
                                  it1 != irule.end();
                                ++it1)
   {
      string rId = it1->first;
      string rul = rule[rId];

      PARKIN::Expression ast = p.parse(rul);

      ++j; k = 0;

      for (IndexList::const_iterator it2 = ispec.begin();
                                     it2 != ispec.end();
                                   ++it2)
      {
         string sId = it2->first;

         if (rul.find(sId) == string::npos) 
         {
             continue;
         }
         
         ++k;

         outs.str( string() );
         outs.clear();
   
         outs << ast.df(sId);

         // char buf[1024];
         // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", j, k, rul.c_str());
         // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
         //              j, k, outs.str().c_str());

         _druledy[rId][sId] = outs.str(); // buf;

         _vdruledy[UnsignedPair(it1->second,it2->second)] =
                                                   StringPair(rId,sId);
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_druledy() const
{
   return _druledy;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdruledy() const
{
   return _vdruledy;
}

//-------------------------------------------------------------------------

void
SbmlJacobian::set_druledp(StringList&      rule,
                          IndexList const& irule,
                          IndexList const& iparm)
{
   PARKIN::Parser p;
   ostringstream  outs;
   int            j = 0;
   int            k = 0;

   _druledp.clear();
   _vdruledp.clear();

   for (IndexList::const_iterator it1 = irule.begin();
                                  it1 != irule.end();
                                ++it1)
   {
      string rId = it1->first;
      string rul = rule[rId];

      PARKIN::Expression ast = p.parse(rul);

      ++j; k = 0;

      vector<string> tok = split(rul);

      for (unsigned n = 0; n < tok.size(); ++n)
      {
         string                    pId = get_pId(tok[n],rId,iparm);
         IndexList::const_iterator it2 = iparm.find(pId);

         if ( it2 != iparm.end() )
         {
            ++k;
   
            outs.str( string() );
            outs.clear();
   
            outs << ast.df(tok[n]);

            // char buf[2048];
            // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
            //              j, k, rul.c_str());
            // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
            //              j, k, outs.str().c_str());

            _druledp[rId][pId] = outs.str(); // buf;

            _vdruledp[UnsignedPair(it1->second,it2->second)] =
                                                   StringPair(rId,pId);
         }
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_druledp() const
{
   return _druledp;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdruledp() const
{
   return _vdruledp;
}

//=========================================================================

void
SbmlJacobian::set_dreacdy(StringList&      reac, 
                          IndexList const& ireac,
                          IndexList const& ispec)
{
   PARKIN::Parser p;
   ostringstream  outs;
   int            j = 0;
   int            k = 0;

   _dreacdy.clear();
   _vdreacdy.clear();

   for (IndexList::const_iterator it1 = ireac.begin();
                                  it1 != ireac.end();
                                ++it1)
   {
      string rId = it1->first;
      string rea = reac[rId];

      PARKIN::Expression ast = p.parse(rea);

      ++j; k = 0;

      for (IndexList::const_iterator it2 = ispec.begin();
                                     it2 != ispec.end();
                                   ++it2)
      {
         string sId = it2->first;

         if (rea.find(sId) == string::npos) 
         {
             continue;
         }
         
         ++k;

         outs.str( string() );
         outs.clear();
   
         outs << ast.df(sId);

         // char buf[1024];
         // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
         //              j, k, rea.c_str());
         // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
         //              j, k, outs.str().c_str());

         _dreacdy[rId][sId] = outs.str(); // buf;

         _vdreacdy[UnsignedPair(it1->second,it2->second)] =
                                                   StringPair(rId,sId);
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_dreacdy() const
{
   return _dreacdy;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdreacdy() const
{
   return _vdreacdy;
}

//-------------------------------------------------------------------------

void
SbmlJacobian::set_dreacdp(StringList&      reac, 
                          IndexList const& ireac,
                          IndexList const& iparm)
{
   PARKIN::Parser p;
   ostringstream  outs;
   int            j = 0;
   int            k = 0;

   _dreacdp.clear();
   _vdreacdp.clear();

   for (IndexList::const_iterator it1 = ireac.begin();
                                  it1 != ireac.end();
                                ++it1)
   {
      string rId = it1->first;
      string rea = reac[rId];

      PARKIN::Expression ast = p.parse(rea);

      ++j; k = 0;

      vector<string> tok = split(rea);

      for (unsigned n = 0; n < tok.size(); ++n)
      {
         string                    pId = get_pId(tok[n],rId,iparm);
         IndexList::const_iterator it2 = iparm.find(pId);

         if ( it2 != iparm.end() )
         {
            ++k;
   
            outs.str( string() );
            outs.clear();
   
            outs << ast.df(tok[n]);

            // char buf[2048];
            // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
            //              j, k, rea.c_str());
            // sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
            //              j, k, outs.str().c_str());

            _dreacdp[rId][pId] = outs.str(); // buf;

            _vdreacdp[UnsignedPair(it1->second,it2->second)] =
                                                   StringPair(rId,pId);
         }
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_dreacdp() const
{
   return _dreacdp;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdreacdp() const
{
   return _vdreacdp;
}

//=========================================================================

/*
void
SbmlJacobian::set_dfdy(StringList const& rate)
{
   _dfdy.clear();
   _vdfdy.clear();
}
*/

void
SbmlJacobian::set_dfdy(StringList&      rate, 
                       IndexList const& ispec)
{
   int j = 0;
   int k = 0;

   _dfdy.clear();
   _vdfdy.clear();

   for (IndexList::const_iterator it1 = ispec.begin();
                                  it1 != ispec.end();
                                ++it1)
   {
      ++j; k = 0;

      for (IndexList::const_iterator it2 = ispec.begin();
                                     it2 != ispec.end();
                                   ++it2)
      {
         char buf[1024];
         // string key = it1->first + "," + it2->first;

         ++k;
   
         sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
                      j, k, rate[it1->first].c_str());

         _dfdy[it1->first][it2->first] = buf;

         _vdfdy[UnsignedPair(it1->second,it2->second)] = 
                                     StringPair(it1->first,it2->first);
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_dfdy() const
{
   return _dfdy;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdfdy() const
{
   return _vdfdy;
}

//-------------------------------------------------------------------------

void
SbmlJacobian::set_dfdp(StringList&      rate, 
                       IndexList const& ispec,
                       IndexList const& iparm)
{
   int j = 0;
   int k = 0;

   _dfdp.clear();
   _vdfdp.clear();

   for (IndexList::const_iterator it1 = ispec.begin();
                                  it1 != ispec.end();
                                ++it1)
   {
      ++j; k = 0;

      for (IndexList::const_iterator it2 = iparm.begin();
                                     it2 != iparm.end();
                                   ++it2)
      {
         char buf[1024];
         // string key = it1->first + "," + it2->first;

         ++k;
   
         sprintf(buf, "\n*** j = %3d ,  k = %3d\n*** %s", 
                      j, k, rate[it1->first].c_str());

         _dfdp[it1->first][it2->first] = buf;

         _vdfdp[UnsignedPair(it1->second,it2->second)] = 
                                     StringPair(it1->first,it2->first);
      }
   }
}

//-------------------------------------------------------------------------

StringArray const&
SbmlJacobian::get_dfdp() const
{
   return _dfdp;
}

//-------------------------------------------------------------------------

ArrayIndex const&
SbmlJacobian::get_vdfdp() const
{
   return _vdfdp;
}

//-------------------------------------------------------------------------

string
SbmlJacobian::get_pId(string const& tok, 
                      string const& id,
                      IndexList const& iparm) const
{
   string pId = tok;
   string  ww = "global_" + tok;
   string  rw = id + "_" + tok;
   string  lw = id + "_local_" + tok;

   if ( iparm.count(lw) > 0 )
   {
      pId = lw;
   } 
   else if ( iparm.count(rw) > 0 )
   {
      pId = rw;
   }
   else if ( iparm.count(ww) > 0 )
   {
      pId = ww;
   }

   return pId;
}

//=========================================================================

//=========================================================================

ostream& operator<< (ostream& os, 
                     pair<StringArray,ArrayIndex> const& parray)
{
   StringArray ary = parray.first;
   ArrayIndex  indx = parray.second;

   for (ArrayIndex::const_iterator it = indx.begin();
                                   it != indx.end();
                                 ++it)
   {
       os << (it->first).first << "," << (it->first).second << " : ";
       os << (it->second).first << "," << (it->second).second << " : ";
       os << ary[(it->second).first][(it->second).second] << endl;
   }

   return os;
}

//-------------------------------------------------------------------------


