/**
 * Copyright (C) 2015 Zuse Institute Berlin
 *
 * @file    sbml2adolc.cpp
 *
 * @brief   Translate an SBML file/model to C/C++ code representing 
 *          the underlying ODE structure: Rule, Reaction, and Event 
 *          formulas in the given SBML Document
 *
 *          Compile with (see `Makefile´ !!)
 *             'g++ -O3 -c -o sbml2adolc.o sbml2adolc.cpp'
 *             'g++ -O3 -c -o SbmlModel.o SbmlModel.cpp'
 *             'g++ -o sbml2adolc sbml2adolc.o SbmlModel.o -lsbml -ladolc' 
 *          (only if libsbml *and* libadolc is installed, of course!)
 *
 * @author  Thomas Dierkes (dierkes at zib dot de)
 * 
 * @date    30.06.2015
 *
 */

#include <string>
#include <iostream>
#include <iterator>
#include "SbmlModel.h"


int
main (int argc, char* argv[])
{
  SBMLDocument const* sbmldoc;
  string              hopt = "-h";
  string              txtopt = "-txt";
  string              prgpath = argv[0];
  string              filename;
  long unsigned       pos = prgpath.rfind("/");
  bool                adolc = true;

  if ( pos != string::npos )
  {
    prgpath.erase(pos); 
  }
  else
  {
    prgpath = ".";
  }


  if (argc == 1)
  {
    cin >> noskipws;

    istream_iterator<char> it(cin);
    istream_iterator<char> end;
    string                 xml(it,end);

    sbmldoc = readSBMLFromString(xml.c_str());
  }
  else if (    (argc == 2) 
            && (hopt.compare(argv[1]) == 0 || txtopt.compare(argv[1]) == 0)
          )
  {
    cerr << endl << "Usage: ./sbml2adolc [-h|-txt] filename{.xml|.sbml| }" << endl;
    cerr << endl;
    return 1;
  }
  else if ( (argc > 2) && (txtopt.compare(argv[1]) == 0) )
  {
    adolc    = false;
    filename = argv[2];
    sbmldoc  = readSBML(filename.c_str());
  }
  else 
  {
    filename = argv[1];
    sbmldoc  = readSBML(filename.c_str());
  }

  if (sbmldoc != 0)
  {
    if (sbmldoc->getNumErrors() > 0)
    {
      cerr << "Encountered the following SBML errors:" << endl;
      sbmldoc->printErrors(cerr);
      delete sbmldoc;
      return -1;
    }

    Model const* m = sbmldoc->getModel();

    if (m == 0)
    {
      cerr << "No model present." << endl;
      delete sbmldoc;
      return -2;
    }


    pos = string::npos;
    pos = filename.rfind("/");

    if ( pos != string::npos )
    {
       filename.erase(filename.begin(),filename.begin()+pos+1);
    }

    SbmlModel model(m, filename, adolc);


    if ( adolc )
    {
      cout << model.toAdolC( prgpath + "/YDOT_LIMEXcpp_TEMPLATE.txt" ) << endl;
    }
    else
    {
      cout << model << endl;
    }

    delete sbmldoc;
  }
  else
  {
    cerr << "Cannot read SBML stream/file." << endl;
    return -3;
  }

  return 0;
}
