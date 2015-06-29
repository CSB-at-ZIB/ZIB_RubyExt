#include <iostream>

extern "C"
{
  #include "ydot_LIMEX.h"
}

using namespace std;

int main()
{
   char    buf[256][MAXIDSTRLEN];
   double*  ptr = &sbmlvariables_.start;
   int     nCom = 0;
   int     nSpe = 0;
   int     nPar = 0;
   double  c[1] = { };
   double  s[1] = { };
   double  p[1] = { };
   int        k = 0;

   init_ode_( c,&k,&nCom, s,&k,&nSpe, p,&k,&nPar );
 
   cout << "nCom = " << nCom << endl;
   for (int j = 0; j < nCom; ++j)
   {
      cout << '\t' << *ptr++;
   }
   cout << endl;

   get_species_ids_( buf, &k );

   cout << "nSpe = " << nSpe << endl;
   for (int j = 0; j < nSpe; ++j)
   {
      cout << '\t' << *ptr++ << " (" << buf[j] << ')';
   }
   cout << endl;

   get_parameter_ids_( buf, &k );

   cout << "nPar = " << nPar << endl;
   for (int j = 0; j < nPar; ++j)
   {
      cout << '\t' << *ptr++ << '\t' << '(' << buf[j] << ')' << endl;
   }

   return 0; 
}
