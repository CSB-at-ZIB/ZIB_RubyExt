#include <iostream>

#define STRLEN 11

extern "C" {
  void get_str_array_( char (*str)[STRLEN], int* n, int _len = STRLEN ); 
}

using namespace std;


int main()
{
  int    n = 9;
  // char astr[][STRLEN] = { "can you", "test", "me", "please" };
  char   astr[n][STRLEN]; // = { };
  double pi = 3.1415;

  cout << " sizeof(astr) = " << sizeof(astr) << endl;


  get_str_array_( astr, &n );


  cout << " pi ~ " << pi << ",  n = " << n << endl;


  for (int j = 0; j < n; ++j)
  {
     astr[j][STRLEN-1] = '\0';
     cout << "astr[" << j << "] = '" << astr[j] << "'" << endl;
     /*
     cout << "astr[" << j << "] = '";
     for (int k = 0; k < STRLEN; ++k) 
     {
        cout << astr[j][k];
     }
     cout << "'" << endl;
     */
  }


  return 0;
}
