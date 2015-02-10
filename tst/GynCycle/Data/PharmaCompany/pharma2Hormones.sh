#! /bin/bash

scriptpath="`dirname $0`"
dataf="pharma_tabelle2.csv"

if [ ! -f "${scriptpath}/${dataf}" ]; then
  echo "`basename ${0}`: ERROR : data file '${dataf}' missing."
  echo
  exit -1 
fi

daymin=`sort -n -t',' -k7 "${scriptpath}/${dataf}" | head -1 | cut -d',' -f7`
daymax=`sort -r -n -t',' -k7 "${scriptpath}/${dataf}" | head -1 | cut -d',' -f7`

awk -F"," -v dmin=${daymin} -v dmax=${daymax} '
          BEGIN {
                    lhfile="Hormones-LH.csv"
                   fshfile="Hormones-FSH.csv"
                    p4file="Hormones-P4.csv"
                  progfile="Hormones-Progesterone.csv"
                    e2file="Hormones-E2.csv"
                  estrfile="Hormones-Estradiol.csv"
                }
           /LH/ {  lh[$1] =  lh[$1] "," $3 ;  lhtp[$1] =  lhtp[$1] "," $7 ; }
          /FSH/ { fsh[$1] = fsh[$1] "," $3 ; fshtp[$1] = fshtp[$1] "," $7 ; }
 /Progesterone/ {  p4[$1] =  p4[$1] "," $3 ;  p4tp[$1] =  p4tp[$1] "," $7 ;
                 prog[$1] =  prog[$1] "," $3 * 1000.0/314.47 ; # P4: 314.47 [g/mol]  
               progtp[$1] =  progtp[$1] "," $7 ; }
    /Estradiol/ {  e2[$1] =  e2[$1] "," $3 ;  e2tp[$1] =  e2tp[$1] "," $7 ;  
                 estr[$1] =  estr[$1] "," $3 * 1000.0/272.39 ;  # E2: 272.39 [g/mol] 
               estrtp[$1] =  estrtp[$1] "," $7 ; }
            END {
                  headline = "Cas,BMI," ;
                  for (j=dmin; j<=dmax; ++j) headline = headline j "," ;

                  print headline >lhfile;
                  process(lh, lhtp, dmin, dmax, lhfile);

                  print headline >fshfile;
                  process(fsh, fshtp, dmin, dmax, fshfile);

                  print headline >p4file;
                  process(p4, p4tp, dmin, dmax, p4file);

                  print headline >progfile;
                  process(prog, progtp, dmin, dmax, progfile);

                  print headline >e2file;
                  process(e2, e2tp, dmin, dmax, e2file);

                  print headline >estrfile;
                  process(estr, estrtp, dmin, dmax, estrfile);

                  # for (x in lh)
                  # {
                  #   split(lh[x],meas) ;
                  #   split(lhtp[x],tp) ;
                  #   k = 2 ; line = x ",9.9," ;
                  #   for (j=dmin; j<=dmax; ++j)
                  #   {
                  #      if ( j==tp[k] )
                  #      {
                  #        line = line meas[k] "," ;
                  #        ++k ;
                  #      }
                  #      else
                  #      {
                  #        line = line "," ;
                  #      }
                  #   }
                  #   print line >>lhfile;
                  # }
                }

function process(data,tpoint,dmin,dmax,outf,  k,line,j)
                {
                  for (x in data)
                  {
                    if ( x==2 || x==7 ) continue ;

                    split(data[x],meas) ;
                    split(tpoint[x],tp) ;
                    k = 2 ; line = x ",9.9," ;
                    for (j=dmin; j<=dmax; ++j)
                    {
                       if ( j==tp[k] )
                       {
                         line = line meas[k] "," ;
                         ++k ;
                       }
                       else
                       {
                         line = line "," ;
                       }
                    }
                    print line >>outf;
                  }
                }
' "${scriptpath}/${dataf}"
