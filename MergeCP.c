/***************************************************************************
                          MergeCP.c  -  description
                             -------------------
    begin                : Mon Feb 2 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <gsl/gsl_sf_gamma.h>
#include "changepoint.h"

void MergeCP( Yi, Ny, Yem, Ng_max, Yim, Nym, Yemm )
double              **Yi;                       // pointer to grouping matrix
       // Yi[:][1], intensity levels between change points
       // Yi[:][2:Ng_max], the group to which it belongs
       // Yi[:][Ng_max+1], number of photons in this group
       // Yi[:][Ng_max+2], time duration of this group
double              **Yim;
int                 Ny;
int                 *Nym;
struct em_group     Yem, Yemm;                  // EM cluster results
int                 Ng_max;                     // maximum number of groups
{
  int               G,i,k,old_k;

  Nym[1] = Ny;
  for (G=2; G<=Ng_max; G++) {
    for (i=1,old_k=0,Nym[G]=0; i<=Ny; i++) {
      k = Yem.c[i][G];
      if ( k != old_k ) {
        Nym[G]++;
        Yim[Nym[G]][1] = Yem.intensity[i][G];
        Yim[Nym[G]][G] = k;
        Yim[Nym[G]][Ng_max+1] = Yi[i][Ng_max+1];
        Yim[Nym[G]][Ng_max+2] = Yi[i][Ng_max+2];
        Yemm.intensity[Nym[G]][G] = Yem.intensity[i][G];
        Yemm.c[Nym[G]][G] = k;
        Yemm.prob[Nym[G]][G] = Yem.prob[i][G];
        }
      else {
        Yim[Nym[G]][Ng_max+1] += Yi[i][Ng_max+1];
        Yim[Nym[G]][Ng_max+2] += Yi[i][Ng_max+2];
        }
      old_k = k;
      } // end of for i
    } // end of for G
  } // end of MergeCP

