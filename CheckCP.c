/***************************************************************************
                          CheckCP.c  -  description
                             -------------------
    begin coding                : Tue Jan 27 14:25:30 PST 2004
    peer-reviewed publication   : J. Phys. Chem. B, 109, 617-628 (2005)
    code initial public release : 2020
    copyright                   : © Haw Yang 2020
    email                       : hawyang@princeton.edu
***************************************************************************/

// Change log (for public release):
// 20200811: (HY) Start code clean up for v2.00 initial public release.

// References:
// 1. "A problem with the likelihood ratio test for a change-point hazard rate model,"
//    Robin Henderson, Biometrika, 77, 835-843 (1990).
//

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>
#include "changepoint.h"

void CheckCP(cp_root, traj, alpha, beta, N, Ncpdlt, trace, ta, ca, Na )
struct changepoint  **cp_root;
struct data         *traj;
double              alpha, beta;
double              *ta, *ca;
int                 *Ncpdlt;
size_t              N, Na;
int                 trace;
{
  size_t            *cps=NULL, *cpsl=NULL, *cpsr=NULL, cp1;
  int               Ncp=0, Ny=0, Ny2=0;
  size_t            i, istart=1, iend=0;
  enum bool         deletion_adjustment=true, h=false;

  // create an array, cps[], to store data
  *Ncpdlt=0;
  while(deletion_adjustment){
    deletion_adjustment=false;
    MakeCPArray( *cp_root, &cps, &cpsl, &cpsr, &Ncp );
    cps = (size_t *) realloc( cps, (Ncp+2)*sizeof(size_t) );
    cps[0] = 0;
    cps[Ncp+1] = N;

    cpsl = (size_t *) realloc( cpsl, (Ncp+2)*sizeof(size_t) );
    cpsl[0] = 0;
    cpsl[Ncp+1] = N;

    cpsr = (size_t *) realloc( cpsr, (Ncp+2)*sizeof(size_t) );
    cpsr[0] = 0;
    cpsr[Ncp+1] = N;

    if(iend>1)
      istart=iend;
    else
      istart=1;

    for(i=istart;i<Ncp+1;i++){
      DeleteCPNode(cps[i],cp_root,&h);
      Ncp--;
      if(cps[i]<cpsr[i-1]){
        deletion_adjustment=true;
        (*Ncpdlt)++;
        break;
        }
      cp1=FindCP(cp_root,traj,cps[i-1]+1,cps[i+1]-1,alpha,beta,&Ncp,ta,ca,Na,3);
      if (cp1==0){
        printf("    %lu/%i Checking node %lu %lu %lu...",(long unsigned int) i,Ncp,(long unsigned int) cpsr[i-1]+1, (long unsigned int) cps[i], (long unsigned int) cpsl[i+1]-1);
        printf("not found! RESTARTING\n");
        deletion_adjustment=true;
        (*Ncpdlt)++;
        iend = (i > 10) ? i : 10; // bug fix credited to Arno van Amersfoort
        break;
        }
      else if(cp1!=cps[i]){
        printf("    %lu/%i Checking node %lu %lu %lu...",(long unsigned int)i,Ncp,(long unsigned int)cpsr[i-1]+1,(long unsigned int)cps[i],(long unsigned int)cpsl[i+1]-1);
        printf("found at %lu.\n",(long unsigned int) cp1);
        deletion_adjustment=true;
        iend=i;
        break;
        }
      }
    free(cps);
    free(cpsl);
    free(cpsr);
    cps=NULL;
    cpsl=NULL;
    cpsr=NULL;
    Ncp=0;
    }
}

