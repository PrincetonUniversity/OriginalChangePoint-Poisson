/***************************************************************************
                          CheckCP.c  -  description
                             -------------------
    begin                : Wed Jan 28 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/

// References:
// 1.
// 2. "A problem with the likelihood ratio test for a change-point hazard rate model,"
//    Robin Henderson, Biometrika, 77, 835-843 (1990).
//
// Change Log:
// 20180220: (HY) Fixed a bug on line 83.
// 20120512: (HY) Implement a bug fix at line ca. 83, credited to Arno van Amersfoort
//                fixed a printf specifier bug for size_t variables
//                replacing %i with %lu
// 20040326: (HY) Implement Henderson's static using exact distribution.
//                As of now, this subroutine still cannot cope with data length longer
//                than critical value length
// 20040206: (HY) implement modifications by Robin Henderson Ref. [2]
// 20040203: (HY) need to add a clause to remove points whose neighboring
//           change points are within its uncertainty.

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

//printf("%i\n",istart);

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
        //iend=((i-10)>1 ? i : 10); // This is a buggy statement
        iend = (i > 10) ? i : 10; // bug fix credited to Arno van Amersfoort
//printf("%i ",iend);
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

/*

  while ( Ncp > 0 && deletion_adjustment ) {
    for (j=1;j<=Ncp;j++) {
      n = cps[j+1] - cps[j-1];
      if ( n < 1 ) {
        // remove this spurious point at index j
        // increase the number of deleted point by 1
        (*Ncpdlt)++;
        // now, delete the cp node
        if (trace==1) printf( "  -=> deleting adjacent transition at index %i\n", cps[j] );
        DeleteCPNode(cps[j],cp_root,&h);
        // The next step may take a lot of time, using a bi-directional chain
        // data structure may help to speed up the process.
        // move up the working array
        for (i=j;i<=Ncp;i++) {
          cps[i] = cps[i+1];
          cpsl[i] = cpsl[i+1];
          cpsr[i] = cpsr[i+1];
          }
        Ncp--;
        deletion_adjustment=true;
        break;
        }
      else{ // check if this point is good
        // create work space
        Tk = (double *) malloc( (n+1) * sizeof(double) );
        lambda1 = (double *) malloc( (n+1) * sizeof(double) );
        lambda2 = (double *) malloc( (n+1) * sizeof(double) );
        llrt = (double *) malloc( (n+1) * sizeof(double) );
        u = (double *) malloc( n * sizeof(double) );
        u1 = (double *) malloc( n * sizeof(double) );
        v2 = (double *) malloc( n * sizeof(double) );
        v21 = (double *) malloc( n * sizeof(double) );
        sk = (double *) malloc( n * sizeof(double) );
        // total time in this trajectory segment
        T = traj[cps[j+1]].time - traj[cps[j-1]].time;
        lambda0 = (double) n / T + DBL_EPSILON;
        // calculate the critical region
        //critical_region = C_ac(alpha,n);
        critical_region = ta[n]; // use tabulated *exact* values
        // compute expection values and standard deviations
        for (i=1; i<n; i++)
          u[i] = u1[i] = v2[i] = v21[i] = 0.0;
        for (i=1; i<n; i++) {
          iL = (double) i;
          for (k=i; k<n; k++) {
            kL = (double) k;
            u[i] -= 1.0/kL;                         // expection value for log(Vk)
            v2[i] += 1.0/kL/kL;                     // variance for [log(Vk)]^2
            }
          }
        for (i=1; i<n; i++) {
          u1[n-i] = u[i];                           // expection value for log(1-Vk)
          v21[n-i] = v2[i];                         // variance for [log(1-Vk)]^2
          }
        nL = (double) n;
        term = M_PI*M_PI/6.0-v2[1];
        for (i=1; i<n; i++) {
          iL = (double) i;
          sk[i] = 4.0*iL*iL*v2[i]+4.0*(nL-iL)*(nL-iL)*v21[i]-8.0*iL*(nL-iL)*term;
          sk[i] = sqrt(sk[i]);
          }
        for ( k=1,k_max=0,llrt_max=0.0; k<n; k++ ) {
          kL = (double) k;
          Tk[k] = traj[cps[j-1]+k].time - traj[cps[j-1]].time;
          // Henderson's standardization
          Vk = Tk[k] / T;
          llrt[k] = (-2.0*kL*log(Vk)+2.0*kL*u[k]-2.0*(nL-kL)*log(1.0-Vk)
                    +2.0*(nL-kL)*u1[k])/sk[k]+0.5*log(4.0*kL*(nL-kL)/nL/nL);
          // Generalized Likelihood Ratio Test

          lambda1[k] = k / Tk[k] + DBL_EPSILON;
          lambda2[k] = (double) (n-k) / (T-Tk[k]) + DBL_EPSILON;
          llrt[k] = (double) k*log(lambda1[k]) + (double) (n-k)*log(lambda2[k])
                  - (double) n*log(lambda0);
          llrt[k] = 2.0 * llrt[k];
                  //+ log((double)4*k*(n-k)/n/n); // this weighting term follows Ref[2]

         if ( llrt[k] > llrt_max ) {
            // find the maximum of the likelihood functions
            llrt_max = llrt[k];
            k_max = k;
            }
          }

        if ( llrt_max > critical_region ) {
          // Is it the same as cps[j] ?
          //if ( (k_max+cps[j-1] <= cpsr[j]) && (k_max+cps[j-1] >= cpsl[j]) )
          if ( k_max+cps[j-1] == cps[j] )
            // within the confidence interval
            deletion_adjustment = false;
          else { // Not the same as cps[j], adjust it
            if (trace==1) printf( "] %i ==> [%u] (%u) [%u] in [%u .. %u]\n",
                                   j, cpsl[j], k_max+cps[j-1], cpsr[j], cps[j-1], cps[j+1] );
            // max{llrt} is the refined change point
            // find the confidence interval, threshold given by ca
            // find the left-hand bound
            k = k_max;
            printf( "k_max = %u, llrt[k]+ca[n]-llrt_max = %le\n", k_max, llrt[k]+ca[n]-llrt_max );
            printf( "ca(n) = %le\n", ca[n] );
            while ( (k>0) && (llrt[k]+ca[n]-llrt_max > 0) ) k--;
            LB = k;
            printf( "LB = %u, n=%u\n", LB, n );
            // find the right-hand bound
            k = k_max;
            while ( (k<n) && (llrt[k]+ca[n]-llrt_max > 0) ) k++;
            RB = k;

            if (trace==1) printf( "  ~=> adjusting cps[%i] = %i to %i,\n", j, cps[j], k_max+cps[j-1] );
            h = false;
            if (trace==1) printf( "   -> deleting cps[%i] = %i\n", j, cps[j] );
            DeleteCPNode(cps[j],cp_root,&h);
            h = false;
            if (trace==1) printf( "   +> adding cps[%i] = %i\n", j, k_max+cps[j-1] );
            if (trace==1) printf( "   o> [%u] (%u) [%u]\n", LB+cps[j-1], k_max+cps[j-1], RB+cps[j-1] );
            AddCPNode( (size_t) k_max+cps[j-1], LB+cps[j-1], RB+cps[j-1], cp_root, &h );
            cps[j] = (size_t) k_max+cps[j-1];
            cpsl[j] = (size_t) LB+cps[j-1];
            cpsr[j] = (size_t) RB+cps[j-1];
            deletion_adjustment = true;
            break;
            }
          }
        else {
          // remove this spurious point at index j
          // increase the number of deleted point by 1
          if (trace==1) printf( "  X=> removing cps[%i] = %i\n", j, cps[j] );
          (*Ncpdlt)++;
          DeleteCPNode(cps[j],cp_root,&h); // now, delete the cp node
          // The next step may take a lot of time, using a bi-directional chain
          // data structure may help to speed up the process.
          for (i=j;i<=Ncp;i++) {
            cps[i] = cps[i+1]; // move up the working array
            cpsl[i] = cpsl[i+1];
            cpsr[i] = cpsr[i+1];
            }
          Ncp--;
          deletion_adjustment=true;
          break;
          }
        // free up work space
        free(Tk);
        free(lambda1);
        free(lambda2);
        free(llrt);
        free(u);
        free(u1);
        free(v2);
        free(v21);
        free(sk);
        } // end of if n > 1
      } // end of for j
    } // end of while
  free(cps);
  free(cpsl);
  free(cpsr);



  }

*/
