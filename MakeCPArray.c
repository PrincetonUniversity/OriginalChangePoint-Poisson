/***************************************************************************
                          MakeCPArray.c  -  description
                             -------------------
    begin                : Wed Jan 28 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/

// 20040326: (HY) modify so that it will also read in confidence bands
 
#include <stdlib.h>
#include <stdio.h>
#include "changepoint.h"
 
void MakeCPArray(p,cps,cpsl,cpsr,Ncp)
struct changepoint  *p;
size_t              **cps,**cpsl,**cpsr;
int                 *Ncp;
{
  if (p != NULL) {
    MakeCPArray(p->left,cps,cpsl,cpsr,Ncp);
    (*Ncp)++;
    (*cps) = realloc(*cps,((*Ncp)+1)*sizeof(size_t));
    (*cpsl) = realloc(*cpsl,((*Ncp)+1)*sizeof(size_t));
    (*cpsr) = realloc(*cpsr,((*Ncp)+1)*sizeof(size_t));
    (*cps)[*Ncp] = p->i;
    (*cpsl)[*Ncp] = p->il;
    (*cpsr)[*Ncp] = p->ir;
    MakeCPArray(p->right,cps,cpsl,cpsr,Ncp);
    }
  }

void SizeCPArray(struct changepoint *p,int *N)
{
  if(p!=NULL){
    SizeCPArray(p->left,N);
    SizeCPArray(p->right,N);
    (*N)++;
    }
}
