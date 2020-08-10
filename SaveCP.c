/***************************************************************************
                          SaveCP.c  -  description
                             -------------------
    begin                : Mon Feb 2 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/

// 20120512: (HY) fix a printf specifier bug for size_t,
//                replacing %8u with %lu

#include <stdlib.h>
#include <stdio.h>
#include "changepoint.h"

void SaveCP(fp,p)
struct changepoint  *p;
FILE                *fp;
{
  if (p != NULL) {
    SaveCP(fp,p->left);
    fprintf( fp, "%lu %lu %lu\n", (long unsigned int) p->i, (long unsigned int) p->il, (long unsigned int) p->ir );
    SaveCP(fp,p->right);
    }
  }
