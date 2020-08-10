/***************************************************************************
                          util.c  -  description
                             -------------------
    begin                : Sat Jan 31 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/
// 20080121: (HY) make the code more compact.
#include <stdlib.h>
 
double **dmatrix( nrows, ncolumns )
int nrows, ncolumns;
{
  int i;
  double **m;

  m = malloc( (nrows+1) * sizeof(double *));
  for(i = 0; i <= nrows; i++) 
    m[i] = malloc( (ncolumns+1) * sizeof(double));

  return m;
}

void free_dmatrix( m, nrows, ncolumns )
double **m;
int nrows, ncolumns;
{
  int i;
  for(i = 0; i <= nrows; i++)
    free((void *) m[i]);
  free((void *) m);
}

int **imatrix( nrows, ncolumns )
int nrows, ncolumns;
{
  int i;
  int    **m;

  m = malloc( (nrows+1) * sizeof(int *));
  for(i = 0; i <= nrows; i++) 
    m[i] = malloc( (ncolumns+1) * sizeof(int));

  return m;
}

void free_imatrix( m, nrows, ncolumns )
int    **m;
int nrows, ncolumns;
{
  int i;
  for(i = 0; i <= nrows; i++)
    free((void *) m[i]);
  free((void *) m);
}

