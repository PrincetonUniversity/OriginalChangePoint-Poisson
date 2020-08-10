/***************************************************************************
                          changepoint.h  -  description
                             -------------------
    begin                : Tue Jan 27 2004
    copyright            : (C) 2004-2018 by Haw Yang
    email                : hawyang@princeton.edu
 ***************************************************************************/

// 20180220: (HY) updated Makefile for recent OSX changes
// 20120512: (HY) revised the Makefile to accomodate new minGW/msys
//                also implemented bux fixes by Arno van Amersfoort
// 20081221: (HY) improve corss platform compilation, refiled as v 1.11
// 20080121: (HY) start develop cross-platform, currently MingW and linux-gnu
// 20040327: (HY) modified struct group to include time uncertainties

#define NA_OVERLAP 200
#define CHANGEPOINT_VERSION "1.13-"
#define COMPILE_DATE "20180220"

#ifdef __MINGW32__
  #ifndef PLATFORM
    #define PLATFORM "MinGW"
  #endif
#endif
#ifdef __APPLE__
  #ifndef PLATFORM
    #define PLATFORM "MacOSX"
  #endif
#endif
#ifdef __GNUC__
  #ifndef PLATFORM
    #define PLATFORM "Linux-GNU"
  #endif
#endif

/*******************************
*         Data Structures      *
*******************************/
enum bool {false=0, true=1};
struct data {
  double        time;       // chronological time stamp
  double        dt;         // inter-photon duration
  double        value;      // numeric value of the data
  };

struct changepoint {
  size_t              i;    // change point index to the trajectory
  size_t              il;   // change point confidence left bound
  size_t              ir;   // change point confidence right bound
  struct changepoint  *left;
  struct changepoint  *right;
  int                 bal;
  };

struct group {
  size_t              n;        // total number of photons in this group
  double              T;        // total time duration of this group
  double              dT;       // uncertainties in time duration
  size_t              m;        // total number of samples
  size_t              *member;  // an array containing the member index of this group
  };

struct em_group {
  int                 **c;          // the classification to which this point belong
  double              **intensity;  // estimated intensity value
  double              **prob;       // posterior probability of this classification
  };

      
/*******************************
*      Function Declaration    *
*******************************/
size_t FindCP();
void CheckCP();
void AddCPNode();
void DeleteCPNode();
void MakeCPArray();
void SizeCPArray();
int FindCPArray();
void AHCluster();
void EMCluster();
void BICCluster();
void SaveCP();
void MergeCP();
double zbrent();
double Vostrikova();
double C_ac();
double **dmatrix();
int **imatrix();
void free_dmatrix();
void free_imatrix();

