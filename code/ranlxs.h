/*******************************************************************************
*
* Random number generator "ranlxs"
*
* See the notes 
*
*   "User's guide for ranlxs and ranlxd [C programs]" (December 1997)
*
*   "Double precision implementation of the random number 
*    generator ranlux" (December 1997)
*
* for a detailed description
*
* The externally accessible functions are 
*
*   void ranlxs(float r[],int n)
*     Computes the next n single-precision random numbers and 
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
* 
*   void rlxs_init(int level,int seed)
*     Initialization of the generator
*
*   void rlxs_get(int state[])
*     Extracts the current state of the generator and stores the 
*     information in the array state[25]
*
*   void rlxs_reset(int state[])
*     Resets the generator to the state defined by the array state[25]
*
* Version: 2.2
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 15.07.1999
*
*******************************************************************************/

#ifndef RANLXS_H
#define RANLXS_H

#include <complex>


void ranlxs(double r[],int n);
void rlxs_init(int level,int seed);
void rlxs_get(int state[]);
void rlxs_reset(int state[]);
void ranlxs_push_back(std::vector<double>& r, int n);
void ranlxs_fill_vector(std::vector<double>& r);
void uranlxs_fill_vector(std::vector<double>& r);
void uranlxs_fill_complex_vector(std::vector<std::complex<double> >& r);

#endif //#ifndef RANLXS_H
