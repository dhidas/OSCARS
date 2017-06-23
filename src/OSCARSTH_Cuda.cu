////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 11:33:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include <cuComplex.h>

#include "OSCARSTH_Cuda.h"

#include "OSCARSTH.h"

#include <cmath>
#include <fstream>
#include <sstream>



#define NTHREADS_PER_BLOCK 512






extern "C" int OSCARSTH_Cuda_GetDeviceCount ()
{
  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);

  return ngpu;
}





