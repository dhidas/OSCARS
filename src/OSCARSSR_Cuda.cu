////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue May 24 11:33:19 EDT 2016
//
////////////////////////////////////////////////////////////////////

#include <cuComplex.h>

#include "OSCARSSR_Cuda.h"

#include "OSCARSSR.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "TVector3DC.h"
#include "TSpectrumContainer.h"


#define NTHREADS_PER_BLOCK 512






extern "C" int OSCARSSR_Cuda_GetDeviceCount ()
{
  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);

  return ngpu;
}





std::string OSCARSSR_Cuda_GetDeviceProperties (int const i)
{
  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);

  char buf[300];

  if (i >= ngpu) {
    sprintf(buf, "ERROR: GPU %i Not available", i);
    return std::string(buf);
  }

  std::string ret = "";

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, i);

  sprintf(buf, "  Device name: %s\n", prop.name);
  ret += std::string(buf);
  sprintf(buf, "  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
  ret += std::string(buf);
  sprintf(buf, "  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
  ret += std::string(buf);
  sprintf(buf, "  Peak Memory Bandwidth (GB/s): %f\n", 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
  ret += std::string(buf);

  return ret;
}





__device__ static __inline__ void Orthogonal(double *a, double *b)
{
  // Return a vector which is orthogonal vector a
  double xx = a[0] < 0.0 ? -a[0] : a[0];
  double yy = a[1] < 0.0 ? -a[1] : a[1];
  double zz = a[2] < 0.0 ? -a[2] : a[2];
  if (xx < yy) {
    if (xx < zz) {
      b[0] = 0;
      b[1] = a[2];
      b[2] = -a[1];
    } else {
      b[0] = a[1];
      b[1] = -a[0];
      b[2] = 0;
    }
  } else {
    if (yy < zz) {
      b[0] = -a[2];
      b[1] = 0;
      b[2] = a[0];
    } else {
      b[0] = a[1];
      b[1] = -a[0];
      b[2] = 0;
    }
  }
  return;
}



__host__ __device__ static __inline__ cuDoubleComplex cuCexp(cuDoubleComplex x)
{
  double factor = exp(x.x);
  return make_cuDoubleComplex(factor * cos(x.y), factor * sin(x.y));
}









__global__ void OSCARSSR_Cuda_FluxGPUMulti (double *x, double *y, double *z, double *bx, double *by, double *bz, double *sx, double *sy, double *sz, double *dt, int *nt, int *ns, double *C0, double *C2, double *C, double *Omega, int *ifirst, double *flux)
{
  // Check that this is within the number of spectrum points requested
  int const ith = threadIdx.x + blockIdx.x * blockDim.x;
  int const is = ith + *ifirst;
  if (is >= *ns) {
    return;
  }


  // Complex i
  cuDoubleComplex I = make_cuDoubleComplex(0, 1);

  cuDoubleComplex ICoverOmega = make_cuDoubleComplex(0, (*C) / (*Omega));

  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);


  // Loop over all points in trajectory
  for (int i = 0; i < *nt; ++i) {

    // Distance to observer
    double const D = sqrt( pow( (ox) - x[i], 2) + pow( (oy) - y[i], 2) + pow((oz) - z[i], 2) );

    // Normal in direction of observer
    double const NX = ((ox) - x[i]) / D;
    double const NY = ((oy) - y[i]) / D;
    double const NZ = ((oz) - z[i]) / D;

    // Exponent for fourier transformed field
    cuDoubleComplex Exponent = make_cuDoubleComplex(0, (*Omega) * ((*dt) * i + D / (*C)));

    cuDoubleComplex X1 = make_cuDoubleComplex((bx[i] - NX) / D, -(*C) * NX / ((*Omega) * D * D));
    cuDoubleComplex Y1 = make_cuDoubleComplex((by[i] - NY) / D, -(*C) * NY / ((*Omega) * D * D));
    cuDoubleComplex Z1 = make_cuDoubleComplex((bz[i] - NZ) / D, -(*C) * NZ / ((*Omega) * D * D));

    cuDoubleComplex MyEXP = cuCexp(Exponent);
    //cuDoubleComplex MyEXP = make_cuDoubleComplex( exp(Exponent.x) * cos(Exponent.y), exp(Exponent.x) * sin(Exponent.y));

    cuDoubleComplex X2 = cuCmul(X1, MyEXP);
    cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
    cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


    SumEX = cuCadd(SumEX, X2);
    SumEY = cuCadd(SumEY, Y2);
    SumEZ = cuCadd(SumEZ, Z2);

    // Sum in fourier transformed field (integral)
    //SumEX += (TVector3DC(B) - (N *     (One + (ICoverOmega / (D)))     )) / D * std::exp(Exponent);
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEZ);


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;

  // Multiply field by Constant C1 and time step
  //SumE *= C1 * DeltaT;

  // Set the flux for this frequency / energy point
  //Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);

  flux[ith] = (*C2) * (EX + EY + EZ);

  return;
}








__global__ void OSCARSSR_Cuda_FluxGPU4 (double *x, double *y, double *z, double *bx, double *by, double *bz, double *sx, double *sy, double *sz, double *dt, int *nt, int *ns, double *C0, double *C2, double *C, double *Omega, double *flux)
{
  // Check that this is within the number of spectrum points requested
  int is = threadIdx.x + blockIdx.x * blockDim.x;
  // Max number for shared memory
  int const NSHAREDMAX = 1000;

  // Number for each thread to copy from global to shared memory
  int const NToCopyPerThread = (int) NSHAREDMAX / NTHREADS_PER_BLOCK;

  // Actual number of elements in shared memory to use
  int const NSHARED = NToCopyPerThread * NTHREADS_PER_BLOCK;

  // Define the shared memory
  __shared__ double sh_x[NSHARED];
  __shared__ double sh_y[NSHARED];
  __shared__ double sh_z[NSHARED];
  __shared__ double sh_bx[NSHARED];
  __shared__ double sh_by[NSHARED];
  __shared__ double sh_bz[NSHARED];


  // observer
  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);



  // Total number of copies
  int const NTotalCopies = *nt / NSHARED + 1;

  // Local offset for this thread in shared memory
  int const ThreadOffset = NToCopyPerThread * threadIdx.x;

  // icp is Copy Number
  for (int icp = 0; icp < NTotalCopies; ++icp) {

    // Offset instarting point in global array
    int const GlobalOffset = icp * NSHARED;


    __syncthreads();

    // icpth is the copy number in a thread
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_x[ThisThreadSharedIndex]  = x[GlobalIndex];
    }
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_y[ThisThreadSharedIndex]  = y[GlobalIndex];
    }
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_z[ThisThreadSharedIndex]  = z[GlobalIndex];
    }
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_bx[ThisThreadSharedIndex] = bx[GlobalIndex];
    }
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_by[ThisThreadSharedIndex] = by[GlobalIndex];
    }
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_bz[ThisThreadSharedIndex] = bz[GlobalIndex];
    }


    __syncthreads();

    if (is >= *ns) {
      continue;
    }

    for (int ish = 0; ish < NSHARED; ++ish) {
      int const i = GlobalOffset + ish;

      if (i >= *nt) {
        break;
      }



      // Distance to observer
      double const D = sqrt( pow( (ox) - sh_x[ish], 2) + pow( (oy) - sh_y[ish], 2) + pow((oz) - sh_z[ish], 2) );

      // Normal in direction of observer
      double const NX = ((ox) - sh_x[ish]) / D;
      double const NY = ((oy) - sh_y[ish]) / D;
      double const NZ = ((oz) - sh_z[ish]) / D;

      // Exponent for fourier transformed field
      cuDoubleComplex Exponent = make_cuDoubleComplex(0, (*Omega) * ((*dt) * i + D / (*C)));

      cuDoubleComplex X1 = make_cuDoubleComplex((sh_bx[ish] - NX) / D, -(*C) * NX / ((*Omega) * D * D));
      cuDoubleComplex Y1 = make_cuDoubleComplex((sh_by[ish] - NY) / D, -(*C) * NY / ((*Omega) * D * D));
      cuDoubleComplex Z1 = make_cuDoubleComplex((sh_bz[ish] - NZ) / D, -(*C) * NZ / ((*Omega) * D * D));

      cuDoubleComplex MyEXP = cuCexp(Exponent);

      cuDoubleComplex X2 = cuCmul(X1, MyEXP);
      cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
      cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


      SumEX = cuCadd(SumEX, X2);
      SumEY = cuCadd(SumEY, Y2);
      SumEZ = cuCadd(SumEZ, Z2);








    }
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEZ);


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;


  flux[is] = (*C2) * (EX + EY + EZ);

  return;
}


























__global__ void OSCARSSR_Cuda_FluxGPU3 (double *x, double *y, double *z, double *bx, double *by, double *bz, double *sx, double *sy, double *sz, double *dt, int *nt, int *ns, double *C0, double *C2, double *C, double *Omega, double *flux)
{
  // Check that this is within the number of spectrum points requested
  int is = threadIdx.x + blockIdx.x * blockDim.x;
  // Max number for shared memory

  int const NSHARED = 1000;

  // Define the shared memory
  __shared__ double sh_x[1000];
  __shared__ double sh_y[1000];
  __shared__ double sh_z[1000];
  __shared__ double sh_bx[1000];
  __shared__ double sh_by[1000];
  __shared__ double sh_bz[1000];


  // observer
  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);



  // Total number of copies
  int const NTotalCopies = *nt / NSHARED + 1;

  // icp is Copy Number
  for (int icp = 0; icp < NTotalCopies; ++icp) {

    // Offset instarting point in global array
    int const GlobalOffset = icp * NSHARED;

    __syncthreads();

    if (threadIdx.x == 0) {
      // icpth is the copy number in a thread
      for (int icpth = 0; icpth < NSHARED; ++icpth) {

        // index of *this* shared memory
        int const ThisThreadSharedIndex = icpth;

        // Global index of *this*
        int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

        // Check if we are within the nt bound
        if (GlobalIndex >= *nt) {
          break;
        }

        // Copy global memory to shared bank
        sh_x[ThisThreadSharedIndex]  = x[GlobalIndex];
        sh_y[ThisThreadSharedIndex]  = y[GlobalIndex];
        sh_z[ThisThreadSharedIndex]  = z[GlobalIndex];
        sh_bx[ThisThreadSharedIndex] = bx[GlobalIndex];
        sh_by[ThisThreadSharedIndex] = by[GlobalIndex];
        sh_bz[ThisThreadSharedIndex] = bz[GlobalIndex];
      }
    }
    __syncthreads();

    if (is >= *ns) {
      continue;
    }

    for (int ish = 0; ish < NSHARED; ++ish) {
      int const i = GlobalOffset + ish;

      if (i >= *nt) {
        break;
      }

      // Distance to observer
      double const D = sqrt( pow( (ox) - sh_x[ish], 2) + pow( (oy) - sh_y[ish], 2) + pow((oz) - sh_z[ish], 2) );

      // Normal in direction of observer
      double const NX = ((ox) - sh_x[ish]) / D;
      double const NY = ((oy) - sh_y[ish]) / D;
      double const NZ = ((oz) - sh_z[ish]) / D;

      // Exponent for fourier transformed field
      cuDoubleComplex Exponent = make_cuDoubleComplex(0, (*Omega) * ((*dt) * i + D / (*C)));

      cuDoubleComplex X1 = make_cuDoubleComplex((sh_bx[ish] - NX) / D, -(*C) * NX / ((*Omega) * D * D));
      cuDoubleComplex Y1 = make_cuDoubleComplex((sh_by[ish] - NY) / D, -(*C) * NY / ((*Omega) * D * D));
      cuDoubleComplex Z1 = make_cuDoubleComplex((sh_bz[ish] - NZ) / D, -(*C) * NZ / ((*Omega) * D * D));

      cuDoubleComplex MyEXP = cuCexp(Exponent);

      cuDoubleComplex X2 = cuCmul(X1, MyEXP);
      cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
      cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


      SumEX = cuCadd(SumEX, X2);
      SumEY = cuCadd(SumEY, Y2);
      SumEZ = cuCadd(SumEZ, Z2);
    }
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEZ);


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;


  flux[is] = (*C2) * (EX + EY + EZ);

  return;
}











__global__ void OSCARSSR_Cuda_FluxGPU2 (double *x, double *y, double *z, double *bx, double *by, double *bz, double *sx, double *sy, double *sz, double *dt, int *nt, int *ns, double *C0, double *C2, double *C, double *Omega, double *flux)
{
  // Check that this is within the number of spectrum points requested
  int is = threadIdx.x + blockIdx.x * blockDim.x;
  // Max number for shared memory
  int const NSHAREDMAX = 1000;

  // Number for each thread to copy from global to shared memory
  int const NToCopyPerThread = (int) NSHAREDMAX / NTHREADS_PER_BLOCK;

  // Actual number of elements in shared memory to use
  int const NSHARED = NToCopyPerThread * NTHREADS_PER_BLOCK;

  // Define the shared memory
  __shared__ double sh_x[NSHARED];
  __shared__ double sh_y[NSHARED];
  __shared__ double sh_z[NSHARED];
  __shared__ double sh_bx[NSHARED];
  __shared__ double sh_by[NSHARED];
  __shared__ double sh_bz[NSHARED];


  // observer
  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);



  // Total number of copies
  int const NTotalCopies = *nt / NSHARED + 1;

  // icp is Copy Number
  for (int icp = 0; icp < NTotalCopies; ++icp) {

    // Offset instarting point in global array
    int const GlobalOffset = icp * NSHARED;

    // Local offset for this thread in shared memory
    int const ThreadOffset = NToCopyPerThread * threadIdx.x;
    //flux[is] = shoffset; return;

    __syncthreads();

    // icpth is the copy number in a thread
    for (int icpth = 0; icpth < NToCopyPerThread; ++icpth) {

      // index of *this* shared memory
      int const ThisThreadSharedIndex = ThreadOffset + icpth;

      // Global index of *this*
      int const GlobalIndex = GlobalOffset + ThisThreadSharedIndex;

      // Check if we are within the nt bound
      if (GlobalIndex >= *nt) {
        break;
      }

      // Copy global memory to shared bank
      sh_x[ThisThreadSharedIndex]  = x[GlobalIndex];
      sh_y[ThisThreadSharedIndex]  = y[GlobalIndex];
      sh_z[ThisThreadSharedIndex]  = z[GlobalIndex];
      sh_bx[ThisThreadSharedIndex] = bx[GlobalIndex];
      sh_by[ThisThreadSharedIndex] = by[GlobalIndex];
      sh_bz[ThisThreadSharedIndex] = bz[GlobalIndex];
    }
    __syncthreads();

    if (is >= *ns) {
      continue;
    }

    for (int ish = 0; ish < NSHARED; ++ish) {
      int const i = GlobalOffset + ish;

      if (i >= *nt) {
        break;
      }



      // Distance to observer
      double const D = sqrt( pow( (ox) - sh_x[ish], 2) + pow( (oy) - sh_y[ish], 2) + pow((oz) - sh_z[ish], 2) );

      // Normal in direction of observer
      double const NX = ((ox) - sh_x[ish]) / D;
      double const NY = ((oy) - sh_y[ish]) / D;
      double const NZ = ((oz) - sh_z[ish]) / D;

      // Exponent for fourier transformed field
      cuDoubleComplex Exponent = make_cuDoubleComplex(0, (*Omega) * ((*dt) * i + D / (*C)));

      cuDoubleComplex X1 = make_cuDoubleComplex((sh_bx[ish] - NX) / D, -(*C) * NX / ((*Omega) * D * D));
      cuDoubleComplex Y1 = make_cuDoubleComplex((sh_by[ish] - NY) / D, -(*C) * NY / ((*Omega) * D * D));
      cuDoubleComplex Z1 = make_cuDoubleComplex((sh_bz[ish] - NZ) / D, -(*C) * NZ / ((*Omega) * D * D));

      cuDoubleComplex MyEXP = cuCexp(Exponent);

      cuDoubleComplex X2 = cuCmul(X1, MyEXP);
      cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
      cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


      SumEX = cuCadd(SumEX, X2);
      SumEY = cuCadd(SumEY, Y2);
      SumEZ = cuCadd(SumEZ, Z2);








    }
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEZ);


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;


  flux[is] = (*C2) * (EX + EY + EZ);

  return;
}








__global__ void OSCARSSR_Cuda_FluxGPU (double *x, double *y, double *z, double *bx, double *by, double *bz, double *sx, double *sy, double *sz, double *dt, int *nt, int *ns, double *C0, double *C2, double *C, double *Omega, double *flux)
{
  // Check that this is within the number of spectrum points requested
  int is = threadIdx.x + blockIdx.x * blockDim.x;
  if (is >= *ns) {
    return;
  }

  // Complex i
  cuDoubleComplex I = make_cuDoubleComplex(0, 1);

  cuDoubleComplex ICoverOmega = make_cuDoubleComplex(0, (*C) / (*Omega));

  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);


  // Loop over all points in trajectory
  for (int i = 0; i < *nt; ++i) {

    // Distance to observer
    double const D = sqrt( pow( (ox) - x[i], 2) + pow( (oy) - y[i], 2) + pow((oz) - z[i], 2) );

    // Normal in direction of observer
    double const NX = ((ox) - x[i]) / D;
    double const NY = ((oy) - y[i]) / D;
    double const NZ = ((oz) - z[i]) / D;

    // Exponent for fourier transformed field
    cuDoubleComplex Exponent = make_cuDoubleComplex(0, (*Omega) * ((*dt) * i + D / (*C)));

    cuDoubleComplex X1 = make_cuDoubleComplex((bx[i] - NX) / D, -(*C) * NX / ((*Omega) * D * D));
    cuDoubleComplex Y1 = make_cuDoubleComplex((by[i] - NY) / D, -(*C) * NY / ((*Omega) * D * D));
    cuDoubleComplex Z1 = make_cuDoubleComplex((bz[i] - NZ) / D, -(*C) * NZ / ((*Omega) * D * D));

    cuDoubleComplex MyEXP = cuCexp(Exponent);
    //cuDoubleComplex MyEXP = make_cuDoubleComplex( exp(Exponent.x) * cos(Exponent.y), exp(Exponent.x) * sin(Exponent.y));

    cuDoubleComplex X2 = cuCmul(X1, MyEXP);
    cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
    cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


    SumEX = cuCadd(SumEX, X2);
    SumEY = cuCadd(SumEY, Y2);
    SumEZ = cuCadd(SumEZ, Z2);

    // Sum in fourier transformed field (integral)
    //SumEX += (TVector3DC(B) - (N *     (One + (ICoverOmega / (D)))     )) / D * std::exp(Exponent);
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * (*Omega) * (*dt)), SumEZ);


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;

  // Multiply field by Constant C1 and time step
  //SumE *= C1 * DeltaT;

  // Set the flux for this frequency / energy point
  //Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);

  flux[is] = (*C2) * (EX + EY + EZ);

  return;
}


















extern "C" void OSCARSSR_Cuda_CalculateFluxGPU (TParticleA& Particle,
                                                TSurfacePoints const& Surface,
                                                double const Energy_eV,
                                                T3DScalarContainer& FluxContainer,
                                                std::string const& Polarization,
                                                double const Angle,
                                                TVector3D const& HorizontalDirection,
                                                TVector3D const& PropogationDirection,
                                                double const Weight)
{
  // Do the setup for and call the GPU calculation of flux.  Your limitation here is only GPU memory.

  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);
  if (ngpu == 0) {
    throw std::invalid_argument("No GPU found");
  }

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  int const NTPoints = (int) T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  double *x     = new double[NTPoints];
  double *y     = new double[NTPoints];
  double *z     = new double[NTPoints];
  double *bx    = new double[NTPoints];
  double *by    = new double[NTPoints];
  double *bz    = new double[NTPoints];


  int const NSPoints = (int) Surface.GetNPoints();

  // Observer
  double *sx     = new double[NSPoints];
  double *sy     = new double[NSPoints];
  double *sz     = new double[NSPoints];

  // Constants
  double const C = TOSCARSSR::C();
  double const Omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);

  // Flux
  double *flux = new double[NSPoints];


  // Set trajectory
  for (size_t i = 0; i < NTPoints; ++i) {
    x[i] = T.GetX(i).GetX();
    y[i] = T.GetX(i).GetY();
    z[i] = T.GetX(i).GetZ();

    bx[i] = T.GetB(i).GetX();
    by[i] = T.GetB(i).GetY();
    bz[i] = T.GetB(i).GetZ();
  }

  // Set the surface points
  for (size_t i = 0; i < NSPoints; ++i) {
    sx[i] = Surface.GetPoint(i).GetX();
    sy[i] = Surface.GetPoint(i).GetY();
    sz[i] = Surface.GetPoint(i).GetZ();
  }




  double *d_x, *d_y, *d_z;
  double *d_bx, *d_by, *d_bz;
  double *d_sx, *d_sy, *d_sz;
  double *d_flux;
  double *d_dt;
  int    *d_nt, *d_ns;

  int const size_x = NTPoints * sizeof(double);
  int const size_s = NSPoints * sizeof(double);

  cudaMalloc((void **) &d_x, size_x);
  cudaMalloc((void **) &d_y, size_x);
  cudaMalloc((void **) &d_z, size_x);

  cudaMalloc((void **) &d_bx, size_x);
  cudaMalloc((void **) &d_by, size_x);
  cudaMalloc((void **) &d_bz, size_x);

  cudaMalloc((void **) &d_sx, size_s);
  cudaMalloc((void **) &d_sy, size_s);
  cudaMalloc((void **) &d_sz, size_s);


  cudaMalloc((void **) &d_dt, sizeof(double));
  cudaMalloc((void **) &d_nt, sizeof(int));
  cudaMalloc((void **) &d_ns, sizeof(int));

  cudaMalloc((void **) &d_flux, size_s);


  cudaMemcpy(d_x, x, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, z, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_bx, bx, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_by, by, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bz, bz, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_sx, sx, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sy, sy, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sz, sz, size_s, cudaMemcpyHostToDevice);



  cudaMemcpy(d_dt, &DeltaT, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nt, &NTPoints, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ns, &NSPoints, sizeof(int), cudaMemcpyHostToDevice);



  // Constant C0 for calculation
  double const C0 = Particle.GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TOSCARSSR::FourPi() * Particle.GetCurrent() / (TOSCARSSR::H() * fabs(Particle.GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;

  // Constants to send in to GPU
  double *d_C0, *d_C2, *d_Omega, *d_C;

  // Allocate memory for constants
  cudaMalloc((void **) &d_C0,        sizeof(double));
  cudaMalloc((void **) &d_C2,        sizeof(double));
  cudaMalloc((void **) &d_Omega,     sizeof(double));
  cudaMalloc((void **) &d_C,         sizeof(double));

  // Copy constants to GPU
  cudaMemcpy(d_C0,        &C0,        sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C2,        &C2,        sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_Omega,     &Omega,     sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C,         &C,         sizeof(double), cudaMemcpyHostToDevice);


  // Send computation to gpu
  int const NBlocks = NSPoints / NTHREADS_PER_BLOCK + 1;
  OSCARSSR_Cuda_FluxGPU<<<NBlocks, NTHREADS_PER_BLOCK>>>(d_x, d_y, d_z, d_bx, d_by, d_bz, d_sx, d_sy, d_sz, d_dt, d_nt, d_ns, d_C0, d_C2, d_C, d_Omega, d_flux);

  // Copy result back from GPU
  cudaMemcpy(flux, d_flux, size_s, cudaMemcpyDeviceToHost);



  // Add result to power density container
  for (size_t i = 0; i < NSPoints; ++i) {
    FluxContainer.AddToPoint(i, flux[i] * Weight);
  }


  // Free all gpu memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);

  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);

  cudaFree(d_sx);
  cudaFree(d_sy);
  cudaFree(d_sz);

  cudaFree(d_dt);
  cudaFree(d_nt);
  cudaFree(d_ns);

  cudaFree(d_flux);

  cudaFree(d_C0);
  cudaFree(d_C2);
  cudaFree(d_Omega);
  cudaFree(d_C);


  // Free all heap memory
  delete [] x;
  delete [] y;
  delete [] z;

  delete [] bx;
  delete [] by;
  delete [] bz;

  delete [] sx;
  delete [] sy;
  delete [] sz;


  delete [] flux;


  return;
}






extern "C" void OSCARSSR_Cuda_CalculateFluxGPU2 (OSCARSSR& OSR,
                                                TSurfacePoints const& Surface,
                                                double const Energy_eV,
                                                T3DScalarContainer& FluxContainer,
                                                std::string const& Polarization,
                                                double const Angle,
                                                TVector3D const& HorizontalDirection,
                                                TVector3D const& PropogationDirection,
                                                int const NParticles,
                                                std::vector<int> const& GPUVector)
{
  // Calculate the flux for NParticles using the GPUs given in GPUVector.  Each particle's
  // trajectory will be sent to all GPUs for processing, meanwhile a new trajectory will
  // be calculated

  // Number of available GPUs
  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);
  if (ngpu == 0) {
    throw std::invalid_argument("No GPU found");
  }

  // Make sure that a gpu listed is within the range and not a duplicate
  std::vector<int> GPUsToUse;
  for (std::vector<int>::const_iterator it = GPUVector.begin(); it != GPUVector.end(); ++it) {
    if ( !(std::find(GPUsToUse.begin(), GPUsToUse.end(), *it) != GPUsToUse.end() && (*it < ngpu)) ) {
      GPUsToUse.push_back(*it);
    }
  }

  // Make sure we have at least one
  if (GPUsToUse.size() == 0) {
    throw std::invalid_argument("GPUs selected do not match hardware");
  }
  int const NGPUsToUse = (int) GPUsToUse.size();

  // Do we calculate for the current particle?
  bool const ThisParticleOnly = NParticles == 0 ? true : false;
  int  const NParticlesReally = ThisParticleOnly ? 1 : NParticles;
  // Type check, new particle if no type


  int *h_nt, *h_ns;
  double *h_dt;
  cudaHostAlloc((void**) &h_nt, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ns, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_dt, sizeof(double), cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt = (int) OSR.GetTrajectory().GetNPoints();
  *h_ns = (int) Surface.GetNPoints();
  *h_dt = (double) OSR.GetTrajectory().GetDeltaT();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = 32*16;
  int const NThreadsRemainder = NThreads % NThreadsPerBlock;
  int const NBlocksTotal = (NThreads - 1) / NThreadsPerBlock + 1;
  int const NBlocksPerGPU = NBlocksTotal / NGPUsToUse;
  int const NRemainderBlocks = NBlocksTotal % NGPUsToUse;
  // UPDATE: To be modified
  int const NFlux = NThreadsPerBlock * (NBlocksPerGPU + (NRemainderBlocks > 0 ? 1 : 0));

  std::vector<int> NBlocksThisGPU(NGPUsToUse, NBlocksPerGPU);
  for (int i = 0; i < NRemainderBlocks; ++i) {
    ++NBlocksThisGPU[i];
  }

  // Memory allocation for Host
  double  *h_x,  *h_y,  *h_z,  *h_bx,  *h_by,  *h_bz,  *h_sx,  *h_sy,  *h_sz,   *h_c0,  *h_c2,  *h_c,  *h_omega;
  int     *h_ifirst;
  double **h_flux;
  cudaHostAlloc((void**) &h_x,           *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,           *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,           *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,          *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,          *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,          *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sx,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sy,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sz,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c0,                  sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c2,                  sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c,                   sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_omega,               sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ifirst, NGPUsToUse * sizeof(int),     cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_flux,   NGPUsToUse * sizeof(double*), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    cudaHostAlloc((void**) &(h_flux[i]), NFlux * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  }

  // First surface point for each gpu
  int NBlocksUsed = 0;
  for (int i = 0; i < NGPUsToUse; ++i) {
    h_ifirst[i] = NBlocksUsed * NThreadsPerBlock;
    std::cout << "ifirst " << h_ifirst[i] << std::endl;
    NBlocksUsed += NBlocksThisGPU[i];
  }


  // Memor allocations for GPU
  int    **d_nt;
  int    **d_ns;
  double **d_dt;
  double **d_x;
  double **d_y;
  double **d_z;
  double **d_bx;
  double **d_by;
  double **d_bz;
  double **d_sx;
  double **d_sy;
  double **d_sz;
  double **d_c0;
  double **d_c2;
  double **d_c;
  double **d_omega;
  int    **d_ifirst;
  double **d_flux;

  cudaHostAlloc((void **) &d_nt,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ns,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_dt,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_x,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_y,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_z,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_by,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sy,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c0,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c2,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_omega,  NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ifirst, NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_flux,   NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaMalloc((void **) &d_nt[i],             sizeof(int));
    cudaMalloc((void **) &d_ns[i],             sizeof(int));
    cudaMalloc((void **) &d_dt[i],             sizeof(double));
    cudaMalloc((void **) &d_x[i],      *h_nt * sizeof(double));
    cudaMalloc((void **) &d_y[i],      *h_nt * sizeof(double));
    cudaMalloc((void **) &d_z[i],      *h_nt * sizeof(double));
    cudaMalloc((void **) &d_bx[i],     *h_nt * sizeof(double));
    cudaMalloc((void **) &d_by[i],     *h_nt * sizeof(double));
    cudaMalloc((void **) &d_bz[i],     *h_nt * sizeof(double));
    cudaMalloc((void **) &d_sx[i],     *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sy[i],     *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sz[i],     *h_ns * sizeof(double));
    cudaMalloc((void **) &d_c0[i],             sizeof(double));
    cudaMalloc((void **) &d_c2[i],             sizeof(double));
    cudaMalloc((void **) &d_c[i],              sizeof(double));
    cudaMalloc((void **) &d_omega[i],          sizeof(double));
    cudaMalloc((void **) &d_ifirst[i],         sizeof(int));
    cudaMalloc((void **) &d_flux[i],   NFlux * sizeof(double));

    // Copy device number to device
    cudaMemcpyAsync(d_ifirst[i], &(h_ifirst[i]), sizeof(int), cudaMemcpyHostToDevice);
  }

  // Compute known host values
  *h_c0    = OSR.GetCurrentParticle().GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());
  *h_c2    = TOSCARSSR::FourPi() * OSR.GetCurrentParticle().GetCurrent() / (TOSCARSSR::H() * fabs(OSR.GetCurrentParticle().GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;
  *h_c     = TOSCARSSR::C();
  *h_omega = TOSCARSSR::EvToAngularFrequency(Energy_eV);
  for (size_t i = 0; i < *h_ns; ++i) {
    h_sx[i] = Surface.GetPoint(i).GetX();
    h_sy[i] = Surface.GetPoint(i).GetY();
    h_sz[i] = Surface.GetPoint(i).GetZ();
  }

  // Copy constants to first device (async)
  int const d0 = GPUsToUse[0];
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_nt[0],    h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_ns[0],    h_ns,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_dt[0],    h_dt,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c0[0],    h_c0,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c2[0],    h_c2,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c[0],     h_c,           sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_omega[0], h_omega,       sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sx[0],    h_sx,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sy[0],    h_sy,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sz[0],    h_sz,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_nt[i+1],     d1, d_nt[i],     d, sizeof(int));
    cudaMemcpyPeerAsync( d_ns[i+1],     d1, d_ns[i],     d, sizeof(int));
    cudaMemcpyPeerAsync( d_dt[i+1],     d1, d_dt[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_c0[i+1],     d1, d_c0[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_c2[i+1],     d1, d_c2[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_c[i+1],      d1, d_c[i],      d, sizeof(double));
    cudaMemcpyPeerAsync( d_omega[i+1],  d1, d_omega[i],  d, sizeof(double));
    cudaMemcpyPeerAsync( d_sx[i+1],     d1, d_sx[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sy[i+1],     d1, d_sy[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sz[i+1],     d1, d_sz[i],     d, *h_ns * sizeof(double));
  }

  // Set first trajectory
  TParticleTrajectoryPoints const& T = OSR.GetTrajectory();
  for (size_t i = 0; i < *h_nt; ++i) {
    h_x[i]  = T.GetX(i).GetX();
    h_y[i]  = T.GetX(i).GetY();
    h_z[i]  = T.GetX(i).GetZ();
    h_bx[i] = T.GetB(i).GetX();
    h_by[i] = T.GetB(i).GetY();
    h_bz[i] = T.GetB(i).GetZ();
  }

  // Set the surface points
  // GPU events
  cudaEvent_t *event_fluxcopy = new cudaEvent_t[NGPUsToUse];
  for (int ig = 0; ig < NGPUsToUse; ++ig) {
    int const d = GPUsToUse[ig];
    cudaSetDevice(d);
    cudaEventCreate(&(event_fluxcopy[ig]));
  }

  // Enable peer (direct gpu-gpu) writes
  for (size_t ig = 0; ig < GPUsToUse.size() - 1; ++ig) {
    // Device number
    int const d  = GPUsToUse[ig];
    int const d1 = GPUsToUse[ig+1];
    int access;
    cudaDeviceCanAccessPeer(&access, d, d1);
    if (access == 1) {
      cudaSetDevice(d);
      cudaDeviceEnablePeerAccess(d1, 0);
    }
  }

  // Loop over number of particles
  for (int ip = 0; ip < NParticlesReally; ++ip) {

    // Copy trajectory to first GPU, then internal async transfers (where possible)
    cudaSetDevice(d0);
    cudaMemcpyAsync(d_x[0],  h_x,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_y[0],  h_y,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_z[0],  h_z,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bx[0], h_bx, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_by[0], h_by, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bz[0], h_bz, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    for (size_t ig = 0; ig < GPUsToUse.size() - 1; ++ig) {
      // Device number
      int const d  = GPUsToUse[ig];
      int const d1 = GPUsToUse[ig+1];
      cudaSetDevice(d);
      cudaMemcpyPeerAsync(d_x[ig+1],  d1, d_x[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_y[ig+1],  d1, d_y[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_z[ig+1],  d1, d_z[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bx[ig+1], d1, d_bx[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_by[ig+1], d1, d_by[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bz[ig+1], d1, d_bz[ig], d, *h_nt * sizeof(double));

    }

    // Wait for previous copy, start next one
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaEventSynchronize(event_fluxcopy[ig]);
      OSCARSSR_Cuda_FluxGPUMulti<<<NBlocksThisGPU[ig], NThreadsPerBlock>>>(d_x[ig], d_y[ig], d_z[ig], d_bx[ig], d_by[ig], d_bz[ig], d_sx[ig], d_sy[ig], d_sz[ig], d_dt[ig], d_nt[ig], d_ns[ig], d_c0[ig], d_c2[ig], d_c[ig], d_omega[ig], d_ifirst[ig], d_flux[ig]);
    }


    // Add result to flux container (from **previous**)
    if (ip > 0) {
      int NBlocksUsed = 0;
      for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
        for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
          if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
            break;
          }
          int iss = ith + NThreadsPerBlock * NBlocksUsed;
          FluxContainer.AddToPoint(iss, h_flux[ig][ith]);
        }
        NBlocksUsed += NBlocksThisGPU[ig];
      }
    }

    // Add copy back to streams
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d  = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaMemcpyAsync(h_flux[ig],  d_flux[ig],  NFlux * sizeof(double), cudaMemcpyDeviceToHost);
      cudaEventRecord(event_fluxcopy[ig]);
    }

    // If it's not the last one, calculate a new trajectory
    if (ip < NParticlesReally - 1) {
      OSR.SetNewParticle();
      OSR.CalculateTrajectory();
      TParticleTrajectoryPoints const& T = OSR.GetTrajectory();

      for (size_t it = 0; it < *h_nt; ++it) {
        h_x[it]  = T.GetX(it).GetX();
        h_y[it]  = T.GetX(it).GetY();
        h_z[it]  = T.GetX(it).GetZ();
        h_bx[it] = T.GetB(it).GetX();
        h_by[it] = T.GetB(it).GetY();
        h_bz[it] = T.GetB(it).GetZ();
      }
    }
  }

  // Wait for last copy
  for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
    cudaEventSynchronize(event_fluxcopy[ig]);
  }

  // Add result to flux container (from **previous**)
  NBlocksUsed = 0;
  for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
    for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
      if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
        break;
      }
      int iss = ith + NThreadsPerBlock * NBlocksUsed;
      FluxContainer.AddToPoint(iss, h_flux[ig][ith]);
    }
    NBlocksUsed += NBlocksThisGPU[ig];
  }

  // Weighting for multi-particle
  double const Weight = 1.0 / (double) NParticlesReally;
  FluxContainer.WeightAll(Weight);

  // Free host memory
  cudaFreeHost(h_nt);
  cudaFreeHost(h_ns);
  cudaFreeHost(h_dt);
  cudaFreeHost(h_x);
  cudaFreeHost(h_y);
  cudaFreeHost(h_z);
  cudaFreeHost(h_bx);
  cudaFreeHost(h_by);
  cudaFreeHost(h_bz);
  cudaFreeHost(h_sx);
  cudaFreeHost(h_sy);
  cudaFreeHost(h_sz);
  cudaFreeHost(h_c0);
  cudaFreeHost(h_c2);
  cudaFreeHost(h_c);
  cudaFreeHost(h_omega);
  cudaFreeHost(h_ifirst);
  // Free host and GPU memory
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {

    cudaFreeHost(h_flux[i]);

    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaFree(d_nt[i]);
    cudaFree(d_ns[i]);
    cudaFree(d_dt[i]);
    cudaFree(d_x[i]);
    cudaFree(d_y[i]);
    cudaFree(d_z[i]);
    cudaFree(d_bx[i]);
    cudaFree(d_by[i]);
    cudaFree(d_bz[i]);
    cudaFree(d_sx[i]);
    cudaFree(d_sy[i]);
    cudaFree(d_sz[i]);
    cudaFree(d_c0[i]);
    cudaFree(d_c2[i]);
    cudaFree(d_c[i]);
    cudaFree(d_omega[i]);
    cudaFree(d_ifirst[i]);
    cudaFree(d_flux[i]);
  }
  cudaFree(h_flux);

  cudaFree(d_nt);
  cudaFree(d_ns);
  cudaFree(d_dt);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);
  cudaFree(d_sx);
  cudaFree(d_sy);
  cudaFree(d_sz);
  cudaFree(d_c0);
  cudaFree(d_c2);
  cudaFree(d_c);
  cudaFree(d_omega);
  cudaFree(h_ifirst);
  cudaFree(d_flux);

  // Delete host gpu pointer arrays
  delete [] event_fluxcopy;


  return;
}



































__global__ void OSCARSSR_Cuda_SpectrumGPU (double *x, double *y, double *z, double *bx, double *by, double *bz, double *ox, double *oy, double *oz, double *dt, int *nt, int *ns, double *C0, double *C2, double *EvToOmega, double *C, double *se, double *sf, cuDoubleComplex* pol, int *pol_state)
{
  // Check that this is within the number of spectrum points requested
  int is = threadIdx.x + blockIdx.x * blockDim.x;
  if (is >= *ns) {
    return;
  }

  // Complex i
  cuDoubleComplex I = make_cuDoubleComplex(0, 1);

  double const Omega = *EvToOmega * se[is];
  cuDoubleComplex ICoverOmega = make_cuDoubleComplex(0, (*C) / Omega);


  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);


  // Loop over all points in trajectory
  for (int i = 0; i < *nt; ++i) {

    // Distance to observer
    double const D = sqrt( pow( (*ox) - x[i], 2) + pow( (*oy) - y[i], 2) + pow((*oz) - z[i], 2) );

    // Normal in direction of observer
    double const NX = ((*ox) - x[i]) / D;
    double const NY = ((*oy) - y[i]) / D;
    double const NZ = ((*oz) - z[i]) / D;

    // Exponent for fourier transformed field
    cuDoubleComplex Exponent = make_cuDoubleComplex(0, Omega * ((*dt) * i + D / (*C)));

    cuDoubleComplex X1 = make_cuDoubleComplex((bx[i] - NX) / D, -(*C) * NX / (Omega * D * D));
    cuDoubleComplex Y1 = make_cuDoubleComplex((by[i] - NY) / D, -(*C) * NY / (Omega * D * D));
    cuDoubleComplex Z1 = make_cuDoubleComplex((bz[i] - NZ) / D, -(*C) * NZ / (Omega * D * D));

    cuDoubleComplex MyEXP = cuCexp(Exponent);
    //cuDoubleComplex MyEXP = make_cuDoubleComplex( exp(Exponent.x) * cos(Exponent.y), exp(Exponent.x) * sin(Exponent.y));

    cuDoubleComplex X2 = cuCmul(X1, MyEXP);
    cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
    cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


    SumEX = cuCadd(SumEX, X2);
    SumEY = cuCadd(SumEY, Y2);
    SumEZ = cuCadd(SumEZ, Z2);

    // Sum in fourier transformed field (integral)
    //SumEX += (TVector3DC(B) - (N *     (One + (ICoverOmega / (D)))     )) / D * std::exp(Exponent);
  }

  SumEX = cuCmul(make_cuDoubleComplex(0, (*C0) * Omega * (*dt)), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex(0, (*C0) * Omega * (*dt)), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex(0, (*C0) * Omega * (*dt)), SumEZ);


  // Check for polarization state
  if (*pol_state == 0) {
    // Do nothing
  } else if (*pol_state == 1) {
    // Linear, just dot with vector and put in direction of vector
    cuDoubleComplex Magnitude = cuCadd(cuCadd(cuCmul(SumEX, pol[0]), cuCmul(SumEY, pol[1])),  cuCmul(SumEZ, pol[2]));
    SumEX = cuCmul(Magnitude, pol[0]);
    SumEY = cuCmul(Magnitude, pol[1]);
    SumEZ = cuCmul(Magnitude, pol[2]);
  } else if (*pol_state == 2) {
    cuDoubleComplex Magnitude = cuCadd(cuCadd(cuCmul(SumEX, cuConj(pol[0])), cuCmul(SumEY, cuConj(pol[1]))),  cuCmul(SumEZ, cuConj(pol[2])));
    SumEX = cuCmul(Magnitude, pol[0]);
    SumEY = cuCmul(Magnitude, pol[1]);
    SumEZ = cuCmul(Magnitude, pol[2]);
  } else {
    // UPDATE: Serious problem
  }


  double const EX = SumEX.x * SumEX.x + SumEX.y * SumEX.y;
  double const EY = SumEY.x * SumEY.x + SumEY.y * SumEY.y;
  double const EZ = SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y;

  // Multiply field by Constant C1 and time step
  //SumE *= C1 * DeltaT;

  // Set the flux for this frequency / energy point
  //Spectrum.AddToFlux(i, C2 *  SumE.Dot( SumE.CC() ).real() * Weight);

  sf[is] = (*C2) * (EX + EY + EZ);

  return;
}








extern "C" void OSCARSSR_Cuda_CalculateSpectrumGPU (TParticleA& Particle,
                                                    TVector3D const& ObservationPoint,
                                                    TSpectrumContainer& Spectrum,
                                                    std::string const& Polarization,
                                                    double const Angle,
                                                    TVector3D const& HorizontalDirection,
                                                    TVector3D const& PropogationDirection,
                                                    double const Weight)
{

  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);
  if (ngpu == 0) {
    throw std::invalid_argument("No GPU found");
  }

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  int const NTPoints = (int) T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  double *x     = new double[NTPoints];
  double *y     = new double[NTPoints];
  double *z     = new double[NTPoints];
  double *bx    = new double[NTPoints];
  double *by    = new double[NTPoints];
  double *bz    = new double[NTPoints];

  int const NSPoints = (int) Spectrum.GetNPoints();

  // Observer
  double ox = ObservationPoint.GetX();
  double oy = ObservationPoint.GetY();
  double oz = ObservationPoint.GetZ();

  // Constants
  double const C = TOSCARSSR::C();
  double const EvToOmega = TOSCARSSR::EvToAngularFrequency(1);

  // Spectrum energy and flux
  double *se     = new double[NSPoints];
  double *sf     = new double[NSPoints];



  // Imaginary "i" and complxe 1+0i
  std::complex<double> const I(0, 1);
  std::complex<double> const One(1, 0);

  // Photon vertical direction and positive and negative helicity
  TVector3D const VerticalDirection = PropogationDirection.Cross(HorizontalDirection).UnitVector();
  TVector3DC const Positive = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) + VerticalDirection * I );
  TVector3DC const Negative = 1. / sqrt(2) * (TVector3DC(HorizontalDirection) - VerticalDirection * I );

  // For polarization input to the gpu
  cuDoubleComplex *pol = new cuDoubleComplex[3];

  // State of polarization: 0 for all, 1 for linear, 2 for circular
  // (requires different threatment of vector pol interally)
  int pol_state = 1;

  if (Polarization == "all") {
    // Do nothing, it is already ALL
    pol_state = 0;
  } else if (Polarization == "linear-horizontal") {
    pol[0] = make_cuDoubleComplex(HorizontalDirection.GetX(), 0);
    pol[1] = make_cuDoubleComplex(HorizontalDirection.GetY(), 0);
    pol[2] = make_cuDoubleComplex(HorizontalDirection.GetZ(), 0);
  } else if (Polarization == "linear-vertical") {
    pol[0] = make_cuDoubleComplex(VerticalDirection.GetX(), 0);
    pol[1] = make_cuDoubleComplex(VerticalDirection.GetY(), 0);
    pol[2] = make_cuDoubleComplex(VerticalDirection.GetZ(), 0);
  } else if (Polarization == "linear") {
    TVector3D PolarizationAngle = HorizontalDirection;
    PolarizationAngle.RotateSelf(Angle, PropogationDirection);
    pol[0] = make_cuDoubleComplex(PolarizationAngle.GetX(), 0);
    pol[1] = make_cuDoubleComplex(PolarizationAngle.GetY(), 0);
    pol[2] = make_cuDoubleComplex(PolarizationAngle.GetZ(), 0);
  } else if (Polarization == "circular-left") {
    //SumE = SumE.Dot(Positive.CC()) * Positive;
    pol_state = 2;
    pol[0] = make_cuDoubleComplex(Positive.CC().GetX().real(), Positive.CC().GetX().imag());
    pol[1] = make_cuDoubleComplex(Positive.CC().GetY().real(), Positive.CC().GetY().imag());
    pol[2] = make_cuDoubleComplex(Positive.CC().GetZ().real(), Positive.CC().GetZ().imag());
  } else if (Polarization == "circular-right") {
    //SumE = SumE.Dot(Negative.CC()) * Negative;
    pol_state = 2;
    pol[0] = make_cuDoubleComplex(Negative.CC().GetX().real(), Negative.CC().GetX().imag());
    pol[1] = make_cuDoubleComplex(Negative.CC().GetY().real(), Negative.CC().GetY().imag());
    pol[2] = make_cuDoubleComplex(Negative.CC().GetZ().real(), Negative.CC().GetZ().imag());
  } else {
    // Throw invalid argument if polarization is not recognized
    //throw std::invalid_argument("Polarization requested not recognized");
  }

  // Set trajectory
  for (size_t i = 0; i < NTPoints; ++i) {
    x[i] = T.GetX(i).GetX();
    y[i] = T.GetX(i).GetY();
    z[i] = T.GetX(i).GetZ();

    bx[i] = T.GetB(i).GetX();
    by[i] = T.GetB(i).GetY();
    bz[i] = T.GetB(i).GetZ();
  }



  // Set energy to value
  for (size_t i = 0; i < NSPoints; ++i) {
    se[i] = Spectrum.GetEnergy(i);
  }



  double *d_x, *d_y, *d_z;
  double *d_bx, *d_by, *d_bz;
  double *d_ox, *d_oy, *d_oz;
  double *d_se, *d_sf;
  double *d_dt;
  int    *d_nt, *d_ns;

  // For polarization
  cuDoubleComplex *d_pol;
  int *d_pol_state;


  int const size_x = NTPoints * sizeof(double);
  int const size_s = NSPoints * sizeof(double);

  cudaMalloc((void **) &d_x, size_x);
  cudaMalloc((void **) &d_y, size_x);
  cudaMalloc((void **) &d_z, size_x);

  cudaMalloc((void **) &d_bx, size_x);
  cudaMalloc((void **) &d_by, size_x);
  cudaMalloc((void **) &d_bz, size_x);

  cudaMalloc((void **) &d_ox, sizeof(double));
  cudaMalloc((void **) &d_oy, sizeof(double));
  cudaMalloc((void **) &d_oz, sizeof(double));

  cudaMalloc((void **) &d_dt, sizeof(double));
  cudaMalloc((void **) &d_nt, sizeof(int));
  cudaMalloc((void **) &d_ns, sizeof(int));

  cudaMalloc((void **) &d_se, size_s);
  cudaMalloc((void **) &d_sf, size_s);

  // Polarization
  cudaMalloc((void **) &d_pol, 3*sizeof(cuDoubleComplex));
  cudaMalloc((void **) &d_pol_state, sizeof(int));


  cudaMemcpy(d_x, x, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, z, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_bx, bx, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_by, by, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bz, bz, size_x, cudaMemcpyHostToDevice);


  cudaMemcpy(d_ox, &ox, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_oy, &oy, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_oz, &oz, sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(d_dt, &DeltaT, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nt, &NTPoints, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ns, &NSPoints, sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(d_se, se, size_s, cudaMemcpyHostToDevice);

  cudaMemcpy(d_pol, pol, 3*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
  cudaMemcpy(d_pol_state, &pol_state, sizeof(int), cudaMemcpyHostToDevice);


  // Constant C0 for calculation
  double const C0 = Particle.GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());

  // Constant for flux calculation at the end
  double const C2 = TOSCARSSR::FourPi() * Particle.GetCurrent() / (TOSCARSSR::H() * fabs(Particle.GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;

  // Constants to send in to GPU
  double *d_C0, *d_C2, *d_EvToOmega, *d_C;

  // Allocate memory for constants
  cudaMalloc((void **) &d_C0,        sizeof(double));
  cudaMalloc((void **) &d_C2,        sizeof(double));
  cudaMalloc((void **) &d_EvToOmega, sizeof(double));
  cudaMalloc((void **) &d_C,         sizeof(double));

  // Copy constants to GPU
  cudaMemcpy(d_C0,        &C0,        sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C2,        &C2,        sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_EvToOmega, &EvToOmega, sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_C,         &C,         sizeof(double), cudaMemcpyHostToDevice);


  // Send computation to gpu
  int const NBlocks = NSPoints / NTHREADS_PER_BLOCK + 1;
  OSCARSSR_Cuda_SpectrumGPU<<<NBlocks, NTHREADS_PER_BLOCK>>>(d_x, d_y, d_z, d_bx, d_by, d_bz, d_ox, d_oy, d_oz, d_dt, d_nt, d_ns, d_C0, d_C2, d_EvToOmega, d_C, d_se, d_sf, d_pol, d_pol_state);

  // Copy result back from GPU
  cudaMemcpy(sf, d_sf, size_s, cudaMemcpyDeviceToHost);



  // Add result to power density container
  for (size_t i = 0; i < NSPoints; ++i) {
    Spectrum.AddToFlux(i, sf[i] * Weight);
  }


  // Free all gpu memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);

  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);

  cudaFree(d_ox);
  cudaFree(d_oy);
  cudaFree(d_oz);

  cudaFree(d_dt);
  cudaFree(d_nt);
  cudaFree(d_ns);

  cudaFree(d_se);
  cudaFree(d_sf);

  cudaFree(d_C0);
  cudaFree(d_C2);
  cudaFree(d_EvToOmega);
  cudaFree(d_C);

  cudaFree(d_pol);
  cudaFree(d_pol_state);


  // Free all heap memory
  delete [] x;
  delete [] y;
  delete [] z;

  delete [] bx;
  delete [] by;
  delete [] bz;


  delete [] se;
  delete [] sf;

  delete [] pol;


  return;
}


















__global__ void OSCARSSR_Cuda_PowerDensityGPU (double *x, double *y, double *z, double *bx, double *by, double *bz, double *aocx, double *aocy, double *aocz, double *sx, double *sy, double *sz, double *snx, double *sny, double *snz, double *dt, int *nt, int *ns, double *power_density)
{
  // Get surface id from block and thread number
  int is = threadIdx.x + blockIdx.x * blockDim.x;

  if (is >= *ns) {
    return;
  }




  // If you could copy int ultra-fast memory, something like this:
  //__shared__ double temp[6144];
  //if (threadIdx.x == 0) {
  //  for (int i = 0; i < *nt; ++i) {
  //    if (i <= 6144) {
  //      break;
  //    }
  //    temp[i] = x[i];
  //  }
  //}
  // __syncthreads();



  // Observation point
  double const OX = sx[is];
  double const OY = sy[is];
  double const OZ = sz[is];

  // Normal vector from input
  double const NormalX = snx[is];
  double const NormalY = sny[is];
  double const NormalZ = snz[is];

  double Sum = 0;

  for (int i = 0; i < *nt; ++i) {

    // Normal vector in direction of observation point
    double const R1 = sqrt( pow(OX - x[i], 2) + pow(OY - y[i], 2) + pow(OZ - z[i], 2) );
    double const N1X = (OX - x[i]) / R1;
    double const N1Y = (OY - y[i]) / R1;
    double const N1Z = (OZ - z[i]) / R1;

    // Surface normal dot with vector normal
    double const N1DotNormal = N1X * NormalX + N1Y * NormalY + N1Z * NormalZ;

    // Orthogonal vector 2 & 3
    double N2X;
    double N2Y;
    double N2Z;

    double const xx = N1X < 0.0 ? -N1X : N1X;
    double const yy = N1Y < 0.0 ? -N1Y : N1Y;
    double const zz = N1Z < 0.0 ? -N1Z : N1Z;
    if (xx < yy) {
      if (xx < zz) {
        N2X = 0;
        N2Y = N1Z;
        N2Z = -N1Y;
      } else {
        N2X = N1Y;
        N2Y = -N1X;
        N2Z = 0;
      }
    } else {
      if (yy < zz) {
        N2X = -N1Z;
        N2Y = 0;
        N2Z = N1X;
      } else {
        N2X = N1Y;
        N2Y = -N1X;
        N2Z = 0;
      }
    }
    double const R2 = sqrt(N2X * N2X + N2Y * N2Y + N2Z * N2Z);
    N2X /= R2;
    N2Y /= R2;
    N2Z /= R2;

    // Ortohgonal vector N3
    double const N3X = N1Y * N2Z - N1Z * N2Y;
    double const N3Y = N1Z * N2X - N1X * N2Z;
    double const N3Z = N1X * N2Y - N1Y * N2X;





    double const x1 = N1X - bx[i];
    double const y1 = N1Y - by[i];
    double const z1 = N1Z - bz[i];

    double const x2 = y1 * aocz[i] - z1 * aocy[i];
    double const y2 = z1 * aocx[i] - x1 * aocz[i];
    double const z2 = x1 * aocy[i] - y1 * aocx[i];

    // Numerator = N1.Cross( ( (N1 - B).Cross((AoverC)) ) );
    double const x3 = N1Y * z2 - N1Z * y2;
    double const y3 = N1Z * x2 - N1X * z2;
    double const z3 = N1X * y2 - N1Y * x2;

    double const BdotN1 = bx[i] * N1X + by[i] * N1Y + bz[i] * N1Z;
    double const Denominator = pow(1. - BdotN1, 5);

    Sum += pow( x3 * N2X + y3 * N2Y + z3 * N2Z, 2) / Denominator / (R1 * R1) * N1DotNormal;
    Sum += pow( x3 * N3X + y3 * N3Y + z3 * N3Z, 2) / Denominator / (R1 * R1) * N1DotNormal;
  }

  power_density[is] = Sum * (*dt);

  return;
}



extern "C" void OSCARSSR_Cuda_CalculatePowerDensityGPU (TParticleA& Particle,
                                                        TSurfacePoints const& Surface,
                                                        T3DScalarContainer& PowerDensityContainer,
                                                        bool const Directional,
                                                        double const Weight)
{

  int ngpu = 0;
  cudaGetDeviceCount(&ngpu);
  if (ngpu == 0) {
    throw std::invalid_argument("No GPU found");
  }

  // Grab the Trajectory
  TParticleTrajectoryPoints& T = Particle.GetTrajectory();

  // Number of points in Trajectory
  int const NTPoints = (int) T.GetNPoints();

  // Timestep from trajectory
  double const DeltaT = T.GetDeltaT();

  double *x     = new double[NTPoints];
  double *y     = new double[NTPoints];
  double *z     = new double[NTPoints];
  double *bx    = new double[NTPoints];
  double *by    = new double[NTPoints];
  double *bz    = new double[NTPoints];
  double *aocx  = new double[NTPoints];
  double *aocy  = new double[NTPoints];
  double *aocz  = new double[NTPoints];

  int const NSPoints = (int) Surface.GetNPoints();

  double *sx     = new double[NSPoints];
  double *sy     = new double[NSPoints];
  double *sz     = new double[NSPoints];

  double *snx    = new double[NSPoints];
  double *sny    = new double[NSPoints];
  double *snz    = new double[NSPoints];

  double *power_density = new double[NSPoints];


  for (size_t i = 0; i < NTPoints; ++i) {
    x[i] = T.GetX(i).GetX();
    y[i] = T.GetX(i).GetY();
    z[i] = T.GetX(i).GetZ();

    bx[i] = T.GetB(i).GetX();
    by[i] = T.GetB(i).GetY();
    bz[i] = T.GetB(i).GetZ();

    aocx[i] = T.GetAoverC(i).GetX();
    aocy[i] = T.GetAoverC(i).GetY();
    aocz[i] = T.GetAoverC(i).GetZ();
  }



  for (size_t i = 0; i < NSPoints; ++i) {
    sx[i] = Surface.GetPoint(i).GetX();
    sy[i] = Surface.GetPoint(i).GetY();
    sz[i] = Surface.GetPoint(i).GetZ();

    snx[i] = Surface.GetPoint(i).GetNormalX();
    sny[i] = Surface.GetPoint(i).GetNormalY();
    snz[i] = Surface.GetPoint(i).GetNormalZ();
  }



  double *d_x, *d_y, *d_z;
  double *d_bx, *d_by, *d_bz;
  double *d_aocx, *d_aocy, *d_aocz;
  double *d_sx, *d_sy, *d_sz;
  double *d_snx, *d_sny, *d_snz;
  double *d_power_density;
  double *d_dt;
  int    *d_nt, *d_ns;

  int const size_x = NTPoints * sizeof(double);
  int const size_s = NSPoints * sizeof(double);

  cudaMalloc((void **) &d_x, size_x);
  cudaMalloc((void **) &d_y, size_x);
  cudaMalloc((void **) &d_z, size_x);

  cudaMalloc((void **) &d_bx, size_x);
  cudaMalloc((void **) &d_by, size_x);
  cudaMalloc((void **) &d_bz, size_x);

  cudaMalloc((void **) &d_aocx, size_x);
  cudaMalloc((void **) &d_aocy, size_x);
  cudaMalloc((void **) &d_aocz, size_x);

  cudaMalloc((void **) &d_sx, size_s);
  cudaMalloc((void **) &d_sy, size_s);
  cudaMalloc((void **) &d_sz, size_s);

  cudaMalloc((void **) &d_snx, size_s);
  cudaMalloc((void **) &d_sny, size_s);
  cudaMalloc((void **) &d_snz, size_s);

  cudaMalloc((void **) &d_power_density, size_s);

  cudaMalloc((void **) &d_dt, sizeof(double));

  cudaMalloc((void **) &d_nt, sizeof(int));
  cudaMalloc((void **) &d_ns, sizeof(int));


  cudaMemcpy(d_x, x, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, y, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, z, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_bx, bx, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_by, by, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_bz, bz, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_aocx, aocx, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_aocy, aocy, size_x, cudaMemcpyHostToDevice);
  cudaMemcpy(d_aocz, aocz, size_x, cudaMemcpyHostToDevice);

  cudaMemcpy(d_sx, sx, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sy, sy, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sz, sz, size_s, cudaMemcpyHostToDevice);

  cudaMemcpy(d_snx, snx, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sny, sny, size_s, cudaMemcpyHostToDevice);
  cudaMemcpy(d_snz, snz, size_s, cudaMemcpyHostToDevice);

  cudaMemcpy(d_dt, &DeltaT, sizeof(double), cudaMemcpyHostToDevice);

  cudaMemcpy(d_nt, &NTPoints, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_ns, &NSPoints, sizeof(int), cudaMemcpyHostToDevice);


  // Send computation to gpu
  int const NBlocks = NSPoints / NTHREADS_PER_BLOCK + 1;
  OSCARSSR_Cuda_PowerDensityGPU<<<NBlocks, NTHREADS_PER_BLOCK>>>(d_x, d_y, d_z, d_bx, d_by, d_bz, d_aocx, d_aocy, d_aocz, d_sx, d_sy, d_sz, d_snx, d_sny, d_snz, d_dt, d_nt, d_ns, d_power_density);

  // Copy result back from GPU
  cudaMemcpy(power_density, d_power_density, size_s, cudaMemcpyDeviceToHost);



  // Add result to power density container
  for (size_t i = 0; i < NSPoints; ++i) {
    PowerDensityContainer.AddToPoint(i, power_density[i] * fabs(Particle.GetQ() * Particle.GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) / 1e6 * Weight);
  }


  // Free all gpu memory
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);

  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);

  cudaFree(d_aocx);
  cudaFree(d_aocy);
  cudaFree(d_aocz);

  cudaFree(d_sx);
  cudaFree(d_sy);
  cudaFree(d_sz);

  cudaFree(d_snx);
  cudaFree(d_sny);
  cudaFree(d_snz);

  cudaFree(d_dt);
  cudaFree(d_nt);
  cudaFree(d_ns);

  cudaFree(d_power_density);





  // Free all heap memory
  delete [] x;
  delete [] y;
  delete [] z;

  delete [] bx;
  delete [] by;
  delete [] bz;

  delete [] aocx;
  delete [] aocy;
  delete [] aocz;

  delete [] sx;
  delete [] sy;
  delete [] sz;

  delete [] snx;
  delete [] sny;
  delete [] snz;

  delete [] power_density;

  return;
}





