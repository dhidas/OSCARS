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
  static int ngpu = 0;
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




__device__ static __inline__ void GetInterpolatingIMinIMax(double *fx, int* nx, double* x, int *imin, int* imax)
{
  // get the klo and khi for interpolation

  int klo=0;
  int khi = *nx - 1;
  int k;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (fx[k] > *x) {
      khi = k;
    } else {
      klo = k;
    }
  }

  *imin = klo;
  *imax = khi;

  return;
}








__device__ static __inline__ double Interpolate (double *fx, double* fy, double *fypp, double* x, int* klo, int* khi)
{
  // Return the Y-value according to spline


  // Distance between points, check that it isn't zero!
  double const h = fx[*khi] - fx[*klo];
  if (h == 0) {
    // UPDATE: supposed to throw CUDA
  }

  // Fractional distance to the points on either side
  double const a = (fx[*khi] - *x) / h;
  double const b = (*x - fx[*klo]) / h;

  // Return the value of Y
  return a * fy[*klo] + b * fy[*khi] + ((a * a * a - a) * fypp[*klo] + (b * b * b - b) * fypp[*khi]) * (h * h) / 6.;
}











__device__ static __inline__ double Interpolate (double *fx, int* nx, double* fy, double *fypp, double* x)
{
  // Return the Y-value according to spline

  int klo=0;
  int khi = *nx - 1;
  int k;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (fx[k] > *x) {
      khi = k;
    } else {
      klo = k;
    }
  }

  // Distance between points, check that it isn't zero!
  double const h = fx[khi] - fx[klo];
  if (h == 0) {
    // UPDATE: supposed to throw CUDA
  }

  // Fractional distance to the points on either side
  double const a = (fx[khi] - *x) / h;
  double const b = (*x - fx[klo]) / h;

  // Return the value of Y
  return a * fy[klo] + b * fy[khi] + ((a * a * a - a) * fypp[klo] + (b * b * b - b) * fypp[khi]) * (h * h) / 6.;
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








__global__ void OSCARSSR_Cuda_FluxGPUMultiWithAInterpolated (double *t,
                                                             double  *x,   double *y,   double *z,
                                                             double  *xp,  double *yp,  double *zp,
                                                             double  *bx,  double *by,  double *bz,
                                                             double  *bxp, double *byp, double *bzp,
                                                             double  *ax,  double *ay,  double *az,
                                                             double  *axp, double *ayp, double *azp,
                                                             double  *sx,  double *sy,  double *sz,
                                                             double  *tstart, double *tstop,
                                                             int *nt,
                                                             int *ns,
                                                             double *C0,   double *C2,  double *C,
                                                             double *Omega,
                                                             int *ifirst,
                                                             int *ml,
                                                             double *prec,
                                                             double *flux)
{
  // Check that this is within the number of spectrum points requested
  int const ith = threadIdx.x + blockIdx.x * blockDim.x;
  int const is = ith + *ifirst;


  double const ox = is >= *ns ? 0 : sx[is];
  double const oy = is >= *ns ? 0 : sy[is];
  double const oz = is >= *ns ? 0 : sz[is];

  // E-field components sum
  cuDoubleComplex SumEX = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEY = make_cuDoubleComplex(0, 0);
  cuDoubleComplex SumEZ = make_cuDoubleComplex(0, 0);



  __shared__ double _t[NTHREADS_PER_BLOCK];
  __shared__ double _x[NTHREADS_PER_BLOCK];
  __shared__ double _y[NTHREADS_PER_BLOCK];
  __shared__ double _z[NTHREADS_PER_BLOCK];
  __shared__ double _bx[NTHREADS_PER_BLOCK];
  __shared__ double _by[NTHREADS_PER_BLOCK];
  __shared__ double _bz[NTHREADS_PER_BLOCK];
  __shared__ double _ax[NTHREADS_PER_BLOCK];
  __shared__ double _ay[NTHREADS_PER_BLOCK];
  __shared__ double _az[NTHREADS_PER_BLOCK];

  __shared__ bool   _done[NTHREADS_PER_BLOCK];
  __shared__ bool   _all_done;


  if (threadIdx.x == 1) {
    _all_done = false;
  }

  // initialize _done.  If not a surface point, we're going to use
  // the thread to do trajectory calculation anyways
  _done[threadIdx.x] = is >= *ns ? true : false;
  
  __syncthreads();

  int this_nt = 1;

  double this_result = 0;
  double last_result = 1;

  double result = -1;
  double dt_total = 0;
  for (int ilevel = 0; !_all_done && (ilevel <= *ml); ++ilevel) {

    // DeltaT inclusive up to this level
    dt_total = (*tstop - *tstart) / pow(2., ilevel+1);//(*tstop - *tstart) / (2 * this_nt);

    // deltaT this level and Time start this level
    double const dt = (*tstop - *tstart) / pow(2., ilevel);//(*tstop - *tstart) / this_nt;
    double const ts = *tstart + (*tstop - *tstart) / pow(2., ilevel + 1);//*tstart + (*tstop - *tstart) / (2. * this_nt);

    int const NTrajectoryBlocks = this_nt / blockDim.x + (this_nt % blockDim.x == 0 ? 0 : 1);

    for (int itb = 0; itb < NTrajectoryBlocks; ++itb) {
      _t[threadIdx.x] = dt * (itb * blockDim.x + threadIdx.x) + ts;

      if (_t[threadIdx.x] < *tstop) {
        int imin, imax;
        GetInterpolatingIMinIMax(t, nt, &(_t[threadIdx.x]), &imin, &imax);

        _x[threadIdx.x]  = Interpolate(t,  x,  xp, &(_t[threadIdx.x]), &imin, &imax);
        _y[threadIdx.x]  = Interpolate(t,  y,  yp, &(_t[threadIdx.x]), &imin, &imax);
        _z[threadIdx.x]  = Interpolate(t,  z,  zp, &(_t[threadIdx.x]), &imin, &imax);
        _bx[threadIdx.x] = Interpolate(t, bx, bxp, &(_t[threadIdx.x]), &imin, &imax);
        _by[threadIdx.x] = Interpolate(t, by, byp, &(_t[threadIdx.x]), &imin, &imax);
        _bz[threadIdx.x] = Interpolate(t, bz, bzp, &(_t[threadIdx.x]), &imin, &imax);
        _ax[threadIdx.x] = Interpolate(t, ax, axp, &(_t[threadIdx.x]), &imin, &imax);
        _ay[threadIdx.x] = Interpolate(t, ay, ayp, &(_t[threadIdx.x]), &imin, &imax);
        _az[threadIdx.x] = Interpolate(t, az, azp, &(_t[threadIdx.x]), &imin, &imax);
      }

      __syncthreads();

      if (!_done[threadIdx.x]) {
      for (int i = 0; i < blockDim.x; ++i) {

        // Check if we are over the limit of trajectory points
        if (is < *ns && (_t[i] < *tstop)) {

        // DO MATH HERE
        // Distance to observer
        double const D = sqrt( pow( (ox) - _x[i], 2) + pow( (oy) - _y[i], 2) + pow((oz) - _z[i], 2) );

        // Normal in direction of observer
        double const NX = ((ox) - _x[i]) / D;
        double const NY = ((oy) - _y[i]) / D;
        double const NZ = ((oz) - _z[i]) / D;

        // Magnitude of Beta squared
        double const One_Minus_BMag2 = 1. -  (_bx[i] * _bx[i] + _by[i] * _by[i] + _bz[i] * _bz[i]);

        // N dot Beta
        double const NDotBeta = NX * _bx[i] + NY * _by[i] + NZ * _bz[i];

        double const FarFieldDenominator =  D * (pow(1. - NDotBeta, 2));
        double const NearFieldDenominator = D * FarFieldDenominator;
        double const NearField_X = One_Minus_BMag2 * (NX - _bx[i]) / NearFieldDenominator;
        double const NearField_Y = One_Minus_BMag2 * (NY - _by[i]) / NearFieldDenominator;
        double const NearField_Z = One_Minus_BMag2 * (NZ - _bz[i]) / NearFieldDenominator;

        double const FFX = (NY - _by[i]) * _az[i] - (NZ - _bz[i]) * _ay[i];
        double const FFY = (NZ - _bz[i]) * _ax[i] - (NX - _bx[i]) * _az[i];
        double const FFZ = (NX - _bx[i]) * _ay[i] - (NY - _by[i]) * _ax[i];

        double const FarField_X = (NY * FFZ - NZ * FFY) / FarFieldDenominator;
        double const FarField_Y = (NZ * FFX - NX * FFZ) / FarFieldDenominator;
        double const FarField_Z = (NX * FFY - NY * FFX) / FarFieldDenominator;


        // Exponent for fourier transformed field
        cuDoubleComplex Exponent = make_cuDoubleComplex(0, -(*Omega) * (_t[i] + D / (*C)));

        cuDoubleComplex X1 = make_cuDoubleComplex(NearField_X + FarField_X, 0);
        cuDoubleComplex Y1 = make_cuDoubleComplex(NearField_Y + FarField_Y, 0);
        cuDoubleComplex Z1 = make_cuDoubleComplex(NearField_Z + FarField_Z, 0);

        cuDoubleComplex MyEXP = cuCexp(Exponent);

        cuDoubleComplex X2 = cuCmul(X1, MyEXP);
        cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
        cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);

        SumEX = cuCadd(SumEX, X2);
        SumEY = cuCadd(SumEY, Y2);
        SumEZ = cuCadd(SumEZ, Z2);

        }
      }
      }

    }

    if (!_done[threadIdx.x]) {
    cuDoubleComplex TSumEX = cuCmul(make_cuDoubleComplex((*C0) * (dt_total), 0), SumEX);
    cuDoubleComplex TSumEY = cuCmul(make_cuDoubleComplex((*C0) * (dt_total), 0), SumEY);
    cuDoubleComplex TSumEZ = cuCmul(make_cuDoubleComplex((*C0) * (dt_total), 0), SumEZ);

    double const EX = (TSumEX.x * TSumEX.x + TSumEX.y * TSumEX.y);
    double const EY = (TSumEY.x * TSumEY.x + TSumEY.y * TSumEY.y);
    double const EZ = (TSumEZ.x * TSumEZ.x + TSumEZ.y * TSumEZ.y);

    this_result = fabs((*C2) * (EX + EY + EZ));

    if (!_done[threadIdx.x] && (ilevel > 8) && (fabs((this_result - last_result) / last_result) < *prec) ) {
      _done[threadIdx.x] = true;
      result = this_result;
    }

    last_result = this_result;
    }


    __syncthreads();
    if (threadIdx.x == 1) {
      for (int ith = 0; ith < NTHREADS_PER_BLOCK; ++ith) {
        _all_done = true;
        if (!_done[ith]) {
          _all_done = false;
        }
      }
    }

    this_nt *= 2;

    __syncthreads();
  }

  if (is >= *ns) {
    return;
  }

  flux[ith] = (double) is; //result;

  return;
}








__global__ void OSCARSSR_Cuda_FluxGPUMultiWithA (double  *x, double  *y, double  *z,
                                                 double *bx, double *by, double *bz,
                                                 double *ax, double *ay, double *az,
                                                 double *sx, double *sy, double *sz,
                                                 double *dt,
                                                 int *nt,
                                                 int *ns,
                                                 double *C0, double *C2, double *C,
                                                 double *Omega,
                                                 int *ifirst,
                                                 double *flux)
{
  // Check that this is within the number of spectrum points requested
  int const ith = threadIdx.x + blockIdx.x * blockDim.x;
  int const is = ith + *ifirst;
  if (is >= *ns) {
    return;
  }


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

    // Magnitude of Beta squared
    double const One_Minus_BMag2 = 1. -  (bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i]);

    // N dot Beta
    double const NDotBeta = NX * bx[i] + NY * by[i] + NZ * bz[i];

    double const FarFieldDenominator =  D * (pow(1. - NDotBeta, 2));
    double const NearFieldDenominator = D * FarFieldDenominator;
    double const NearField_X = One_Minus_BMag2 * (NX - bx[i]) / NearFieldDenominator;
    double const NearField_Y = One_Minus_BMag2 * (NY - by[i]) / NearFieldDenominator;
    double const NearField_Z = One_Minus_BMag2 * (NZ - bz[i]) / NearFieldDenominator;

    double const FFX = (NY - by[i]) * az[i] - (NZ - bz[i]) * ay[i];
    double const FFY = (NZ - bz[i]) * ax[i] - (NX - bx[i]) * az[i];
    double const FFZ = (NX - bx[i]) * ay[i] - (NY - by[i]) * ax[i];

    double const FarField_X = (NY * FFZ - NZ * FFY) / FarFieldDenominator;
    double const FarField_Y = (NZ * FFX - NX * FFZ) / FarFieldDenominator;
    double const FarField_Z = (NX * FFY - NY * FFX) / FarFieldDenominator;
    

    // Exponent for fourier transformed field
    cuDoubleComplex Exponent = make_cuDoubleComplex(0, -(*Omega) * ((*dt) * ((double) i) + D / (*C)));

    cuDoubleComplex X1 = make_cuDoubleComplex(NearField_X + FarField_X, 0);
    cuDoubleComplex Y1 = make_cuDoubleComplex(NearField_Y + FarField_Y, 0);
    cuDoubleComplex Z1 = make_cuDoubleComplex(NearField_Z + FarField_Z, 0);




    cuDoubleComplex MyEXP = cuCexp(Exponent);

    cuDoubleComplex X2 = cuCmul(X1, MyEXP);
    cuDoubleComplex Y2 = cuCmul(Y1, MyEXP);
    cuDoubleComplex Z2 = cuCmul(Z1, MyEXP);


    SumEX = cuCadd(SumEX, X2);
    SumEY = cuCadd(SumEY, Y2);
    SumEZ = cuCadd(SumEZ, Z2);

  }

  SumEX = cuCmul(make_cuDoubleComplex((*C0) * (*dt), 0), SumEX);
  SumEY = cuCmul(make_cuDoubleComplex((*C0) * (*dt), 0), SumEY);
  SumEZ = cuCmul(make_cuDoubleComplex((*C0) * (*dt), 0), SumEZ);


  double const EX = (SumEX.x * SumEX.x + SumEX.y * SumEX.y);
  double const EY = (SumEY.x * SumEY.x + SumEY.y * SumEY.y);
  double const EZ = (SumEZ.x * SumEZ.x + SumEZ.y * SumEZ.y);

  flux[ith] = (*C2) * (EX + EY + EZ);

  return;
}








extern "C" void OSCARSSR_Cuda_CalculateFluxGPU (OSCARSSR& OSR,
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


  int *h_nt, *h_nt_max, *h_ns;
  double *h_dt;
  cudaHostAlloc((void**) &h_nt_max, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nt,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ns,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_dt,     sizeof(double), cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt_max = (int) OSR.GetNPointsTrajectory();
  *h_nt     = (int) OSR.GetTrajectory().GetNPoints();
  *h_ns     = (int) Surface.GetNPoints();
  *h_dt     = (double) OSR.GetTrajectory().GetDeltaT();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = NTHREADS_PER_BLOCK;
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
  cudaHostAlloc((void**) &h_x,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
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
    cudaMalloc((void **) &d_nt[i],                 sizeof(int));
    cudaMalloc((void **) &d_ns[i],                 sizeof(int));
    cudaMalloc((void **) &d_dt[i],                 sizeof(double));
    cudaMalloc((void **) &d_x[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_y[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_z[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bx[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_by[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bz[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_sx[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sy[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sz[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_c0[i],                 sizeof(double));
    cudaMalloc((void **) &d_c2[i],                 sizeof(double));
    cudaMalloc((void **) &d_c[i],                  sizeof(double));
    cudaMalloc((void **) &d_omega[i],              sizeof(double));
    cudaMalloc((void **) &d_ifirst[i],             sizeof(int));
    cudaMalloc((void **) &d_flux[i],       NFlux * sizeof(double));

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
  int const NPointsThisTrajectory = T.GetNPoints();
  *h_nt = 0;
  for (size_t i = 0; i < NPointsThisTrajectory; ++i) {
    if (T.GetA(i).Mag() < 1e-100) {
      continue;
    }
    h_x[*h_nt]  = T.GetX(i).GetX();
    h_y[*h_nt]  = T.GetX(i).GetY();
    h_z[*h_nt]  = T.GetX(i).GetZ();
    h_bx[*h_nt] = T.GetB(i).GetX();
    h_by[*h_nt] = T.GetB(i).GetY();
    h_bz[*h_nt] = T.GetB(i).GetZ();
    ++(*h_nt);
  }
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_nt[0],    h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_nt[i+1],     d1, d_nt[i],     d, sizeof(int));
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
    cudaMemcpyAsync(d_nt[0], h_nt,         sizeof(int),    cudaMemcpyHostToDevice);
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
      cudaMemcpyPeerAsync(d_nt[ig+1], d1, d_nt[ig], d,         sizeof(int));
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
      int const NPointsThisTrajectory = T.GetNPoints();

      *h_nt = 0;
      for (size_t it = 0; it < NPointsThisTrajectory; ++it) {
        if (T.GetA(it).Mag() < 1e-100) {
          continue;
        }
        h_x[*h_nt]  = T.GetX(it).GetX();
        h_y[*h_nt]  = T.GetX(it).GetY();
        h_z[*h_nt]  = T.GetX(it).GetZ();
        h_bx[*h_nt] = T.GetB(it).GetX();
        h_by[*h_nt] = T.GetB(it).GetY();
        h_bz[*h_nt] = T.GetB(it).GetZ();
        ++(*h_nt);
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
  cudaFreeHost(h_nt_max);
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

















extern "C" void OSCARSSR_Cuda_CalculateFluxGPUWithA (OSCARSSR& OSR,
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


  int *h_nt, *h_nt_max, *h_ns;
  double *h_dt;
  cudaHostAlloc((void**) &h_nt_max, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nt,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ns,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_dt,     sizeof(double), cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt_max = (int) OSR.GetNPointsTrajectory();
  *h_nt     = (int) OSR.GetTrajectory().GetNPoints();
  *h_ns     = (int) Surface.GetNPoints();
  *h_dt     = (double) OSR.GetTrajectory().GetDeltaT();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = NTHREADS_PER_BLOCK;
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
  double  *h_x,  *h_y,  *h_z,  *h_bx,  *h_by,  *h_bz,  *h_ax, *h_ay, *h_az, *h_sx,  *h_sy,  *h_sz,   *h_c0,  *h_c2,  *h_c,  *h_omega;
  int     *h_ifirst;
  double **h_flux;
  cudaHostAlloc((void**) &h_x,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ax,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ay,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_az,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
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
  double **d_ax;
  double **d_ay;
  double **d_az;
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
  cudaHostAlloc((void **) &d_ax,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ay,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_az,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
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
    cudaMalloc((void **) &d_nt[i],                 sizeof(int));
    cudaMalloc((void **) &d_ns[i],                 sizeof(int));
    cudaMalloc((void **) &d_dt[i],                 sizeof(double));
    cudaMalloc((void **) &d_x[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_y[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_z[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bx[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_by[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bz[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_ax[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_ay[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_az[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_sx[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sy[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sz[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_c0[i],                 sizeof(double));
    cudaMalloc((void **) &d_c2[i],                 sizeof(double));
    cudaMalloc((void **) &d_c[i],                  sizeof(double));
    cudaMalloc((void **) &d_omega[i],              sizeof(double));
    cudaMalloc((void **) &d_ifirst[i],             sizeof(int));
    cudaMalloc((void **) &d_flux[i],       NFlux * sizeof(double));

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
  int const NPointsThisTrajectory = T.GetNPoints();
  *h_nt = 0;
  for (size_t i = 0; i < NPointsThisTrajectory; ++i) {
    h_x[*h_nt]  = T.GetX(i).GetX();
    h_y[*h_nt]  = T.GetX(i).GetY();
    h_z[*h_nt]  = T.GetX(i).GetZ();
    h_bx[*h_nt] = T.GetB(i).GetX();
    h_by[*h_nt] = T.GetB(i).GetY();
    h_bz[*h_nt] = T.GetB(i).GetZ();
    h_ax[*h_nt] = T.GetAoverC(i).GetX();
    h_ay[*h_nt] = T.GetAoverC(i).GetY();
    h_az[*h_nt] = T.GetAoverC(i).GetZ();
    ++(*h_nt);
  }
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_nt[0],    h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_nt[i+1],     d1, d_nt[i],     d, sizeof(int));
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
    cudaMemcpyAsync(d_nt[0], h_nt,         sizeof(int),    cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_x[0],  h_x,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_y[0],  h_y,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_z[0],  h_z,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bx[0], h_bx, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_by[0], h_by, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bz[0], h_bz, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ax[0], h_ax, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ay[0], h_ay, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_az[0], h_az, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    for (size_t ig = 0; ig < GPUsToUse.size() - 1; ++ig) {
      // Device number
      int const d  = GPUsToUse[ig];
      int const d1 = GPUsToUse[ig+1];
      cudaSetDevice(d);
      cudaMemcpyPeerAsync(d_nt[ig+1], d1, d_nt[ig], d,         sizeof(int));
      cudaMemcpyPeerAsync(d_x[ig+1],  d1, d_x[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_y[ig+1],  d1, d_y[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_z[ig+1],  d1, d_z[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bx[ig+1], d1, d_bx[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_by[ig+1], d1, d_by[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bz[ig+1], d1, d_bz[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ax[ig+1], d1, d_ax[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ay[ig+1], d1, d_ay[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_az[ig+1], d1, d_az[ig], d, *h_nt * sizeof(double));

    }

    // Wait for previous copy, start next one
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaEventSynchronize(event_fluxcopy[ig]);
      OSCARSSR_Cuda_FluxGPUMultiWithA<<<NBlocksThisGPU[ig], NThreadsPerBlock>>>( d_x[ig],  d_y[ig],  d_z[ig],
                                                                           d_bx[ig], d_by[ig], d_bz[ig],
                                                                           d_ax[ig], d_ay[ig], d_az[ig],
                                                                           d_sx[ig], d_sy[ig], d_sz[ig],
                                                                           d_dt[ig],
                                                                           d_nt[ig],
                                                                           d_ns[ig],
                                                                           d_c0[ig], d_c2[ig], d_c[ig],
                                                                           d_omega[ig],
                                                                           d_ifirst[ig],
                                                                           d_flux[ig]);
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
      int const NPointsThisTrajectory = T.GetNPoints();

      *h_nt = 0;
      for (size_t it = 0; it < NPointsThisTrajectory; ++it) {
        h_x[*h_nt]  = T.GetX(it).GetX();
        h_y[*h_nt]  = T.GetX(it).GetY();
        h_z[*h_nt]  = T.GetX(it).GetZ();
        h_bx[*h_nt] = T.GetB(it).GetX();
        h_by[*h_nt] = T.GetB(it).GetY();
        h_bz[*h_nt] = T.GetB(it).GetZ();
        h_ax[*h_nt] = T.GetAoverC(it).GetX();
        h_ay[*h_nt] = T.GetAoverC(it).GetY();
        h_az[*h_nt] = T.GetAoverC(it).GetZ();
        ++(*h_nt);
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
  cudaFreeHost(h_nt_max);
  cudaFreeHost(h_nt);
  cudaFreeHost(h_ns);
  cudaFreeHost(h_dt);
  cudaFreeHost(h_x);
  cudaFreeHost(h_y);
  cudaFreeHost(h_z);
  cudaFreeHost(h_bx);
  cudaFreeHost(h_by);
  cudaFreeHost(h_bz);
  cudaFreeHost(h_ax);
  cudaFreeHost(h_ay);
  cudaFreeHost(h_az);
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
    cudaFree(d_ax[i]);
    cudaFree(d_ay[i]);
    cudaFree(d_az[i]);
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
  cudaFree(d_ax);
  cudaFree(d_ay);
  cudaFree(d_az);
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
























extern "C" void OSCARSSR_Cuda_CalculateFluxGPUWithAInterpolated (OSCARSSR& OSR,
                                                TSurfacePoints const& Surface,
                                                double const Energy_eV,
                                                T3DScalarContainer& FluxContainer,
                                                std::string const& Polarization,
                                                double const Angle,
                                                TVector3D const& HorizontalDirection,
                                                TVector3D const& PropogationDirection,
                                                int const NParticles,
                                                std::vector<int> const& GPUVector,
                                                double const Precision,
                                                int const MaxLevel)
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

  // nt - number of interpolated track points
  // t  - array of timestamps in trajectory
  // x, y, z; xp, yp, zp - position and derivs
  // beta x, y, z; xp, yp, zp
  // aOc x, y, z; xp, yp, zp
  // max_level_extended or max_level
  // surface points
  // nthreads_per_block

  int *h_ns;
  int *h_nt;
  cudaHostAlloc((void**) &h_ns,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nt,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt     = (int) OSR.GetCurrentParticle().GetTrajectoryInterpolated().GetNPoints();
  *h_ns     = (int) Surface.GetNPoints();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = NTHREADS_PER_BLOCK;
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

  // Interpolating structure
  double  *h_t;
  double  *h_tstart, *h_tstop;
  double  *h_x,   *h_y,   *h_z;
  double  *h_xp,  *h_yp,  *h_zp;
  double  *h_bx,  *h_by,  *h_bz;
  double  *h_bxp, *h_byp, *h_bzp;
  double  *h_ax,  *h_ay,  *h_az;
  double  *h_axp, *h_ayp, *h_azp;

  // Surface points
  double  *h_sx,  *h_sy,  *h_sz;

  // Constants and photon frequency
  double *h_c0,  *h_c2,  *h_c,  *h_omega;

  // first point for each thread, max level
  int     *h_ifirst;
  int     *h_ml;

  double *h_prec;

  // Results
  double **h_flux;

  // Allocate host memory
  cudaHostAlloc((void**) &h_t,       *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_tstart,          sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_tstop,           sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_x,       *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,       *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,       *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_xp,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_yp,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_zp,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bxp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_byp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bzp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ax,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ay,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_az,      *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_axp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ayp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_azp,     *h_nt * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_sx,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sy,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sz,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_c0,                  sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c2,                  sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c,                   sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_omega,               sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ifirst, NGPUsToUse * sizeof(int),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ml,                  sizeof(int),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_prec,                sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_flux,   NGPUsToUse * sizeof(double*), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    cudaHostAlloc((void**) &(h_flux[i]), NFlux * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  }

  // First surface point for each gpu
  int NBlocksUsed = 0;
  for (int i = 0; i < NGPUsToUse; ++i) {
    h_ifirst[i] = NBlocksUsed * NThreadsPerBlock;
    NBlocksUsed += NBlocksThisGPU[i];
  }

  *h_ml = MaxLevel; //UPDATE: max level should be an input

  // Precision
  *h_prec = Precision;

  // Memor allocations for GPU
  int    **d_nt;
  double **d_tstart;
  double **d_tstop;
  int    **d_ns;

  double **d_t;

  double **d_x;
  double **d_y;
  double **d_z;
  double **d_xp;
  double **d_yp;
  double **d_zp;

  double **d_bx;
  double **d_by;
  double **d_bz;
  double **d_bxp;
  double **d_byp;
  double **d_bzp;

  double **d_ax;
  double **d_ay;
  double **d_az;
  double **d_axp;
  double **d_ayp;
  double **d_azp;

  double **d_sx;
  double **d_sy;
  double **d_sz;

  double **d_c0;
  double **d_c2;
  double **d_c;
  double **d_omega;

  int    **d_ifirst;
  int    **d_ml;
  double **d_prec;
  double **d_flux;

  cudaHostAlloc((void **) &d_nt,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_tstart, NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_tstop,  NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ns,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_t,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_x,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_y,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_z,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_xp,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_yp,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_zp,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_bx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_by,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bxp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_byp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bzp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_ax,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ay,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_az,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_axp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ayp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_azp,    NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_sx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sy,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_c0,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c2,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_omega,  NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void **) &d_ifirst, NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ml,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_prec,   NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_flux,   NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaMalloc((void **) &d_nt[i],                 sizeof(int));
    cudaMalloc((void **) &d_tstart[i],             sizeof(double));
    cudaMalloc((void **) &d_tstop[i],              sizeof(double));
    cudaMalloc((void **) &d_ns[i],                 sizeof(int));

    cudaMalloc((void **) &d_t[i],          *h_nt * sizeof(double));

    cudaMalloc((void **) &d_x[i],          *h_nt * sizeof(double));
    cudaMalloc((void **) &d_y[i],          *h_nt * sizeof(double));
    cudaMalloc((void **) &d_z[i],          *h_nt * sizeof(double));
    cudaMalloc((void **) &d_xp[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_yp[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_zp[i],         *h_nt * sizeof(double));

    cudaMalloc((void **) &d_bx[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_by[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_bz[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_bxp[i],        *h_nt * sizeof(double));
    cudaMalloc((void **) &d_byp[i],        *h_nt * sizeof(double));
    cudaMalloc((void **) &d_bzp[i],        *h_nt * sizeof(double));

    cudaMalloc((void **) &d_ax[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_ay[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_az[i],         *h_nt * sizeof(double));
    cudaMalloc((void **) &d_axp[i],        *h_nt * sizeof(double));
    cudaMalloc((void **) &d_ayp[i],        *h_nt * sizeof(double));
    cudaMalloc((void **) &d_azp[i],        *h_nt * sizeof(double));

    cudaMalloc((void **) &d_sx[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sy[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sz[i],         *h_ns * sizeof(double));

    cudaMalloc((void **) &d_c0[i],                 sizeof(double));
    cudaMalloc((void **) &d_c2[i],                 sizeof(double));
    cudaMalloc((void **) &d_c[i],                  sizeof(double));
    cudaMalloc((void **) &d_omega[i],              sizeof(double));

    cudaMalloc((void **) &d_ifirst[i],             sizeof(int));
    cudaMalloc((void **) &d_ml[i],                 sizeof(int));
    cudaMalloc((void **) &d_prec[i],               sizeof(double));
    cudaMalloc((void **) &d_flux[i],       NFlux * sizeof(double));

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
  cudaMemcpyAsync(d_ns[0],    h_ns,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c0[0],    h_c0,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c2[0],    h_c2,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c[0],     h_c,           sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_omega[0], h_omega,       sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sx[0],    h_sx,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sy[0],    h_sy,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sz[0],    h_sz,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_ml[0],    h_ml,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_prec[0],  h_prec,        sizeof(double), cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_ns[i+1],     d1, d_ns[i],     d, sizeof(int));
    cudaMemcpyPeerAsync( d_c0[i+1],     d1, d_c0[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_c2[i+1],     d1, d_c2[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_c[i+1],      d1, d_c[i],      d, sizeof(double));
    cudaMemcpyPeerAsync( d_omega[i+1],  d1, d_omega[i],  d, sizeof(double));
    cudaMemcpyPeerAsync( d_sx[i+1],     d1, d_sx[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sy[i+1],     d1, d_sy[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sz[i+1],     d1, d_sz[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_ml[i+1],     d1, d_ml[i],     d, sizeof(int));
    cudaMemcpyPeerAsync( d_prec[i+1],   d1, d_prec[i],   d, sizeof(double));
  }

  // Set first trajectory
  TOMATH::TSpline1D3<TParticleTrajectoryPoint> const& T = OSR.GetCurrentParticle().GetTrajectoryInterpolated().GetSpline();
  int const NPointsThisTrajectory = T.GetNPoints();
  *h_nt = NPointsThisTrajectory;
  *h_tstart = T.GetXStart();
  *h_tstop  = T.GetXStop();
  for (size_t i = 0; i < NPointsThisTrajectory; ++i) {
    TParticleTrajectoryPoint const& P  = T.GetY(i);
    TParticleTrajectoryPoint const& PP = T.GetYPP(i);

    h_t[i]   = T.GetX(i);

    h_x[i]   =  P.GetX().GetX();
    h_y[i]   =  P.GetX().GetY();
    h_z[i]   =  P.GetX().GetZ();
    h_xp[i]  = PP.GetX().GetX();
    h_yp[i]  = PP.GetX().GetY();
    h_zp[i]  = PP.GetX().GetZ();

    h_bx[i]  =  P.GetB().GetX();
    h_by[i]  =  P.GetB().GetY();
    h_bz[i]  =  P.GetB().GetZ();
    h_bxp[i] = PP.GetB().GetX();
    h_byp[i] = PP.GetB().GetY();
    h_bzp[i] = PP.GetB().GetZ();

    h_ax[i]  =  P.GetAoverC().GetX();
    h_ay[i]  =  P.GetAoverC().GetY();
    h_az[i]  =  P.GetAoverC().GetZ();
    h_axp[i] = PP.GetAoverC().GetX();
    h_ayp[i] = PP.GetAoverC().GetY();
    h_azp[i] = PP.GetAoverC().GetZ();
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
    cudaMemcpyAsync(d_nt[0],  h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_t[0],   h_t,   *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_tstart[0], h_tstart,   sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_tstop[0],  h_tstop,    sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_x[0],   h_x,   *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_y[0],   h_y,   *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_z[0],   h_z,   *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_xp[0],  h_xp,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_yp[0],  h_yp,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_zp[0],  h_zp,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bx[0],  h_bx,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_by[0],  h_by,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bz[0],  h_bz,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bxp[0], h_bxp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_byp[0], h_byp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bzp[0], h_bzp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ax[0],  h_ax,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ay[0],  h_ay,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_az[0],  h_az,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_axp[0], h_axp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ayp[0], h_ayp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_azp[0], h_azp, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    for (size_t ig = 0; ig < GPUsToUse.size() - 1; ++ig) {
      // Device number
      int const d  = GPUsToUse[ig];
      int const d1 = GPUsToUse[ig+1];
      cudaSetDevice(d);
      cudaMemcpyPeerAsync(d_nt[ig+1],  d1, d_nt[ig],  d,         sizeof(int));
      cudaMemcpyPeerAsync(d_t[ig+1],   d1, d_t[ig],   d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_tstart[ig+1], d1, d_tstart[ig], d,   sizeof(double));
      cudaMemcpyPeerAsync(d_tstop[ig+1],  d1, d_tstop[ig],  d,    sizeof(double));
      cudaMemcpyPeerAsync(d_x[ig+1],   d1, d_x[ig],   d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_y[ig+1],   d1, d_y[ig],   d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_z[ig+1],   d1, d_z[ig],   d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_xp[ig+1],  d1, d_xp[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_yp[ig+1],  d1, d_yp[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_zp[ig+1],  d1, d_zp[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bx[ig+1],  d1, d_bx[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_by[ig+1],  d1, d_by[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bz[ig+1],  d1, d_bz[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bxp[ig+1], d1, d_bxp[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_byp[ig+1], d1, d_byp[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bzp[ig+1], d1, d_bzp[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ax[ig+1],  d1, d_ax[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ay[ig+1],  d1, d_ay[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_az[ig+1],  d1, d_az[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_axp[ig+1], d1, d_axp[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ayp[ig+1], d1, d_ayp[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_azp[ig+1], d1, d_azp[ig], d, *h_nt * sizeof(double));

    }

    // Wait for previous copy, start next one
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaEventSynchronize(event_fluxcopy[ig]);
      OSCARSSR_Cuda_FluxGPUMultiWithAInterpolated<<<NBlocksThisGPU[ig], NThreadsPerBlock>>>(
                                                                           d_t[ig],
                                                                           d_x[ig],   d_y[ig],   d_z[ig],
                                                                           d_xp[ig],  d_yp[ig],  d_zp[ig],
                                                                           d_bx[ig],  d_by[ig],  d_bz[ig],
                                                                           d_bxp[ig], d_byp[ig], d_bzp[ig],
                                                                           d_ax[ig],  d_ay[ig],  d_az[ig],
                                                                           d_axp[ig], d_ayp[ig], d_azp[ig],
                                                                           d_sx[ig],  d_sy[ig],  d_sz[ig],
                                                                           d_tstart[ig], d_tstop[ig],
                                                                           d_nt[ig],
                                                                           d_ns[ig],
                                                                           d_c0[ig], d_c2[ig], d_c[ig],
                                                                           d_omega[ig],
                                                                           d_ifirst[ig],
                                                                           d_ml[ig],
                                                                           d_prec[ig],
                                                                           d_flux[ig]);
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

      TOMATH::TSpline1D3<TParticleTrajectoryPoint> const& T = OSR.GetCurrentParticle().GetTrajectoryInterpolated().GetSpline();
      int const NPointsThisTrajectory = T.GetNPoints();
      *h_nt = NPointsThisTrajectory;
      for (size_t it = 0; it < NPointsThisTrajectory; ++it) {
        TParticleTrajectoryPoint const& P  = T.GetY(it);
        TParticleTrajectoryPoint const& PP = T.GetYPP(it);

        h_t[it]   = T.GetX(it);

        h_x[it]   =  P.GetX().GetX();
        h_y[it]   =  P.GetX().GetY();
        h_z[it]   =  P.GetX().GetZ();
        h_xp[it]  = PP.GetX().GetX();
        h_yp[it]  = PP.GetX().GetY();
        h_zp[it]  = PP.GetX().GetZ();

        h_bx[it]  =  P.GetB().GetX();
        h_by[it]  =  P.GetB().GetY();
        h_bz[it]  =  P.GetB().GetZ();
        h_bxp[it] = PP.GetB().GetX();
        h_byp[it] = PP.GetB().GetY();
        h_bzp[it] = PP.GetB().GetZ();

        h_ax[it]  =  P.GetAoverC().GetX();
        h_ay[it]  =  P.GetAoverC().GetY();
        h_az[it]  =  P.GetAoverC().GetZ();
        h_axp[it] = PP.GetAoverC().GetX();
        h_ayp[it] = PP.GetAoverC().GetY();
        h_azp[it] = PP.GetAoverC().GetZ();
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
  cudaFreeHost(h_t);
  cudaFreeHost(h_ns);
  cudaFreeHost(h_x);
  cudaFreeHost(h_y);
  cudaFreeHost(h_z);
  cudaFreeHost(h_xp);
  cudaFreeHost(h_yp);
  cudaFreeHost(h_zp);
  cudaFreeHost(h_bx);
  cudaFreeHost(h_by);
  cudaFreeHost(h_bz);
  cudaFreeHost(h_bxp);
  cudaFreeHost(h_byp);
  cudaFreeHost(h_bzp);
  cudaFreeHost(h_ax);
  cudaFreeHost(h_ay);
  cudaFreeHost(h_az);
  cudaFreeHost(h_axp);
  cudaFreeHost(h_ayp);
  cudaFreeHost(h_azp);
  cudaFreeHost(h_sx);
  cudaFreeHost(h_sy);
  cudaFreeHost(h_sz);
  cudaFreeHost(h_c0);
  cudaFreeHost(h_c2);
  cudaFreeHost(h_c);
  cudaFreeHost(h_omega);
  cudaFreeHost(h_ifirst);
  cudaFreeHost(h_ml);
  cudaFreeHost(h_prec);
  // Free host and GPU memory
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {

    cudaFreeHost(h_flux[i]);

    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaFree(d_nt[i]);
    cudaFree(d_t[i]);
    cudaFree(d_ns[i]);
    cudaFree(d_x[i]);
    cudaFree(d_y[i]);
    cudaFree(d_z[i]);
    cudaFree(d_xp[i]);
    cudaFree(d_yp[i]);
    cudaFree(d_zp[i]);
    cudaFree(d_bx[i]);
    cudaFree(d_by[i]);
    cudaFree(d_bz[i]);
    cudaFree(d_bxp[i]);
    cudaFree(d_byp[i]);
    cudaFree(d_bzp[i]);
    cudaFree(d_ax[i]);
    cudaFree(d_ay[i]);
    cudaFree(d_az[i]);
    cudaFree(d_axp[i]);
    cudaFree(d_ayp[i]);
    cudaFree(d_azp[i]);
    cudaFree(d_sx[i]);
    cudaFree(d_sy[i]);
    cudaFree(d_sz[i]);
    cudaFree(d_c0[i]);
    cudaFree(d_c2[i]);
    cudaFree(d_c[i]);
    cudaFree(d_omega[i]);
    cudaFree(d_ifirst[i]);
    cudaFree(d_ml[i]);
    cudaFree(d_prec[i]);
    cudaFree(d_flux[i]);
  }
  cudaFree(h_flux);

  cudaFree(d_nt);
  cudaFree(d_t);
  cudaFree(d_ns);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_xp);
  cudaFree(d_yp);
  cudaFree(d_zp);
  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);
  cudaFree(d_bxp);
  cudaFree(d_byp);
  cudaFree(d_bzp);
  cudaFree(d_ax);
  cudaFree(d_ay);
  cudaFree(d_az);
  cudaFree(d_axp);
  cudaFree(d_ayp);
  cudaFree(d_azp);
  cudaFree(d_sx);
  cudaFree(d_sy);
  cudaFree(d_sz);
  cudaFree(d_c0);
  cudaFree(d_c2);
  cudaFree(d_c);
  cudaFree(d_omega);
  cudaFree(h_ifirst);
  cudaFree(d_ml);
  cudaFree(d_prec);
  cudaFree(d_flux);

  // Delete host gpu pointer arrays
  delete [] event_fluxcopy;


  return;
}




























































__global__ void OSCARSSR_Cuda_SpectrumGPUMulti (double *x, double *y, double *z, double *bx, double *by, double *bz, double *obs, double *dt, int *nt, int *ns, double *C0, double *C2, double *EvToOmega, double *C, double *se, double *sf, cuDoubleComplex* pol, int *pol_state, int const *ifirst)
{
  // Check that this is within the number of spectrum points requested
  int const ith = threadIdx.x + blockIdx.x * blockDim.x;
  int const is = ith + *ifirst;
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
    double const D = sqrt( pow( obs[0] - x[i], 2) + pow( obs[1] - y[i], 2) + pow(obs[2] - z[i], 2) );

    // Normal in direction of observer
    double const NX = (obs[0] - x[i]) / D;
    double const NY = (obs[1] - y[i]) / D;
    double const NZ = (obs[2] - z[i]) / D;

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

  sf[ith] = (*C2) * (EX + EY + EZ);

  return;
}








extern "C" void OSCARSSR_Cuda_CalculateSpectrumGPU (OSCARSSR& OSR,
                                                    TParticleA& Particle,
                                                    TVector3D const& ObservationPoint,
                                                    TSpectrumContainer& Spectrum,
                                                    std::string const& Polarization,
                                                    double const Angle,
                                                    TVector3D const& HorizontalDirection,
                                                    TVector3D const& PropogationDirection,
                                                    int const NParticles,
                                                    std::vector<int> const& GPUVector)
{
  // Calculate the spectrum for NParticles using the GPUs given in GPUVector.  Each particle's
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

  int *h_nt, *h_nt_max, *h_ns;
  double *h_dt;
  cudaHostAlloc((void**) &h_nt_max, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nt,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ns,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_dt,     sizeof(double), cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt_max = (int) OSR.GetNPointsTrajectory();
  *h_nt     = (int) OSR.GetTrajectory().GetNPoints();
  *h_ns     = (int) Spectrum.GetNPoints();
  *h_dt     = (double) OSR.GetTrajectory().GetDeltaT();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = NTHREADS_PER_BLOCK;
  int const NThreadsRemainder = NThreads % NThreadsPerBlock;
  int const NBlocksTotal = (NThreads - 1) / NThreadsPerBlock + 1;
  int const NBlocksPerGPU = NBlocksTotal / NGPUsToUse;
  int const NRemainderBlocks = NBlocksTotal % NGPUsToUse;
  // UPDATE: To be modified
  int const NSpectrum = NThreadsPerBlock * (NBlocksPerGPU + (NRemainderBlocks > 0 ? 1 : 0));

  std::vector<int> NBlocksThisGPU(NGPUsToUse, NBlocksPerGPU);
  for (int i = 0; i < NRemainderBlocks; ++i) {
    ++NBlocksThisGPU[i];
  }

  // Memory allocation for Host
  double  *h_x,  *h_y,  *h_z,  *h_bx,  *h_by,  *h_bz,  *h_obs, *h_se,   *h_c0,  *h_c2,  *h_c,  *h_ev2omega;
  int     *h_ifirst, *h_pol_state;
  double **h_spectrum;
  cuDoubleComplex *h_pol;

  cudaHostAlloc((void**) &h_x,       *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,       *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,       *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,      *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,      *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,      *h_nt_max * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_obs,             3 * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_se,          *h_ns * sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c0,                  sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c2,                  sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_c,                   sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ev2omega,            sizeof(double),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_pol,                 sizeof(int),             cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_pol_state,       3 * sizeof(cuDoubleComplex), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ifirst, NGPUsToUse * sizeof(int),             cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_spectrum, NGPUsToUse * sizeof(double*),       cudaHostAllocWriteCombined | cudaHostAllocMapped);
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    cudaHostAlloc((void**) &(h_spectrum[i]), NSpectrum * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  }

  // First surface point for each gpu
  int NBlocksUsed = 0;
  for (int i = 0; i < NGPUsToUse; ++i) {
    h_ifirst[i] = NBlocksUsed * NThreadsPerBlock;
    NBlocksUsed += NBlocksThisGPU[i];
  }


  // Memor allocations for GPU
  int             **d_nt;
  int             **d_ns;
  double          **d_dt;
  double          **d_x;
  double          **d_y;
  double          **d_z;
  double          **d_bx;
  double          **d_by;
  double          **d_bz;
  double          **d_obs;
  double          **d_se;
  double          **d_c0;
  double          **d_c2;
  double          **d_c;
  double          **d_ev2omega;
  int             **d_ifirst;
  cuDoubleComplex **d_pol;
  int             **d_pol_state;
  double          **d_spectrum;

  cudaHostAlloc((void **) &d_nt,        NGPUsToUse * sizeof(int*),             cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ns,        NGPUsToUse * sizeof(int*),             cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_dt,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_x,         NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_y,         NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_z,         NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bx,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_by,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bz,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_obs,       NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_se,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c0,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c2,        NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_c,         NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ev2omega,  NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ifirst,    NGPUsToUse * sizeof(int*),             cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_pol,       NGPUsToUse * sizeof(cuDoubleComplex*), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_pol_state, NGPUsToUse * sizeof(int*),             cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_spectrum,  NGPUsToUse * sizeof(double*),          cudaHostAllocWriteCombined | cudaHostAllocMapped);

  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaMalloc((void **) &d_nt[i],                   sizeof(int));
    cudaMalloc((void **) &d_ns[i],                   sizeof(int));
    cudaMalloc((void **) &d_dt[i],                   sizeof(double));
    cudaMalloc((void **) &d_x[i],        *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_y[i],        *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_z[i],        *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bx[i],       *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_by[i],       *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bz[i],       *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_obs[i],              3 * sizeof(double));
    cudaMalloc((void **) &d_se[i],           *h_ns * sizeof(double));
    cudaMalloc((void **) &d_c0[i],                   sizeof(double));
    cudaMalloc((void **) &d_c2[i],                   sizeof(double));
    cudaMalloc((void **) &d_c[i],                    sizeof(double));
    cudaMalloc((void **) &d_ev2omega[i],             sizeof(double));
    cudaMalloc((void **) &d_ifirst[i],               sizeof(int));
    cudaMalloc((void **) &d_pol[i],              3 * sizeof(cuDoubleComplex));
    cudaMalloc((void **) &d_pol_state[i],            sizeof(int));
    cudaMalloc((void **) &d_spectrum[i], NSpectrum * sizeof(double));

    // Copy device number to device
    cudaMemcpyAsync(d_ifirst[i], &(h_ifirst[i]), sizeof(int), cudaMemcpyHostToDevice);
  }

  // Compute known host values
  h_obs[0]     = ObservationPoint[0];
  h_obs[1]     = ObservationPoint[1];
  h_obs[2]     = ObservationPoint[2];
  *h_c0        = OSR.GetCurrentParticle().GetQ() / (TOSCARSSR::FourPi() * TOSCARSSR::C() * TOSCARSSR::Epsilon0() * TOSCARSSR::Sqrt2Pi());
  *h_c2        = TOSCARSSR::FourPi() * OSR.GetCurrentParticle().GetCurrent() / (TOSCARSSR::H() * fabs(OSR.GetCurrentParticle().GetQ()) * TOSCARSSR::Mu0() * TOSCARSSR::C()) * 1e-6 * 0.001;
  *h_c         = TOSCARSSR::C();
  *h_ev2omega  = TOSCARSSR::EvToAngularFrequency(1);
  h_pol[0]     = pol[0];
  h_pol[1]     = pol[1];
  h_pol[2]     = pol[2];
  *h_pol_state = pol_state;
  for (size_t i = 0; i < *h_ns; ++i) {
    h_se[i] = Spectrum.GetEnergy(i);
  }

  // Copy constants to first device (async)
  int const d0 = GPUsToUse[0];
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_nt[0],       h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_ns[0],       h_ns,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_dt[0],       h_dt,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_obs[0],      h_obs,     3 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c0[0],       h_c0,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c2[0],       h_c2,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_c[0],        h_c,           sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_ev2omega[0], h_ev2omega,    sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_se[0],       h_se,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_pol[0],      h_pol,     3 * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_pol_state[0],h_pol_state,   sizeof(int),    cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_nt[i+1],        d1, d_nt[i],       d,         sizeof(int));
    cudaMemcpyPeerAsync( d_ns[i+1],        d1, d_ns[i],       d,         sizeof(int));
    cudaMemcpyPeerAsync( d_dt[i+1],        d1, d_dt[i],       d,         sizeof(double));
    cudaMemcpyPeerAsync( d_obs[i+1],       d1, d_obs[i],      d,     3 * sizeof(double));
    cudaMemcpyPeerAsync( d_c0[i+1],        d1, d_c0[i],       d,         sizeof(double));
    cudaMemcpyPeerAsync( d_c2[i+1],        d1, d_c2[i],       d,         sizeof(double));
    cudaMemcpyPeerAsync( d_c[i+1],         d1, d_c[i],        d,         sizeof(double));
    cudaMemcpyPeerAsync( d_ev2omega[i+1],  d1, d_ev2omega[i], d,         sizeof(double));
    cudaMemcpyPeerAsync( d_se[i+1],        d1, d_se[i],       d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_pol[i+1],       d1, d_pol[i],      d,     3 * sizeof(cuDoubleComplex));
    cudaMemcpyPeerAsync( d_pol_state[i+1], d1, d_pol_state[i],d,         sizeof(int));
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
  cudaEvent_t *event_spectrumcopy = new cudaEvent_t[NGPUsToUse];
  for (int ig = 0; ig < NGPUsToUse; ++ig) {
    int const d = GPUsToUse[ig];
    cudaSetDevice(d);
    cudaEventCreate(&(event_spectrumcopy[ig]));
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
    cudaMemcpyAsync(d_nt[0], h_nt,         sizeof(int),    cudaMemcpyHostToDevice);
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
      cudaMemcpyPeerAsync(d_nt[ig+1], d1, d_nt[ig], d,         sizeof(int));
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
      cudaEventSynchronize(event_spectrumcopy[ig]);
      OSCARSSR_Cuda_SpectrumGPUMulti<<<NBlocksThisGPU[ig], NThreadsPerBlock>>>(d_x[ig], d_y[ig], d_z[ig], d_bx[ig], d_by[ig], d_bz[ig], d_obs[ig], d_dt[ig], d_nt[ig], d_ns[ig], d_c0[ig], d_c2[ig], d_ev2omega[ig], d_c[ig], d_se[ig], d_spectrum[ig], d_pol[ig], d_pol_state[ig], d_ifirst[ig]);
    }


    // Add result to spectrum container (from **previous**)
    if (ip > 0) {
      int NBlocksUsed = 0;
      for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
        for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
          if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
            break;
          }
          int iss = ith + NThreadsPerBlock * NBlocksUsed;
          Spectrum.AddToFlux(iss, h_spectrum[ig][ith]);
        }
        NBlocksUsed += NBlocksThisGPU[ig];
      }
    }

    // Add copy back to streams
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d  = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaMemcpyAsync(h_spectrum[ig],  d_spectrum[ig],  NSpectrum * sizeof(double), cudaMemcpyDeviceToHost);
      cudaEventRecord(event_spectrumcopy[ig]);
    }

    // If it's not the last one, calculate a new trajectory
    if (ip < NParticlesReally - 1) {
      OSR.SetNewParticle();
      OSR.CalculateTrajectory();
      TParticleTrajectoryPoints const& T = OSR.GetTrajectory();
      *h_nt = T.GetNPoints();

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
    cudaEventSynchronize(event_spectrumcopy[ig]);
  }

  // Add result to spectrum container (from **previous**)
  NBlocksUsed = 0;
  for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
    for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
      if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
        break;
      }
      int iss = ith + NThreadsPerBlock * NBlocksUsed;
      Spectrum.AddToFlux(iss, h_spectrum[ig][ith]);
    }
    NBlocksUsed += NBlocksThisGPU[ig];
  }

  // Weighting for multi-particle
  double const Weight = 1.0 / (double) NParticlesReally;
  Spectrum.Scale(Weight);

  // Free host memory
  cudaFreeHost(h_nt_max);
  cudaFreeHost(h_nt);
  cudaFreeHost(h_ns);
  cudaFreeHost(h_dt);
  cudaFreeHost(h_x);
  cudaFreeHost(h_y);
  cudaFreeHost(h_z);
  cudaFreeHost(h_bx);
  cudaFreeHost(h_by);
  cudaFreeHost(h_bz);
  cudaFreeHost(h_obs);
  cudaFreeHost(h_se);
  cudaFreeHost(h_c0);
  cudaFreeHost(h_c2);
  cudaFreeHost(h_c);
  cudaFreeHost(h_ev2omega);
  cudaFreeHost(h_ifirst);
  cudaFreeHost(h_pol);
  cudaFreeHost(h_pol_state);
  // Free host and GPU memory
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {

    cudaFreeHost(h_spectrum[i]);

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
    cudaFree(d_obs[i]);
    cudaFree(d_se[i]);
    cudaFree(d_c0[i]);
    cudaFree(d_c2[i]);
    cudaFree(d_c[i]);
    cudaFree(d_ev2omega[i]);
    cudaFree(d_ifirst[i]);
    cudaFree(d_pol[i]);
    cudaFree(d_pol_state[i]);
    cudaFree(d_spectrum[i]);
  }
  cudaFree(h_spectrum);

  cudaFree(d_nt);
  cudaFree(d_ns);
  cudaFree(d_dt);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);
  cudaFree(d_obs);
  cudaFree(d_se);
  cudaFree(d_c0);
  cudaFree(d_c2);
  cudaFree(d_c);
  cudaFree(d_ev2omega);
  cudaFree(h_ifirst);
  cudaFree(h_pol);
  cudaFree(h_pol_state);
  cudaFree(d_spectrum);

  // Delete host gpu pointer arrays
  delete [] event_spectrumcopy;


  return;
}








__global__ void OSCARSSR_Cuda_PowerDensityGPUMultiWithA (double  *x, double  *y, double  *z,  // position
                                                         double *bx, double *by, double *bz,  // beta
                                                         double *ax, double *ay, double *az,  // a / c
                                                         double *sx, double *sy, double *sz,  // surface coordinates
                                                         double *nx, double *ny, double *nz,  // surface normal vectors
                                                         double *dt,                          // DeltaT
                                                         int *nt,                             // number of trajectory points
                                                         int *ns,                             // number of surface elements
                                                         int *shn,                            // use normal
                                                         int *ifirst,                         // first index for this thread
                                                         double *power_density)
{
  // Check that this is within the number of spectrum points requested
  int const ith = threadIdx.x + blockIdx.x * blockDim.x;
  int const is = ith + *ifirst;
  if (is >= *ns) {
    return;
  }


  double const ox = sx[is];
  double const oy = sy[is];
  double const oz = sz[is];

  // Normal vector from input
  double const NormalX = nx[is];
  double const NormalY = ny[is];
  double const NormalZ = nz[is];

  double Sum = 0;

  for (int i = 0; i < *nt; ++i) {

    // Normal vector in direction of observation point
    double const R1 = sqrt( pow(ox - x[i], 2) + pow(oy - y[i], 2) + pow(oz - z[i], 2) );
    double const N1X = (ox - x[i]) / R1;
    double const N1Y = (oy - y[i]) / R1;
    double const N1Z = (oz - z[i]) / R1;

    // Surface normal dot with vector normal
    double const N1DotNormal = *shn == 1 ? N1X * NormalX + N1Y * NormalY + N1Z * NormalZ : 1;

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

    double const x2 = y1 * az[i] - z1 * ay[i];
    double const y2 = z1 * ax[i] - x1 * az[i];
    double const z2 = x1 * ay[i] - y1 * ax[i];

    // Numerator = N1.Cross( ( (N1 - B).Cross((AoverC)) ) );
    double const x3 = N1Y * z2 - N1Z * y2;
    double const y3 = N1Z * x2 - N1X * z2;
    double const z3 = N1X * y2 - N1Y * x2;

    double const BdotN1 = bx[i] * N1X + by[i] * N1Y + bz[i] * N1Z;
    double const Denominator = pow(1. - BdotN1, 5);

    Sum += pow( x3 * N2X + y3 * N2Y + z3 * N2Z, 2) / Denominator / (R1 * R1) * N1DotNormal;
    Sum += pow( x3 * N3X + y3 * N3Y + z3 * N3Z, 2) / Denominator / (R1 * R1) * N1DotNormal;
  }

  power_density[ith] = Sum * (*dt);

  return;
}








__global__ void OSCARSSR_Cuda_PowerDensityGPU (double  *x, double  *y, double  *z,
                                               double *bx, double *by, double *bz,
                                               double *aocx, double *aocy, double *aocz,
                                               double *sx, double *sy, double *sz,
                                               double *snx, double *sny, double *snz,
                                               double *dt,
                                               int *nt,
                                               int *ns,
                                               int *shn,
                                               double *power_density)
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
    double const N1DotNormal = *shn == 1 ? N1X * NormalX + N1Y * NormalY + N1Z * NormalZ : 1;

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



extern "C" void OSCARSSR_Cuda_CalculatePowerDensityGPUWithA (OSCARSSR& OSR,
                                                             TSurfacePoints const& Surface,
                                                             T3DScalarContainer& PowerDensityContainer,
                                                             int const NParticles,
                                                             std::vector<int> const& GPUVector)
{
  // Calculate the pd for NParticles using the GPUs given in GPUVector.  Each particle's
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


  int *h_nt, *h_nt_max, *h_ns;
  double *h_dt;
  cudaHostAlloc((void**) &h_nt_max, sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nt,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ns,     sizeof(int),    cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_dt,     sizeof(double), cudaHostAllocWriteCombined | cudaHostAllocMapped);

  // First one, set particle and trajectory
  if (!ThisParticleOnly) {
    OSR.SetNewParticle();
  }
  if (OSR.GetTrajectory().GetNPoints() == 0) {
    OSR.CalculateTrajectory();
  }

  // Needed number of points in the track and time step
  *h_nt_max = (int) OSR.GetNPointsTrajectory();
  *h_nt     = (int) OSR.GetTrajectory().GetNPoints();
  *h_ns     = (int) Surface.GetNPoints();
  *h_dt     = (double) OSR.GetTrajectory().GetDeltaT();


  int const NThreads = *h_ns;
  int const NThreadsPerBlock = NTHREADS_PER_BLOCK;
  int const NThreadsRemainder = NThreads % NThreadsPerBlock;
  int const NBlocksTotal = (NThreads - 1) / NThreadsPerBlock + 1;
  int const NBlocksPerGPU = NBlocksTotal / NGPUsToUse;
  int const NRemainderBlocks = NBlocksTotal % NGPUsToUse;
  // UPDATE: To be modified
  int const NPowerDensity = NThreadsPerBlock * (NBlocksPerGPU + (NRemainderBlocks > 0 ? 1 : 0));

  std::vector<int> NBlocksThisGPU(NGPUsToUse, NBlocksPerGPU);
  for (int i = 0; i < NRemainderBlocks; ++i) {
    ++NBlocksThisGPU[i];
  }

  // Memory allocation for Host
  double  *h_x,  *h_y,  *h_z,  *h_bx,  *h_by,  *h_bz,  *h_ax, *h_ay, *h_az, *h_sx,  *h_sy,  *h_sz, *h_nx, *h_ny, *h_nz;
  int     *h_ifirst, *h_shn;
  double **h_pd;
  cudaHostAlloc((void**) &h_x,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_y,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_z,       *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bx,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_by,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_bz,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ax,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ay,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_az,      *h_nt_max * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sx,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sy,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_sz,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nx,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ny,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_nz,          *h_ns * sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_ifirst, NGPUsToUse * sizeof(int),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void**) &h_shn,                 sizeof(int),     cudaHostAllocWriteCombined | cudaHostAllocMapped);

  cudaHostAlloc((void**) &h_pd,     NGPUsToUse * sizeof(double*), cudaHostAllocWriteCombined | cudaHostAllocMapped);
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    cudaHostAlloc((void**) &(h_pd[i]), NPowerDensity* sizeof(double),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  }

  // First surface point for each gpu
  int NBlocksUsed = 0;
  for (int i = 0; i < NGPUsToUse; ++i) {
    h_ifirst[i] = NBlocksUsed * NThreadsPerBlock;
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
  double **d_ax;
  double **d_ay;
  double **d_az;
  double **d_sx;
  double **d_sy;
  double **d_sz;
  double **d_nx;
  double **d_ny;
  double **d_nz;
  int    **d_ifirst;
  int    **d_shn;
  double **d_pd;

  cudaHostAlloc((void **) &d_nt,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ns,     NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_dt,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_x,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_y,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_z,      NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_by,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_bz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ax,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ay,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_az,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sy,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_sz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_nx,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ny,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_nz,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_ifirst, NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_shn,    NGPUsToUse * sizeof(int*),     cudaHostAllocWriteCombined | cudaHostAllocMapped);
  cudaHostAlloc((void **) &d_pd,     NGPUsToUse * sizeof(double*),  cudaHostAllocWriteCombined | cudaHostAllocMapped);

  for (size_t i = 0; i < GPUsToUse.size(); ++i) {
    // Device number
    int const d = GPUsToUse[i];

    cudaSetDevice(d);
    cudaMalloc((void **) &d_nt[i],                 sizeof(int));
    cudaMalloc((void **) &d_ns[i],                 sizeof(int));
    cudaMalloc((void **) &d_dt[i],                 sizeof(double));
    cudaMalloc((void **) &d_x[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_y[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_z[i],      *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bx[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_by[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_bz[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_ax[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_ay[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_az[i],     *h_nt_max * sizeof(double));
    cudaMalloc((void **) &d_sx[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sy[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_sz[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_nx[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_ny[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_nz[i],         *h_ns * sizeof(double));
    cudaMalloc((void **) &d_ifirst[i],             sizeof(int));
    cudaMalloc((void **) &d_shn[i],                sizeof(int));
    cudaMalloc((void **) &d_pd[i], NPowerDensity * sizeof(double));

    // Copy device number to device
    cudaMemcpyAsync(d_ifirst[i], &(h_ifirst[i]), sizeof(int), cudaMemcpyHostToDevice);
  }

  // Compute known host values
  *h_shn   = Surface.HasNormal() ? 1 : 0;
  for (size_t i = 0; i < *h_ns; ++i) {
    h_sx[i] = Surface.GetPoint(i).GetX();
    h_sy[i] = Surface.GetPoint(i).GetY();
    h_sz[i] = Surface.GetPoint(i).GetZ();

    h_nx[i] = Surface.GetPoint(i).GetNormalX();
    h_ny[i] = Surface.GetPoint(i).GetNormalY();
    h_nz[i] = Surface.GetPoint(i).GetNormalZ();
  }

  // Copy constants to first device (async)
  int const d0 = GPUsToUse[0];
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_ns[0],    h_ns,          sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_dt[0],    h_dt,          sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_shn[0],   h_shn,         sizeof(int),    cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sx[0],    h_sx,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sy[0],    h_sy,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_sz[0],    h_sz,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_nx[0],    h_nx,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_ny[0],    h_ny,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpyAsync(d_nz[0],    h_nz,  *h_ns * sizeof(double), cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_ns[i+1],     d1, d_ns[i],     d, sizeof(int));
    cudaMemcpyPeerAsync( d_dt[i+1],     d1, d_dt[i],     d, sizeof(double));
    cudaMemcpyPeerAsync( d_shn[i+1],    d1, d_shn[i],    d, sizeof(int));
    cudaMemcpyPeerAsync( d_sx[i+1],     d1, d_sx[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sy[i+1],     d1, d_sy[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_sz[i+1],     d1, d_sz[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_nx[i+1],     d1, d_nx[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_ny[i+1],     d1, d_ny[i],     d, *h_ns * sizeof(double));
    cudaMemcpyPeerAsync( d_nz[i+1],     d1, d_nz[i],     d, *h_ns * sizeof(double));
  }

  // Set first trajectory
  TParticleTrajectoryPoints const& T = OSR.GetTrajectory();
  int const NPointsThisTrajectory = T.GetNPoints();
  *h_nt = 0;
  for (size_t i = 0; i < NPointsThisTrajectory; ++i) {
    h_x[*h_nt]  = T.GetX(i).GetX();
    h_y[*h_nt]  = T.GetX(i).GetY();
    h_z[*h_nt]  = T.GetX(i).GetZ();
    h_bx[*h_nt] = T.GetB(i).GetX();
    h_by[*h_nt] = T.GetB(i).GetY();
    h_bz[*h_nt] = T.GetB(i).GetZ();
    h_ax[*h_nt] = T.GetAoverC(i).GetX();
    h_ay[*h_nt] = T.GetAoverC(i).GetY();
    h_az[*h_nt] = T.GetAoverC(i).GetZ();
    ++(*h_nt);
  }
  cudaSetDevice(d0);
  cudaMemcpyAsync(d_nt[0],    h_nt,          sizeof(int),    cudaMemcpyHostToDevice);
  for (size_t i = 0; i < GPUsToUse.size() - 1; ++i) {
    // Device number
    int const d  = GPUsToUse[i];
    int const d1 = GPUsToUse[i+1];
    cudaSetDevice(d);
    cudaMemcpyPeerAsync( d_nt[i+1],     d1, d_nt[i],     d, sizeof(int));
  }


  // Set the surface points
  // GPU events
  cudaEvent_t *event_pdcopy = new cudaEvent_t[NGPUsToUse];
  for (int ig = 0; ig < NGPUsToUse; ++ig) {
    int const d = GPUsToUse[ig];
    cudaSetDevice(d);
    cudaEventCreate(&(event_pdcopy[ig]));
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
    cudaMemcpyAsync(d_nt[0], h_nt,         sizeof(int),    cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_x[0],  h_x,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_y[0],  h_y,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_z[0],  h_z,  *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bx[0], h_bx, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_by[0], h_by, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_bz[0], h_bz, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ax[0], h_ax, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_ay[0], h_ay, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpyAsync(d_az[0], h_az, *h_nt * sizeof(double), cudaMemcpyHostToDevice);
    for (size_t ig = 0; ig < GPUsToUse.size() - 1; ++ig) {
      // Device number
      int const d  = GPUsToUse[ig];
      int const d1 = GPUsToUse[ig+1];
      cudaSetDevice(d);
      cudaMemcpyPeerAsync(d_nt[ig+1], d1, d_nt[ig], d,         sizeof(int));
      cudaMemcpyPeerAsync(d_x[ig+1],  d1, d_x[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_y[ig+1],  d1, d_y[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_z[ig+1],  d1, d_z[ig],  d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bx[ig+1], d1, d_bx[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_by[ig+1], d1, d_by[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_bz[ig+1], d1, d_bz[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ax[ig+1], d1, d_ax[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_ay[ig+1], d1, d_ay[ig], d, *h_nt * sizeof(double));
      cudaMemcpyPeerAsync(d_az[ig+1], d1, d_az[ig], d, *h_nt * sizeof(double));

    }

    // Wait for previous copy, start next one
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaEventSynchronize(event_pdcopy[ig]);
      OSCARSSR_Cuda_PowerDensityGPUMultiWithA<<<NBlocksThisGPU[ig], NThreadsPerBlock>>>( d_x[ig],  d_y[ig],  d_z[ig],
                                                                           d_bx[ig], d_by[ig], d_bz[ig],
                                                                           d_ax[ig], d_ay[ig], d_az[ig],
                                                                           d_sx[ig], d_sy[ig], d_sz[ig],
                                                                           d_nx[ig], d_ny[ig], d_nz[ig],
                                                                           d_dt[ig],
                                                                           d_nt[ig],
                                                                           d_ns[ig],
                                                                           d_shn[ig],
                                                                           d_ifirst[ig],
                                                                           d_pd[ig]);
    }


    // Add result to pd container (from **previous**)
    if (ip > 0) {
      int NBlocksUsed = 0;
      for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
        for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
          if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
            break;
          }
          int iss = ith + NThreadsPerBlock * NBlocksUsed;
          PowerDensityContainer.AddToPoint(iss, h_pd[ig][ith] * fabs(OSR.GetCurrentParticle().GetQ() * OSR.GetCurrentParticle().GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) / 1e6);
        }
        NBlocksUsed += NBlocksThisGPU[ig];
      }
    }

    // Add copy back to streams
    for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
      int const d  = GPUsToUse[ig];
      cudaSetDevice(d);
      cudaMemcpyAsync(h_pd[ig],  d_pd[ig],  NPowerDensity * sizeof(double), cudaMemcpyDeviceToHost);
      cudaEventRecord(event_pdcopy[ig]);
    }

    // If it's not the last one, calculate a new trajectory
    if (ip < NParticlesReally - 1) {
      OSR.SetNewParticle();
      OSR.CalculateTrajectory();
      TParticleTrajectoryPoints const& T = OSR.GetTrajectory();
      int const NPointsThisTrajectory = T.GetNPoints();

      *h_nt = 0;
      for (size_t it = 0; it < NPointsThisTrajectory; ++it) {
        h_x[*h_nt]  = T.GetX(it).GetX();
        h_y[*h_nt]  = T.GetX(it).GetY();
        h_z[*h_nt]  = T.GetX(it).GetZ();
        h_bx[*h_nt] = T.GetB(it).GetX();
        h_by[*h_nt] = T.GetB(it).GetY();
        h_bz[*h_nt] = T.GetB(it).GetZ();
        h_ax[*h_nt] = T.GetAoverC(it).GetX();
        h_ay[*h_nt] = T.GetAoverC(it).GetY();
        h_az[*h_nt] = T.GetAoverC(it).GetZ();
        ++(*h_nt);
      }
    }
  }

  // Wait for last copy
  for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
    cudaEventSynchronize(event_pdcopy[ig]);
  }

  // Add result to pd container (from **previous**)
  NBlocksUsed = 0;
  for (size_t ig = 0; ig < GPUsToUse.size(); ++ig) {
    for (size_t ith = 0; ith < NBlocksThisGPU[ig] * NThreadsPerBlock; ++ith) {
      if (ith + NThreadsPerBlock * NBlocksUsed >= *h_ns) {
        break;
      }
      int iss = ith + NThreadsPerBlock * NBlocksUsed;
      PowerDensityContainer.AddToPoint(iss, h_pd[ig][ith] * fabs(OSR.GetCurrentParticle().GetQ() * OSR.GetCurrentParticle().GetCurrent()) / (16 * TOSCARSSR::Pi2() * TOSCARSSR::Epsilon0() * TOSCARSSR::C()) / 1e6);
    }
    NBlocksUsed += NBlocksThisGPU[ig];
  }

  // Weighting for multi-particle
  double const Weight = 1.0 / (double) NParticlesReally;
  PowerDensityContainer.WeightAll(Weight);

  // Free host memory
  cudaFreeHost(h_nt_max);
  cudaFreeHost(h_nt);
  cudaFreeHost(h_ns);
  cudaFreeHost(h_dt);
  cudaFreeHost(h_x);
  cudaFreeHost(h_y);
  cudaFreeHost(h_z);
  cudaFreeHost(h_bx);
  cudaFreeHost(h_by);
  cudaFreeHost(h_bz);
  cudaFreeHost(h_ax);
  cudaFreeHost(h_ay);
  cudaFreeHost(h_az);
  cudaFreeHost(h_sx);
  cudaFreeHost(h_sy);
  cudaFreeHost(h_sz);
  cudaFreeHost(h_nx);
  cudaFreeHost(h_ny);
  cudaFreeHost(h_nz);
  cudaFreeHost(h_shn);
  cudaFreeHost(h_ifirst);
  // Free host and GPU memory
  for (size_t i = 0; i < GPUsToUse.size(); ++i) {

    cudaFreeHost(h_pd[i]);

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
    cudaFree(d_ax[i]);
    cudaFree(d_ay[i]);
    cudaFree(d_az[i]);
    cudaFree(d_sx[i]);
    cudaFree(d_sy[i]);
    cudaFree(d_sz[i]);
    cudaFree(d_nx[i]);
    cudaFree(d_ny[i]);
    cudaFree(d_nz[i]);
    cudaFree(d_shn[i]);
    cudaFree(d_ifirst[i]);
    cudaFree(d_pd[i]);
  }
  cudaFree(h_pd);

  cudaFree(d_nt);
  cudaFree(d_ns);
  cudaFree(d_dt);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_bx);
  cudaFree(d_by);
  cudaFree(d_bz);
  cudaFree(d_ax);
  cudaFree(d_ay);
  cudaFree(d_az);
  cudaFree(d_sx);
  cudaFree(d_sy);
  cudaFree(d_sz);
  cudaFree(d_nx);
  cudaFree(d_ny);
  cudaFree(d_nz);
  cudaFree(d_shn);
  cudaFree(h_ifirst);
  cudaFree(d_pd);

  // Delete host gpu pointer arrays
  delete [] event_pdcopy;


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

  int     shn    = Surface.HasNormal() ? 1 : 0;
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
  int    *d_nt, *d_ns, *d_shn;

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
  cudaMalloc((void **) &d_shn, sizeof(int));


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
  cudaMemcpy(d_shn, &shn, sizeof(int), cudaMemcpyHostToDevice);


  // Send computation to gpu
  int const NBlocks = NSPoints / NTHREADS_PER_BLOCK + 1;
  OSCARSSR_Cuda_PowerDensityGPU<<<NBlocks, NTHREADS_PER_BLOCK>>>(d_x, d_y, d_z, d_bx, d_by, d_bz, d_aocx, d_aocy, d_aocz, d_sx, d_sy, d_sz, d_snx, d_sny, d_snz, d_dt, d_nt, d_ns, d_shn, d_power_density);

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





