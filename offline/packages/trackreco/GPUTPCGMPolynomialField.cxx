// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMPolynomialField.cxx
/// \author Sergey Gorbunov, David Rohr

#include "GPUTPCGMPolynomialField.h"

#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

void GPUTPCGMPolynomialField::Print() const
{
  const double kCLight = 0.000299792458;
  typedef std::numeric_limits<float> flt;
  cout << std::scientific;
#if __cplusplus >= 201103L
  cout << std::setprecision(flt::max_digits10 + 2);
#endif
  cout << " nominal field " << mNominalBz << " [kG * (2.99792458E-4 GeV/c/kG/cm)]"
       << " == " << mNominalBz / kCLight << " [kG]" << endl;

  cout << " TpcBx[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBx[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TpcBy[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBy[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TpcBz[NTPCM] = { ";
  for (int i = 0; i < NTPCM; i++) {
    cout << mTpcBz[i];
    if (i < NTPCM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << "TRD field: \n"
       << endl;

  cout << " TrdBx[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBx[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TrdBy[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBy[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " TrdBz[NTRDM] = { ";
  for (int i = 0; i < NTRDM; i++) {
    cout << mTrdBz[i];
    if (i < NTRDM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << "ITS field: \n"
       << endl;

  cout << " ItsBx[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBx[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " ItsBy[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBy[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }

  cout << " ItsBz[NITSM] = { ";
  for (int i = 0; i < NITSM; i++) {
    cout << mItsBz[i];
    if (i < NITSM - 1) {
      cout << ", ";
    } else {
      cout << " };" << endl;
    }
  }
}

void GPUTPCGMPolynomialField::GetPolynomsTpc(float x, float y, float z, float f[NTPCM])
{
  f[0] = 1.f;
  f[1] = x;
  f[2] = y;
  f[3] = z;
  f[4] = x * x;
  f[5] = x * y;
  f[6] = x * z;
  f[7] = y * y;
  f[8] = y * z;
  f[9] = z * z;
}

void GPUTPCGMPolynomialField::GetField(float x, float y, float z, float B[3]) const
{
  const float* fBxS = &mTpcBx[1];
  const float* fByS = &mTpcBy[1];
  const float* fBzS = &mTpcBz[1];

  const float f[NTPCM - 1] = {x, y, z, x * x, x * y, x * z, y * y, y * z, z * z};
  float bx = mTpcBx[0], by = mTpcBy[0], bz = mTpcBz[0];
  for (int i = NTPCM - 1; i--;) {
    // for (int i=0;i<NTPCM-1; i++){
    bx += fBxS[i] * f[i];
    by += fByS[i] * f[i];
    bz += fBzS[i] * f[i];
  }
  B[0] = bx;
  B[1] = by;
  B[2] = bz;
}

float GPUTPCGMPolynomialField::GetFieldBz(float x, float y, float z) const
{
  const float* fBzS = &mTpcBz[1];

  const float f[NTPCM - 1] = {x, y, z, x * x, x * y, x * z, y * y, y * z, z * z};
  float bz = mTpcBz[0];
  for (int i = NTPCM - 1; i--;) {
    bz += fBzS[i] * f[i];
  }
  return bz;
}

void GPUTPCGMPolynomialField::GetPolynomsTrd(float x, float y, float z, float f[NTRDM])
{
  float xx = x * x, xy = x * y, xz = x * z, yy = y * y, yz = y * z, zz = z * z;
  f[0] = 1.f;
  f[1] = x;
  f[2] = y;
  f[3] = z;
  f[4] = xx;
  f[5] = xy;
  f[6] = xz;
  f[7] = yy;
  f[8] = yz;
  f[9] = zz;
  f[10] = x * xx;
  f[11] = x * xy;
  f[12] = x * xz;
  f[13] = x * yy;
  f[14] = x * yz;
  f[15] = x * zz;
  f[16] = y * yy;
  f[17] = y * yz;
  f[18] = y * zz;
  f[19] = z * zz;
}

void GPUTPCGMPolynomialField::GetFieldTrd(float x, float y, float z, float B[3]) const
{
  float f[NTRDM];
  GetPolynomsTrd(x, y, z, f);
  float bx = 0.f, by = 0.f, bz = 0.f;
  for (int i = 0; i < NTRDM; i++) {
    bx += mTrdBx[i] * f[i];
    by += mTrdBy[i] * f[i];
    bz += mTrdBz[i] * f[i];
  }
  B[0] = bx;
  B[1] = by;
  B[2] = bz;
}

float GPUTPCGMPolynomialField::GetFieldTrdBz(float x, float y, float z) const
{
  float f[NTRDM];
  GetPolynomsTrd(x, y, z, f);
  float bz = 0.f;
  for (int i = 0; i < NTRDM; i++) {
    bz += mTrdBz[i] * f[i];
  }
  return bz;
}

void GPUTPCGMPolynomialField::GetPolynomsIts(float x, float y, float z, float f[NITSM])
{
  float xx = x * x, xy = x * y, xz = x * z, yy = y * y, yz = y * z, zz = z * z;
  f[0] = 1.f;
  f[1] = x;
  f[2] = y;
  f[3] = z;
  f[4] = xx;
  f[5] = xy;
  f[6] = xz;
  f[7] = yy;
  f[8] = yz;
  f[9] = zz;
  /*
        f[10]=x*xx; f[11]=x*xy; f[12]=x*xz; f[13]=x*yy; f[14]=x*yz; f[15]=x*zz;
        f[16]=y*yy; f[17]=y*yz; f[18]=y*zz;
        f[19]=z*zz;
   */
}

void GPUTPCGMPolynomialField::GetFieldIts(float x, float y, float z, float B[3]) const
{
  const float* fBxS = &mItsBx[1];
  const float* fByS = &mItsBy[1];
  const float* fBzS = &mItsBz[1];

  const float f[NITSM - 1] = {x, y, z, x * x, x * y, x * z, y * y, y * z, z * z};
  float bx = mItsBx[0], by = mItsBy[0], bz = mItsBz[0];
  for (int i = NITSM - 1; i--;) {
    bx += fBxS[i] * f[i];
    by += fByS[i] * f[i];
    bz += fBzS[i] * f[i];
  }
  B[0] = bx;
  B[1] = by;
  B[2] = bz;
}

float GPUTPCGMPolynomialField::GetFieldItsBz(float x, float y, float z) const
{
  const float* fBzS = &mItsBz[1];

  const float f[NITSM - 1] = {x, y, z, x * x, x * y, x * z, y * y, y * z, z * z};
  float bz = mItsBz[0];
  for (int i = NITSM - 1; i--;) {
    bz += fBzS[i] * f[i];
  }
  return bz;
}
