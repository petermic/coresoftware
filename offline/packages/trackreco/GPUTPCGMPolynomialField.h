// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMPolynomialField.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMPOLYNOMIALFIELD_H
#define GPUTPCGMPOLYNOMIALFIELD_H

#include "GPUTPCDef.h"

/**
 * @class GPUTPCGMPolynomialField
 *
 */

class GPUTPCGMPolynomialField
{
 public:
  GPUTPCGMPolynomialField() : mNominalBz(0.f)
  {
    Reset();
  }

  inline void Reset()
  {
    mNominalBz = 0.f;
    for (int i = 0; i < NTPCM; i++) {
      mTpcBx[i] = 0.f;
      mTpcBy[i] = 0.f;
      mTpcBz[i] = 0.f;
    }
    for (int i = 0; i < NTRDM; i++) {
      mTrdBx[i] = 0.f;
      mTrdBy[i] = 0.f;
      mTrdBz[i] = 0.f;
    }
    for (int i = 0; i < NITSM; i++) {
      mItsBx[i] = 0.f;
      mItsBy[i] = 0.f;
      mItsBz[i] = 0.f;
    }
  }


  inline void SetFieldNominal(float nominalBz) { mNominalBz = nominalBz; }

  inline void SetFieldTpc(const float* Bx, const float* By, const float* Bz)
  {
    if (Bx && By && Bz) {
      for (int i = 0; i < NTPCM; i++) {
        mTpcBx[i] = Bx[i];
        mTpcBy[i] = By[i];
        mTpcBz[i] = Bz[i];
      }
    }
  }

  inline void SetFieldTrd(const float* Bx, const float* By, const float* Bz)
  {
    if (Bx && By && Bz) {
      for (int i = 0; i < NTRDM; i++) {
        mTrdBx[i] = Bx[i];
        mTrdBy[i] = By[i];
        mTrdBz[i] = Bz[i];
      }
    }
  }

  inline void SetFieldIts(const float* Bx, const float* By, const float* Bz)
  {
    if (Bx && By && Bz) {
      for (int i = 0; i < NITSM; i++) {
        mItsBx[i] = Bx[i];
        mItsBy[i] = By[i];
        mItsBz[i] = Bz[i];
      }
    }
  }

  float GetNominalBz() const { return mNominalBz; }

  void GetField(float x, float y, float z, float B[3]) const;
  float GetFieldBz(float x, float y, float z) const;

  void GetFieldTrd(float x, float y, float z, float B[3]) const;
  float GetFieldTrdBz(float x, float y, float z) const;

  void GetFieldIts(float x, float y, float z, float B[3]) const;
  float GetFieldItsBz(float x, float y, float z) const;

  void Print() const;

  static constexpr int NTPCM = 10; // number of coefficients
  static constexpr int NTRDM = 20; // number of coefficients for the TRD field
  static constexpr int NITSM = 10; // number of coefficients for the ITS field

  static void GetPolynomsTpc(float x, float y, float z, float f[NTPCM]);
  static void GetPolynomsTrd(float x, float y, float z, float f[NTRDM]);
  static void GetPolynomsIts(float x, float y, float z, float f[NITSM]);

  const float* GetCoefmTpcBx() const { return mTpcBx; }
  const float* GetCoefmTpcBy() const { return mTpcBy; }
  const float* GetCoefmTpcBz() const { return mTpcBz; }

  const float* GetCoefmTrdBx() const { return mTrdBx; }
  const float* GetCoefmTrdBy() const { return mTrdBy; }
  const float* GetCoefmTrdBz() const { return mTrdBz; }

  const float* GetCoefmItsBx() const { return mItsBx; }
  const float* GetCoefmItsBy() const { return mItsBy; }
  const float* GetCoefmItsBz() const { return mItsBz; }

 private:
  float mNominalBz;    // nominal constant field value in [kG * 2.99792458E-4 GeV/c/cm]
  float mTpcBx[NTPCM]; // polynomial coefficients
  float mTpcBy[NTPCM];
  float mTpcBz[NTPCM];
  float mTrdBx[NTRDM]; // polynomial coefficients
  float mTrdBy[NTRDM];
  float mTrdBz[NTRDM];
  float mItsBx[NITSM]; // polynomial coefficients
  float mItsBy[NITSM];
  float mItsBz[NITSM];
};


#endif
