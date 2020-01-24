// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMTrackParam.cxx
/// \author David Rohr, Sergey Gorbunov

#define GPUCA_CADEBUG 0
#define DEBUG_SINGLE_TRACK -1

#define CADEBUG(expr) expr

#include "GPUTPCDef.h"
#include "GPUTPCGMTrackParam.h"
#include "GPUTPCGMPhysicalTrackModel.h"
#include "GPUTPCGMPropagator.h"
#include "GPUTPCGMPolynomialField.h"
#include "GPUParam.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>

static constexpr float kRho = 1.025e-3f;  // 0.9e-3;
static constexpr float kRadLen = 29.532f; // 28.94;
static constexpr float kDeg2Rad = M_PI / 180.f;
static constexpr float kSectAngle = 2 * M_PI / 18.f;

void GPUTPCGMTrackParam::MirrorTo(GPUTPCGMPropagator& prop, float toY, float toZ, bool inFlyDirection, const GPUParam& param, unsigned char row, unsigned char clusterState, bool mirrorParameters)
{
  if (mirrorParameters) {
    prop.Mirror(inFlyDirection);
  }
  float err2Y, err2Z;
  prop.GetErr2(err2Y, err2Z, param, toZ, row, clusterState);
  prop.Model().Y() = mP[0] = toY;
  prop.Model().Z() = mP[1] = toZ;
  if (mC[0] < err2Y) {
    mC[0] = err2Y;
  }
  if (mC[2] < err2Z) {
    mC[2] = err2Z;
  }
  if (fabs(mC[5]) < 0.1f) {
    mC[5] = mC[5] > 0 ? 0.1f : -0.1f;
  }
  if (mC[9] < 1.f) {
    mC[9] = 1.f;
  }
  mC[1] = mC[4] = mC[6] = mC[8] = mC[11] = mC[13] = 0;
  prop.SetTrack(this, prop.GetAlpha());
  mNDF = -3;
  mChi2 = 0;
}

bool GPUTPCGMTrackParam::FollowCircleChk(float lrFactor, float toY, float toX, bool up, bool right)
{
  return fabs(mX * lrFactor - toY) > 1.f &&                                                                       // transport further in Y
         fabs(mP[2]) < 0.7f &&                                                                                    // rotate back
         (up ? (-mP[0] * lrFactor > toX || (right ^ (mP[2] > 0))) : (-mP[0] * lrFactor < toX || (right ^ (mP[2] < 0)))); // don't overshoot in X
}

bool GPUTPCGMTrackParam::CheckCov() const
{
  const float* c = mC;
  bool ok = c[0] >= 0 && c[2] >= 0 && c[5] >= 0 && c[9] >= 0 && c[14] >= 0 && (c[1] * c[1] <= c[2] * c[0]) && (c[3] * c[3] <= c[5] * c[0]) && (c[4] * c[4] <= c[5] * c[2]) && (c[6] * c[6] <= c[9] * c[0]) && (c[7] * c[7] <= c[9] * c[2]) && (c[8] * c[8] <= c[9] * c[5]) &&
            (c[10] * c[10] <= c[14] * c[0]) && (c[11] * c[11] <= c[14] * c[2]) && (c[12] * c[12] <= c[14] * c[5]) && (c[13] * c[13] <= c[14] * c[9]);
  return ok;
}

bool GPUTPCGMTrackParam::CheckNumericalQuality(float overrideCovYY) const
{
  //* Check that the track parameters and covariance matrix are reasonable
  bool ok = std::isfinite(mX) && std::isfinite(mChi2);
  CADEBUG(
    printf("OK %d - ", (int)ok); for (int i = 0; i < 5; i++) { printf("%f ", mP[i]); } printf(" - "); for (int i = 0; i < 15; i++) { printf("%f ", mC[i]); } printf("\n"));
  const float* c = mC;
  for (int i = 0; i < 15; i++) {
    ok = ok && std::isfinite(c[i]);
  }
  CADEBUG(printf("OK1 %d\n", (int)ok));
  for (int i = 0; i < 5; i++) {
    ok = ok && std::isfinite(mP[i]);
  }
  CADEBUG(printf("OK2 %d\n", (int)ok));
  if ((overrideCovYY > 0 ? overrideCovYY : c[0]) > 4.f * 4.f || c[2] > 4.f * 4.f || c[5] > 2.f * 2.f || c[9] > 2.f * 2.f) {
    ok = 0;
  }
  CADEBUG(printf("OK3 %d\n", (int)ok));
  if (fabs(mP[2]) > GPUCA_MAX_SIN_PHI) {
    ok = 0;
  }
  CADEBUG(printf("OK4 %d\n", (int)ok));
  if (!CheckCov()) {
    ok = false;
  }
  CADEBUG(printf("OK5 %d\n", (int)ok));
  return ok;
}

#if defined(GPUCA_ALIROOT_LIB) & !defined(GPUCA_GPUCODE)
bool GPUTPCGMTrackParam::GetExtParam(AliExternalTrackParam& T, double alpha) const
{
  //* Convert from GPUTPCGMTrackParam to AliExternalTrackParam parameterisation,
  //* the angle alpha is the global angle of the local X axis

  bool ok = CheckNumericalQuality();

  double par[5], cov[15];
  for (int i = 0; i < 5; i++) {
    par[i] = mP[i];
  }
  for (int i = 0; i < 15; i++) {
    cov[i] = mC[i];
  }

  if (par[2] > GPUCA_MAX_SIN_PHI) {
    par[2] = GPUCA_MAX_SIN_PHI;
  }
  if (par[2] < -GPUCA_MAX_SIN_PHI) {
    par[2] = -GPUCA_MAX_SIN_PHI;
  }

  if (fabs(par[4]) < 1.e-5) {
    par[4] = 1.e-5; // some other software will crash if q/Pt==0
  }
  if (fabs(par[4]) > 1. / 0.08) {
    ok = 0; // some other software will crash if q/Pt is too big
  }
  T.Set((double)mX, alpha, par, cov);
  return ok;
}

void GPUTPCGMTrackParam::SetExtParam(const AliExternalTrackParam& T)
{
  //* Convert from AliExternalTrackParam parameterisation

  for (int i = 0; i < 5; i++) {
    mP[i] = T.GetParameter()[i];
  }
  for (int i = 0; i < 15; i++) {
    mC[i] = T.GetCovariance()[i];
  }
  mX = T.GetX();
  if (mP[2] > GPUCA_MAX_SIN_PHI) {
    mP[2] = GPUCA_MAX_SIN_PHI;
  }
  if (mP[2] < -GPUCA_MAX_SIN_PHI) {
    mP[2] = -GPUCA_MAX_SIN_PHI;
  }
}
#endif

bool GPUTPCGMTrackParam::Rotate(float alpha)
{
  float cA = cos(alpha);
  float sA = sin(alpha);
  float x0 = mX;
  float sinPhi0 = mP[2], cosPhi0 = sqrt(1 - mP[2] * mP[2]);
  float cosPhi = cosPhi0 * cA + sinPhi0 * sA;
  float sinPhi = -cosPhi0 * sA + sinPhi0 * cA;
  float j0 = cosPhi0 / cosPhi;
  float j2 = cosPhi / cosPhi0;
  mX = x0 * cA + mP[0] * sA;
  mP[0] = -x0 * sA + mP[0] * cA;
  mP[2] = sinPhi + j2;
  mC[0] *= j0 * j0;
  mC[1] *= j0;
  mC[3] *= j0;
  mC[6] *= j0;
  mC[10] *= j0;

  mC[3] *= j2;
  mC[4] *= j2;
  mC[5] *= j2 * j2;
  mC[8] *= j2;
  mC[12] *= j2;
  if (cosPhi < 0) { // change direction ( t0 direction is already changed in t0.UpdateValues(); )
    SinPhi() = -SinPhi();
    DzDs() = -DzDs();
    QPt() = -QPt();
    mC[3] = -mC[3];
    mC[4] = -mC[4];
    mC[6] = -mC[6];
    mC[7] = -mC[7];
    mC[10] = -mC[10];
    mC[11] = -mC[11];
  }
  return true;
}

int GPUTPCGMTrackParam::initResetT0()
{
  const float absQPt = fabs(mP[4]);
  if (absQPt < (150.f / 40.f)) {
    return 150.f / 40.f;
  }
  return std::max(10.f, 150.f / mP[4]);
}

void GPUTPCGMTrackParam::ResetCovariance()
{
  mC[0] = 100.f;
  mC[1] = 0.f;
  mC[2] = 100.f;
  mC[3] = 0.f;
  mC[4] = 0.f;
  mC[5] = 1.f;
  mC[6] = 0.f;
  mC[7] = 0.f;
  mC[8] = 0.f;
  mC[9] = 10.f;
  mC[10] = 0.f;
  mC[11] = 0.f;
  mC[12] = 0.f;
  mC[13] = 0.f;
  mC[14] = 10.f;
  mChi2 = 0;
  mNDF = -5;
}

float GPUTPCGMTrackParam::GetMirroredY(float Bz) const
{
  // get Y of the point which has the same X, but located on the other side of trajectory
  float qptBz = GetQPt() * Bz;
  float cosPhi2 = 1.f - GetSinPhi() * GetSinPhi();
  if (fabs(qptBz) < 1.e-8f) {
    qptBz = 1.e-8f;
  }
  if (cosPhi2 < 0.f) {
    cosPhi2 = 0.f;
  }
  return GetY() - 2.f * sqrt(cosPhi2) / qptBz;
}
