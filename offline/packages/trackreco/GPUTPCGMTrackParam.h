// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMTrackParam.h
/// \author David Rohr, Sergey Gorbunov

#ifndef GPUTPCGMTRACKPARAM_H
#define GPUTPCGMTRACKPARAM_H

#include "GPUTPCDef.h"
#include "GPUTPCGMMergedTrackHit.h"
#include "GPUdEdxInfo.h"

#include <cmath>
#include <algorithm>

class AliExternalTrackParam;

class GPUTPCGMMerger;
class GPUTPCGMBorderTrack;
struct GPUParam;
class GPUTPCGMPhysicalTrackModel;
class GPUTPCGMPolynomialField;
class GPUTPCGMMergedTrack;
class GPUTPCGMPropagator;

/**
 * @class GPUTPCGMTrackParam
 *
 * GPUTPCGMTrackParam class describes the track parametrisation
 * which is used by the GPUTPCGMTracker slice tracker.
 *
 */
class GPUTPCGMTrackParam
{
 public:
  struct GPUTPCOuterParam {
    float X, alpha;
    float P[5];
    float C[15];
  };

  float& X()
  {
    return mX;
  }
  float& Y()
  {
    return mP[0];
  }
  float& Z()
  {
    return mP[1];
  }
  float& SinPhi()
  {
    return mP[2];
  }
  float& DzDs()
  {
    return mP[3];
  }
  float& QPt()
  {
    return mP[4];
  }
  float& ZOffset()
  {
    return mZOffset;
  }

  float GetX() const { return mX; }
  float GetY() const { return mP[0]; }
  float GetZ() const { return mP[1]; }
  float GetSinPhi() const { return mP[2]; }
  float GetDzDs() const { return mP[3]; }
  float GetQPt() const { return mP[4]; }
  float GetZOffset() const { return mZOffset; }

  float GetKappa(float Bz) const { return -mP[4] * Bz; }

  void SetX(float v) { mX = v; }

  float* Par()
  {
    return mP;
  }
  const float* GetPar() const { return mP; }
  float GetPar(int i) const { return (mP[i]); }
  void SetPar(int i, float v) { mP[i] = v; }

  float& Chi2()
  {
    return mChi2;
  }
  int& NDF()
  {
    return mNDF;
  }

  float Err2Y() const { return mC[0]; }
  float Err2Z() const { return mC[2]; }
  float Err2SinPhi() const { return mC[5]; }
  float Err2DzDs() const { return mC[9]; }
  float Err2QPt() const { return mC[14]; }

  float GetChi2() const { return mChi2; }
  int GetNDF() const { return mNDF; }

  float GetCosPhi() const { return sqrt(float(1.f) - GetSinPhi() * GetSinPhi()); }

  float GetErr2Y() const { return mC[0]; }
  float GetErr2Z() const { return mC[2]; }
  float GetErr2SinPhi() const { return mC[5]; }
  float GetErr2DzDs() const { return mC[9]; }
  float GetErr2QPt() const { return mC[14]; }

  float* Cov()
  {
    return mC;
  }

  const float* GetCov() const { return mC; }
  float GetCov(int i) const { return mC[i]; }

  void SetCov(int i, float v) { mC[i] = v; }
  void SetChi2(float v) { mChi2 = v; }
  void SetNDF(int v) { mNDF = v; }

  float GetMirroredY(float Bz) const;

  void ResetCovariance();

  bool CheckNumericalQuality(float overrideCovYY = -1.f) const;
  bool CheckCov() const;

  bool Fit(const GPUTPCGMMerger* merger, int iTrk, GPUTPCGMMergedTrackHit* clusters, int& N, int& NTolerated, float& Alpha, int attempt = 0, float maxSinPhi = GPUCA_MAX_SIN_PHI, GPUTPCOuterParam* outerParam = nullptr, GPUdEdxInfo* dEdxOut = nullptr);
  void MirrorTo(GPUTPCGMPropagator& prop, float toY, float toZ, bool inFlyDirection, const GPUParam& param, unsigned char row, unsigned char clusterState, bool mirrorParameters);
  int MergeDoubleRowClusters(int ihit, int wayDirection, GPUTPCGMMergedTrackHit* clusters, const GPUTPCGMMerger* merger, GPUTPCGMPropagator& prop, float& xx, float& yy, float& zz, int maxN, float clAlpha, unsigned char& clusterState, bool rejectChi2);

  void AttachClustersMirror(const GPUTPCGMMerger* Merger, int slice, int iRow, int iTrack, float toY, GPUTPCGMPropagator& prop);
  void AttachClustersPropagate(const GPUTPCGMMerger* Merger, int slice, int lastRow, int toRow, int iTrack, bool goodLeg, GPUTPCGMPropagator& prop, bool inFlyDirection, float maxSinPhi = GPUCA_MAX_SIN_PHI);
  void AttachClusters(const GPUTPCGMMerger* Merger, int slice, int iRow, int iTrack, bool goodLeg);
  void AttachClusters(const GPUTPCGMMerger* Merger, int slice, int iRow, int iTrack, bool goodLeg, float Y, float Z);

  int FollowCircle(const GPUTPCGMMerger* Merger, GPUTPCGMPropagator& prop, int slice, int iRow, int iTrack, bool goodLeg, float toAlpha, float toX, float toY, int toSlice, int toRow, bool inFlyDirection);

  void MarkClusters(GPUTPCGMMergedTrackHit* clusters, int ihitFirst, int ihitLast, int wayDirection, unsigned char state)
  {
    clusters[ihitFirst].state |= state;
    while (ihitFirst != ihitLast) {
      ihitFirst += wayDirection;
      clusters[ihitFirst].state |= state;
    }
  }
  void UnmarkClusters(GPUTPCGMMergedTrackHit* clusters, int ihitFirst, int ihitLast, int wayDirection, unsigned char state)
  {
    clusters[ihitFirst].state &= ~state;
    while (ihitFirst != ihitLast) {
      ihitFirst += wayDirection;
      clusters[ihitFirst].state &= ~state;
    }
  }

  bool Rotate(float alpha);
  void ShiftZ(const GPUTPCGMMergedTrackHit* clusters, const GPUTPCGMMerger* merger, int N);

  static float Reciprocal(float x) { return 1.f / x; }
  static void Assign(float& x, bool mask, float v)
  {
    if (mask) {
      x = v;
    }
  }

  static void Assign(int& x, bool mask, int v)
  {
    if (mask) {
      x = v;
    }
  }

  static void RefitTrack(GPUTPCGMMergedTrack& track, int iTrk, const GPUTPCGMMerger* merger, GPUTPCGMMergedTrackHit* clusters);

#if defined(GPUCA_ALIROOT_LIB) & !defined(GPUCA_GPUCODE)
  bool GetExtParam(AliExternalTrackParam& T, double alpha) const;
  void SetExtParam(const AliExternalTrackParam& T);
#endif

  void ConstrainSinPhi(float limit = GPUCA_MAX_SIN_PHI)
  {
    if (mP[2] > limit) {
      mP[2] = limit;
    } else if (mP[2] < -limit) {
      mP[2] = -limit;
    }
  }

 private:
  bool FollowCircleChk(float lrFactor, float toY, float toX, bool up, bool right);
  int initResetT0();

  float mX; // x position
  float mZOffset;
  float mP[5];  // 'active' track parameters: Y, Z, SinPhi, DzDs, q/Pt
  float mC[15]; // the covariance matrix for Y,Z,SinPhi,..
  float mChi2;  // the chi^2 value
  int mNDF;     // the Number of Degrees of Freedom
};

int GPUTPCGMTrackParam::initResetT0()
{
  const float absQPt = abs(mP[4]);
  if (absQPt < (150.f / 40.f)) {
    return 150.f / 40.f;
  }
  return max(10.f, 150.f / mP[4]);
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
  if (abs(qptBz) < 1.e-8f) {
    qptBz = 1.e-8f;
  }
  if (cosPhi2 < 0.f) {
    cosPhi2 = 0.f;
  }
  return GetY() - 2.f * sqrt(cosPhi2) / qptBz;
}

#endif
