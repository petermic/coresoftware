// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMPropagator.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMPROPAGATOR_H
#define GPUTPCGMPROPAGATOR_H

#include "GPUTPCGMOfflineStatisticalErrors.h"
#include "GPUTPCGMPhysicalTrackModel.h"
#include "GPUTPCGMPolynomialField.h"

namespace o2
{
namespace base
{
class MatLayerCylSet;
} // namespace base
} // namespace o2

class GPUTPCGMTrackParam;
struct GPUParam;

/**
 * @class GPUTPCGMPropagator
 *
 */

class GPUTPCGMPropagator
{
 public:
  /// Enumeration of field regions
  enum FieldRegion {
    TPC = 0, ///< TPC
    ITS = 1, ///< ITS
    TRD = 2  ///< outer TPC -> outer TRD
  };

  GPUTPCGMPropagator() = default;

  struct MaterialCorrection {
    MaterialCorrection() : radLen(29.532f), rho(1.025e-3f), rhoOverRadLen(rho / radLen), DLMax(0.f), EP2(0.f), sigmadE2(0.f), k22(0.f), k33(0.f), k43(0.f), k44(0.f) {}

    float radLen, rho, rhoOverRadLen, DLMax, EP2, sigmadE2, k22, k33, k43, k44; // precalculated values for MS and EnergyLoss correction
  };

  void SetMaterial(float radLen, float rho);

  void SetPolynomialField(const GPUTPCGMPolynomialField* field) { mField = field; }

  void SelectFieldRegion(FieldRegion region) { mFieldRegion = region; }

  void SetFitInProjections(bool Flag) { mFitInProjections = Flag; }
  void SetToyMCEventsFlag(bool Flag) { mToyMCEvents = Flag; }
  void SetSeedingErrors(bool Flag) { mSeedingErrors = Flag; }
  void SetMatLUT(const o2::base::MatLayerCylSet* lut) { mMatLUT = lut; }

  void SetMaxSinPhi(float maxSinPhi) { mMaxSinPhi = maxSinPhi; }

  void SetTrack(GPUTPCGMTrackParam* track, float Alpha);
  void ResetT0()
  {
    if (!mT) {
      return;
    }
    mT0.Set(*mT);
  }

  int RotateToAlpha(float newAlpha);

  int PropagateToXAlpha(float posX, float posAlpha, bool inFlyDirection);

  int PropagateToXAlphaBz(float posX, float posAlpha, bool inFlyDirection);

  int Update(float posY, float posZ, int iRow, const GPUParam& param, short clusterState, bool rejectChi2, bool refit);
  int Update(float posY, float posZ, short clusterState, bool rejectChi2, float err2Y, float err2Z);
  float PredictChi2(float posY, float posZ, int iRow, const GPUParam& param, short clusterState) const;
  float PredictChi2(float posY, float posZ, float err2Y, float err2Z) const;
  int RejectCluster(float chiY, float chiZ, unsigned char clusterState)
  {
    if (chiY > 9.f || chiZ > 9.f) {
      return 2;
    }
    if ((chiY > 6.25f || chiZ > 6.25f) && (clusterState & (GPUTPCGMMergedTrackHit::flagSplit | GPUTPCGMMergedTrackHit::flagShared))) {
      return 2;
    }
    if ((chiY > 1.f || chiZ > 6.25f) && (clusterState & (GPUTPCGMMergedTrackHit::flagEdge | GPUTPCGMMergedTrackHit::flagSingle))) {
      return 2;
    }
    return 0;
  }

  float GetBz(float Alpha, float X, float Y, float Z) const;
  void GetBxByBz(float Alpha, float X, float Y, float Z, float B[3]) const;

  void GetErr2(float& err2Y, float& err2Z, const GPUParam& param, float posZ, int iRow, short clusterState) const;

  float GetAlpha() const { return mAlpha; }
  float GetQPt0() const { return mT0.GetQPt(); }
  float GetSinPhi0() const { return mT0.GetSinPhi(); }
  float GetCosPhi0() const { return mT0.GetCosPhi(); }
  void Mirror(bool inFlyDirection);
  void Rotate180();
  void ChangeDirection();
  float GetMirroredYModel() const;
  float GetMirroredYTrack() const;
  int GetPropagatedYZ(float x, float& projY, float& projZ);
  bool GetFitInProjections() const { return mFitInProjections; }

  GPUTPCGMPhysicalTrackModel& Model()
  {
    return mT0;
  }
  void CalculateMaterialCorrection();
  void SetStatErrorCurCluster(GPUTPCGMMergedTrackHit* c) { mStatErrors.SetCurCluster(c); }

 private:
  static float ApproximateBetheBloch(float beta2);
  int FollowLinearization(const GPUTPCGMPhysicalTrackModel& t0e, float Bz, float dLp, bool inFlyDirection);

  const GPUTPCGMPolynomialField* mField = nullptr;
  FieldRegion mFieldRegion = TPC;

  GPUTPCGMTrackParam* mT = nullptr;
  float mAlpha = 0; // rotation angle of the track coordinate system
  GPUTPCGMPhysicalTrackModel mT0;
  MaterialCorrection mMaterial;
  bool mSeedingErrors = 0;
  bool mFitInProjections = 1; // fit (Y,SinPhi,QPt) and (Z,DzDs) paramteres separatelly
  bool mToyMCEvents = 0;      // events are simulated with simple home-made simulation
  float mMaxSinPhi = GPUCA_MAX_SIN_PHI;

  GPUTPCGMOfflineStatisticalErrors mStatErrors;
  const o2::base::MatLayerCylSet* mMatLUT = nullptr;
};


#endif
