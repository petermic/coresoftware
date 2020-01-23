// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMPhysicalTrackModel.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMPHYSICALTRACKMODEL_H
#define GPUTPCGMPHYSICALTRACKMODEL_H

#include "GPUTPCGMTrackParam.h"

/**
 * @class GPUTPCGMPhysicalTrackModel
 *
 * GPUTPCGMPhysicalTrackModel class is a trajectory in physical parameterisation (X,Y,Z,Px,PY,Pz,Q)
 * without covariance matrix. Px>0 and Q is {-1,+1} (no uncharged tracks).
 *
 * It is used to linearise transport equations for GPUTPCGMTrackParam trajectory during (re)fit.
 *
 */

class GPUTPCGMPhysicalTrackModel
{
 public:
  GPUTPCGMPhysicalTrackModel();
  GPUTPCGMPhysicalTrackModel(const GPUTPCGMTrackParam& t);

  void Set(const GPUTPCGMTrackParam& t);
  void Set(float X, float Y, float Z, float Px, float Py, float Pz, float Q);

  float& X()
  {
    return mX;
  }
  float& Y()
  {
    return mY;
  }
  float& Z()
  {
    return mZ;
  }
  float& Px()
  {
    return mPx;
  }
  float& Py()
  {
    return mPy;
  }
  float& Pz()
  {
    return mPz;
  }
  float& Q()
  {
    return mQ;
  }

  float& SinPhi()
  {
    return mSinPhi;
  }
  float& CosPhi()
  {
    return mCosPhi;
  }
  float& SecPhi()
  {
    return mSecPhi;
  }
  float& DzDs()
  {
    return mDzDs;
  }
  float& DlDs()
  {
    return mDlDs;
  }
  float& QPt()
  {
    return mQPt;
  }
  float& P()
  {
    return mP;
  }
  float& Pt()
  {
    return mPt;
  }

  const float& SinPhi() const { return mSinPhi; }
  const float& DzDs() const { return mDzDs; }

  float GetX() const { return mX; }
  float GetY() const { return mY; }
  float GetZ() const { return mZ; }
  float GetPx() const { return mPx; }
  float GetPy() const { return mPy; }
  float GetPz() const { return mPz; }
  float GetQ() const { return mQ; }

  float GetSinPhi() const { return mSinPhi; }
  float GetCosPhi() const { return mCosPhi; }
  float GetSecPhi() const { return mSecPhi; }
  float GetDzDs() const { return mDzDs; }
  float GetDlDs() const { return mDlDs; }
  float GetQPt() const { return mQPt; }
  float GetP() const { return mP; }
  float GetPt() const { return mPt; }

  int PropagateToXBzLightNoUpdate(float x, float Bz, float& dLp);
  int PropagateToXBzLight(float x, float Bz, float& dLp);

  int PropagateToXBxByBz(float x, float Bx, float By, float Bz, float& dLp);

  int PropagateToLpBz(float Lp, float Bz);

  bool SetDirectionAlongX();

  void UpdateValues();

  void Print() const;

  float GetMirroredY(float Bz) const;

  void Rotate(float alpha);
  void RotateLight(float alpha);

 private:
  // physical parameters of the trajectory

  float mX = 0.f;    // X
  float mY = 0.f;    // Y
  float mZ = 0.f;    // Z
  float mPx = 1.e4f; // Px, >0
  float mPy = 0.f;   // Py
  float mPz = 0.f;   // Pz
  float mQ = 1.f;    // charge, +-1

  // some additional variables needed for GMTrackParam transport

  float mSinPhi = 0.f; // SinPhi = Py/Pt
  float mCosPhi = 1.f; // CosPhi = abs(Px)/Pt
  float mSecPhi = 1.f; // 1/cos(phi) = Pt/abs(Px)
  float mDzDs = 0.f;   // DzDs = Pz/Pt
  float mDlDs = 0.f;   // DlDs = P/Pt
  float mQPt = 0.f;    // QPt = q/Pt
  float mP = 1.e4f;    // momentum
  float mPt = 1.e4f;   // Pt momentum
};

GPUTPCGMPhysicalTrackModel::GPUTPCGMPhysicalTrackModel(const GPUTPCGMTrackParam& t) { Set(t); }

void GPUTPCGMPhysicalTrackModel::Set(const GPUTPCGMTrackParam& t)
{
  float pti = abs(t.GetQPt());
  if (pti < 1.e-4f) {
    pti = 1.e-4f; // set 10000 GeV momentum for straight track
  }
  mQ = (t.GetQPt() >= 0) ? 1.f : -1.f; // only charged tracks are considered
  mX = t.GetX();
  mY = t.GetY();
  mZ = t.GetZ();

  mPt = 1.f / pti;
  mSinPhi = t.GetSinPhi();
  if (mSinPhi > GPUCA_MAX_SIN_PHI) {
    mSinPhi = GPUCA_MAX_SIN_PHI;
  }
  if (mSinPhi < -GPUCA_MAX_SIN_PHI) {
    mSinPhi = -GPUCA_MAX_SIN_PHI;
  }
  mCosPhi = sqrt((1.f - mSinPhi) * (1.f + mSinPhi));
  mSecPhi = 1.f / mCosPhi;
  mDzDs = t.GetDzDs();
  mDlDs = sqrt(1.f + mDzDs * mDzDs);
  mP = mPt * mDlDs;

  mPy = mPt * mSinPhi;
  mPx = mPt * mCosPhi;
  mPz = mPt * mDzDs;
  mQPt = mQ * pti;
}

void GPUTPCGMPhysicalTrackModel::Set(float X, float Y, float Z, float Px, float Py, float Pz, float Q)
{
  mX = X;
  mY = Y;
  mZ = Z;
  mPx = Px;
  mPy = Py;
  mPz = Pz;
  mQ = (Q >= 0) ? 1 : -1;
  UpdateValues();
}

void GPUTPCGMPhysicalTrackModel::UpdateValues()
{
  float px = mPx;
  if (abs(px) < 1.e-4f) {
    px = copysign(1.e-4f, px);
  }

  mPt = sqrt(px * px + mPy * mPy);
  float pti = 1.f / mPt;
  mP = sqrt(px * px + mPy * mPy + mPz * mPz);
  mSinPhi = mPy * pti;
  mCosPhi = px * pti;
  mSecPhi = mPt / px;
  mDzDs = mPz * pti;
  mDlDs = mP * pti;
  mQPt = mQ * pti;
}

bool GPUTPCGMPhysicalTrackModel::SetDirectionAlongX()
{
  //
  // set direction of movenment collinear to X axis
  // return value is true when direction has been changed
  //
  if (mPx >= 0) {
    return 0;
  }

  mPx = -mPx;
  mPy = -mPy;
  mPz = -mPz;
  mQ = -mQ;
  UpdateValues();
  return 1;
}

float GPUTPCGMPhysicalTrackModel::GetMirroredY(float Bz) const
{
  // get Y of the point which has the same X, but located on the other side of trajectory
  if (abs(Bz) < 1.e-8f) {
    Bz = 1.e-8f;
  }
  return mY - 2.f * mQ * mPx / Bz;
}

void GPUTPCGMPhysicalTrackModel::RotateLight(float alpha)
{
  //* Rotate the coordinate system in XY on the angle alpha

  float cA = cos(alpha);
  float sA = sin(alpha);
  float x = mX, y = mY, px = mPx, py = mPy;
  mX = x * cA + y * sA;
  mY = -x * sA + y * cA;
  mPx = px * cA + py * sA;
  mPy = -px * sA + py * cA;
}

void GPUTPCGMPhysicalTrackModel::Rotate(float alpha)
{
  //* Rotate the coordinate system in XY on the angle alpha
  RotateLight(alpha);
  UpdateValues();
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
