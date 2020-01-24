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
  GPUTPCGMPhysicalTrackModel() = default;
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


#endif
