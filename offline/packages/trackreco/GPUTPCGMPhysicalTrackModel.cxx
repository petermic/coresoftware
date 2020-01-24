// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMPhysicalTrackModel.cxx
/// \author Sergey Gorbunov, David Rohr

#include "GPUTPCGMPhysicalTrackModel.h"

int GPUTPCGMPhysicalTrackModel::PropagateToXBzLight(float x, float Bz, float& dLp)
{
  GPUTPCGMPhysicalTrackModel t = *this;
  if (fabs(x - t.X()) < 1.e-8f) {
    return 0;
  }
  int err = t.PropagateToXBzLightNoUpdate(x, Bz, dLp);
  if (err) {
    return (err);
  }
  t.UpdateValues();
  *this = t;
  return 0;
}

int GPUTPCGMPhysicalTrackModel::PropagateToXBzLightNoUpdate(float x, float Bz, float& dLp)
{
  //
  // transport the track to X=x in magnetic field B = ( 0, 0, Bz[kG*0.000299792458] )
  // dLp is a return value == path length / track momentum [cm/(GeV/c)]
  // the method returns error code (0 == no error)
  //
  // Additional values are not recalculated, UpdateValues() has to be called afterwards!!
  //
  float b = mQ * Bz;
  float pt2 = mPx * mPx + mPy * mPy;
  float dx = x - mX;
  float pye = mPy - dx * b; // extrapolated py
  float pxe2 = pt2 - pye * pye;

  if (mPx < (1.f - GPUCA_MAX_SIN_PHI) || pxe2 < (1.f - GPUCA_MAX_SIN_PHI) * (1.f - GPUCA_MAX_SIN_PHI)) {
    return -1; // can not transport to x=x
  }
  float pxe = sqrt(pxe2); // extrapolated px
  float pti = 1.f / sqrt(pt2);

  float ty = (mPy + pye) / (mPx + pxe);
  float dy = dx * ty;
  float dS; // path in XY
  {
    float chord = dx * sqrt(1.f + ty * ty); // chord to the extrapolated point == sqrt(dx^2+dy^2)*sign(dx)
    float sa = 0.5f * chord * b * pti;              //  sin( half of the rotation angle ) ==  (chord/2) / radius

    // dS = (Pt/b)*2*arcsin( sa )
    //    = (Pt/b)*2*sa*(1 + 1/6 sa^2 + 3/40 sa^4 + 5/112 sa^6 +... )
    //    =       chord*(1 + 1/6 sa^2 + 3/40 sa^4 + 5/112 sa^6 +... )

    float sa2 = sa * sa;
    const float k2 = 1.f / 6.f;
    const float k4 = 3.f / 40.f;
    // const float k6 = 5.f/112.f;
    dS = chord + chord * sa2 * (k2 + k4 * sa2);
    // dS = sqrt(pt2)/b*2.*asin( sa );
  }

  dLp = pti * dS; // path in XYZ / p == path in XY / pt

  float dz = mPz * dLp;

  mX = x;
  mY += dy;
  mZ += dz;
  mPx = pxe;
  mPy = pye;
  // mPz = mPz;
  // mQ = mQ;
  return 0;
}

int GPUTPCGMPhysicalTrackModel::PropagateToXBxByBz(float x, float Bx, float By, float Bz, float& dLp)
{
  //
  // transport the track to X=x in magnetic field B = ( Bx, By, Bz )[kG*0.000299792458]
  // xyzPxPyPz as well as all the additional values will change. No need to call UpdateValues() afterwards.
  // the method returns error code (0 == no error)
  //

  if (0) { // simple transport in Bz for test proposes
    return PropagateToXBzLight(x, Bz, dLp);
  }
  dLp = 0.f;

  GPUTPCGMPhysicalTrackModel t = *this;

  // Rotate to the system where Bx=By=0.

  float bt = sqrt(Bz * Bz + By * By);
  float bb = sqrt(Bx * Bx + By * By + Bz * Bz);

  float c1 = 1.f, s1 = 0.f;
  float c2 = 1.f, s2 = 0.f;

  if (bt > 1.e-4f) {
    c1 = Bz / bt;
    s1 = By / bt;
    c2 = bt / bb;
    s2 = -Bx / bb;
  }

  // rotation matrix: first around x, then around y'
  // after the first rotation: Bx'==Bx, By'==0, Bz'==Bt, X'==X
  // after the second rotation: Bx''==0, By''==0, Bz''==B, X'' axis is as close as possible to the original X

  //
  //     ( c2 0 s2 )   ( 1  0   0 )
  // R = (  0 1 0  ) X ( 0 c1 -s1 )
  //     (-s2 0 c2 )   ( 0 s1  c1 )
  //

  float R0[3] = {c2, s1 * s2, c1 * s2};
  float R1[3] = {0, c1, -s1};
  float R2[3] = {-s2, s1 * c2, c1 * c2};

  // parameters and the extrapolation point in the rotated coordinate system
  {
    float lx = t.X(), ly = t.Y(), lz = t.Z(), lpx = t.Px(), lpy = t.Py(), lpz = t.Pz();

    t.X() = R0[0] * lx + R0[1] * ly + R0[2] * lz;
    t.Y() = R1[0] * lx + R1[1] * ly + R1[2] * lz;
    t.Z() = R2[0] * lx + R2[1] * ly + R2[2] * lz;

    t.Px() = R0[0] * lpx + R0[1] * lpy + R0[2] * lpz;
    t.Py() = R1[0] * lpx + R1[1] * lpy + R1[2] * lpz;
    t.Pz() = R2[0] * lpx + R2[1] * lpy + R2[2] * lpz;
  }

  float dx = x - mX;
  float xe = t.X() + dx; // propagate on same dx in rotated system

  // transport in rotated coordinate system to X''=xe:

  if (t.Px() < (1.f - GPUCA_MAX_SIN_PHI)) {
    t.Px() = 1.f - GPUCA_MAX_SIN_PHI;
  }
  if (t.PropagateToXBzLightNoUpdate(xe, bb, dLp) != 0) {
    return -1;
  }

  // rotate coordinate system back to the original R{-1}==R{T}
  {
    float lx = t.X(), ly = t.Y(), lz = t.Z(), lpx = t.Px(), lpy = t.Py(), lpz = t.Pz();

    t.X() = R0[0] * lx + R1[0] * ly + R2[0] * lz;
    t.Y() = R0[1] * lx + R1[1] * ly + R2[1] * lz;
    t.Z() = R0[2] * lx + R1[2] * ly + R2[2] * lz;

    t.Px() = R0[0] * lpx + R1[0] * lpy + R2[0] * lpz;
    t.Py() = R0[1] * lpx + R1[1] * lpy + R2[1] * lpz;
    t.Pz() = R0[2] * lpx + R1[2] * lpy + R2[2] * lpz;
  }

  // a small (hopefully) additional step to X=x. Perhaps it may be replaced by linear extrapolation.

  float ddLp = 0;
  if (t.Px() < (1.f - GPUCA_MAX_SIN_PHI)) {
    t.Px() = 1.f - GPUCA_MAX_SIN_PHI;
  }
  if (t.PropagateToXBzLightNoUpdate(x, Bz, ddLp) != 0) {
    return -1;
  }

  dLp += ddLp;

  t.UpdateValues();
  *this = t;
  return 0;
}

int GPUTPCGMPhysicalTrackModel::PropagateToLpBz(float Lp, float Bz)
{
  // Lp is path length L over track momentum p in [cm/GeV], Bz in kG*clight
  //
  // it is a copy of AliExternalTrackParam: ghelix3 routine.
  //
  // the method returns error code (0 == no error)
  //

  float qfield = mQ * Bz;

  float step = Lp;

  const float kOvSqSix = sqrt(1.f / 6.f);

  float px = mPx;
  float py = mPy;
  float pz = mPz;

  float tet = qfield * step;

  float tsint, sintt, sint, cos1t;
  if (fabs(tet) > 0.03f) {
    sint = sin(tet);
    sintt = sint / tet;
    tsint = (tet - sint) / tet;
    float t = sin(0.5f * tet);
    cos1t = 2.f * t * t / tet;
  } else {
    tsint = tet * tet / 6.f;
    sintt = (1.f - tet * kOvSqSix) * (1.f + tet * kOvSqSix); // 1.- tsint;
    sint = tet * sintt;
    cos1t = 0.5f * tet;
  }

  float f1 = step * sintt;
  float f2 = step * cos1t;
  float f3 = step * tsint;
  float f4 = -tet * cos1t;
  float f5 = sint;

  mX += f1 * px - f2 * py;
  mY += f1 * py + f2 * px;
  mZ += f1 * pz + f3 * pz;

  mPx += f4 * px - f5 * py;
  mPy += f4 * py + f5 * px;

  UpdateValues();

  return 0;
}

#include <iostream>

void GPUTPCGMPhysicalTrackModel::Print() const
{
  std::cout << "GPUTPCGMPhysicalTrackModel:  x " << mX << " y " << mY << " z " << mZ << " px " << mPx << " py " << mPy << " pz " << mPz << " q " << mQ << std::endl;
}

GPUTPCGMPhysicalTrackModel::GPUTPCGMPhysicalTrackModel(const GPUTPCGMTrackParam& t) { Set(t); }

void GPUTPCGMPhysicalTrackModel::Set(const GPUTPCGMTrackParam& t)
{
  float pti = fabs(t.GetQPt());
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
  if (fabs(px) < 1.e-4f) {
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
  if (fabs(Bz) < 1.e-8f) {
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
