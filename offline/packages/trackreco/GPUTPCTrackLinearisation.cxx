#include "GPUTPCTrackLinearisation.h"

GPUTPCTrackLinearisation::GPUTPCTrackLinearisation(const GPUTPCTrackParam& t) : mSinPhi(t.SinPhi()), mCosPhi(0), mDzDs(t.DzDs()), mQPt(t.QPt())
{
  if (mSinPhi > GPUCA_MAX_SIN_PHI) {
    mSinPhi = GPUCA_MAX_SIN_PHI;
  } else if (mSinPhi < -GPUCA_MAX_SIN_PHI) {
    mSinPhi = -GPUCA_MAX_SIN_PHI;
  }
  mCosPhi = sqrt(1 - mSinPhi * mSinPhi);
  if (t.SignCosPhi() < 0) {
    mCosPhi = -mCosPhi;
  }
}

void GPUTPCTrackLinearisation::Set(float SinPhi1, float CosPhi1, float DzDs1, float QPt1)
{
  SetSinPhi(SinPhi1);
  SetCosPhi(CosPhi1);
  SetDzDs(DzDs1);
  SetQPt(QPt1);
}
