// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCGMOfflineStatisticalErrors.h
/// \author David Rohr

#ifndef GPUTPCGMOFFLINESTATISTICALERRORS
#define GPUTPCGMOFFLINESTATISTICALERRORS

#include "GPUTPCGMMergedTrackHit.h"

struct GPUTPCGMMergedTrackHit;

struct GPUTPCGMOfflineStatisticalErrors {
  GPUd() void SetCurCluster(GPUTPCGMMergedTrackHit* /*c*/) {}
  GPUd() void GetOfflineStatisticalErrors(float& /*err2Y*/, float& /*err2Z*/, float /*sinPhi*/, float /*dzds*/, unsigned char /*clusterState*/) const {}
};

#endif
